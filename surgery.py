import os
from warnings import warn
import pandas as pd
import scanpy as sc
import sparse

import scarches as sca
import numpy as np
from scipy import sparse
import gdown
import anndata as ad
import torch




def run_surgery_process(adata: ad.AnnData, name: str, batch_key: str = 'sample', 
                        ref_model_dir_prefix: str = '/root/datos/maestria/netopaas/lung_scRNA/HCA_Lung/',
                        surgery_epochs: int = 350,
                        surgery_save=False,
                        uncertainty_threshold = 0.2, # This uncertainty threshold limits the false positive rate to <0.5 (as per Sikkema et al., bioRxiv 2022)
                        n_neighbors = 30, 
                        early_stopping_kwargs_surgery: dict = None) -> ad.AnnData:
    """
    Process the given AnnData object using a surgical method based on provided parameters.

    Parameters:
    - adata (ad.AnnData): The input data object containing scRNA-seq expression values.
    - name (str): A name or label associated with the dataset.
    - batch_key (str, optional): The key in adata object to be used for batch effect. Defaults to 'sample'.
    - ref_model_dir_prefix (str, optional): Directory prefix for the reference model and embeddings. 
        Defaults to '/root/datos/maestria/netopaas/lung_scRNA/HCA_Lung/'.
    - surgery_epochs (int, optional): Number of epochs for the surgical process. Defaults to 350.
    - surgery_save (bool, optional): Flag to determine if the surgery model should be saved. Defaults to False.
    - uncertainty_threshold (float, optional): Threshold for determining uncertainty during surgery process. 
        Defaults to 0.2.
    - n_neighbors (int, optional): Number of neighbors for uncertainty determination. Defaults to 30.
        Note: As per Sikkema et al., bioRxiv 2022, this uncertainty threshold limits the false positive rate to <0.5.
    - early_stopping_kwargs_surgery (dict, optional): Dictionary containing early stopping parameters specific 
        for the surgical process. If None, default values will be applied. Defaults to None.

    Returns:
    - ad.AnnData: Processed AnnData with the atals combiened

    Notes:
    - Requires GPU for optimal performance. A warning will be issued if GPU is not detected.
    - This function is a part of a surgical procedure on scRNA-seq data and may require additional
      context and parameters for a comprehensive understanding and execution.

    Examples:
    ```python
    processed_adata = run_surgery_process(adata, "sample_name", batch_key="sample_id")
    ```

    Raises:
    - Warning: If GPU is not detected.

    """
    
    if not torch.cuda.is_available():
        warn("You are not using GPU stop before it gets too far")
    
    path_reference_emb = os.path.join(ref_model_dir_prefix, "HLCA_emb_and_metadata.h5ad")
    ref_model_dir = os.path.join(ref_model_dir_prefix, "HLCA_reference_model")
    surgery_model_dir_prefix = ref_model_dir_prefix
    surgery_model_dir = os.path.join(surgery_model_dir_prefix, "surgery_model")

    if not early_stopping_kwargs_surgery:
        early_stopping_kwargs_surgery = {
            "early_stopping_monitor": "elbo_train",
            "early_stopping_patience": 10,
            "early_stopping_min_delta": 0.001,
            "plan_kwargs": {"weight_decay": 0.0},
        }

    adata_ref = sc.read_h5ad(path_reference_emb)


    adata_query_unprep = adata.copy()

    # We need RAW count data for the model we are using so we get it from adata.layers
    adata_query_unprep.X = adata_query_unprep.layers['counts']


    # For faster computation convert the matrix to sparse.
    # First check with `adata.X` if it is sparse
    if not(type(adata_query_unprep.X) in [sparse._csc.csc_matrix, sparse._csr.csr_matrix]):
        warning.warn("Your matrix is not in sparse format, so computation \
                     will take up much more RAM, might even be slower. \
                     Converting to sparse to spare RAM")
        adata_query_unprep.X = sparse.csr_matrix(adata_query_unprep.X)


    # #### Change to Ensembl IDs

    # Revert to ensemble ids using the assembly reference gch38 (or hg19 in the comments) used in the experiment.
    # Doing this via THE PACKAGE Biomart always leaves us with some gaps because there are many manually annotated genes.

    # If your query feature naming does not match your reference model feature naming, you will need to add the right feature names. For the HLCA reference, the mapping of the 2000 input gene IDs to their gene names is stored on the HLCA Zenodo page, so you can add gene ids using that mapping. Alternatively, you can map your gene IDs to gene names (or the reverse) using BioMart mapping tables. In most cases your raw data includes both gene IDs and names, in which case mapping is not necessary.
    # 
    # Let’s download the HLCA-specific gene mapping:
    path_gene_mapping_df = f'{ref_model_dir_prefix}/HLCA_reference_model_gene_order_ids_and_symbols.csv'


    # Download gene information from HLCA github:
    if not os.path.exists(path_gene_mapping_df):
        url = "https://zenodo.org/record/7599104/files/HLCA_reference_model_gene_order_ids_and_symbols.csv"
        gdown.download(url, path_gene_mapping_df, quiet=True)


    gene_id_to_gene_name_df = pd.read_csv(path_gene_mapping_df, index_col=0)

    # Store your gene names in an adata.var.column if they are currently the index:
    # if gene names are in .var.index:
    adata_query_unprep.var["gene_names"] = adata_query_unprep.var.index
    gene_name_column_name = "gene_names"


    # Map gene names to gene ids for all of the 2000 reference model genes that we can find in our data:
    # 
    # Check number of detected genes:
    n_overlap = (
        adata_query_unprep.var[gene_name_column_name]
        .isin(gene_id_to_gene_name_df.gene_symbol)
        .sum()
    )
    n_genes_model = gene_id_to_gene_name_df.shape[0]
    print(
        f"Number of model input genes detected: {n_overlap} out of {n_genes_model} ({round(n_overlap/n_genes_model*100)}%)"
    )


    # Subset query data to only the genes that are part of the modeling input, then map gene names to gene ids using the table above. Store the resulting ids both in the .var.index (for scArches) and in a .var[gene_ids] (for merging duplicate genes).
    adata_query_unprep = adata_query_unprep[
        :,
        adata_query_unprep.var[gene_name_column_name].isin(
            gene_id_to_gene_name_df.gene_symbol
        ),
    ].copy()  # subset your data to genes used in the reference model
    adata_query_unprep.var.index = adata_query_unprep.var[gene_name_column_name].map(
        dict(zip(gene_id_to_gene_name_df.gene_symbol, gene_id_to_gene_name_df.index))
    )  # add gene ids for the gene names, and store in .var.index
    # remove index name to prevent bugs later on
    adata_query_unprep.var.index.name = None
    adata_query_unprep.var["gene_ids"] = adata_query_unprep.var.index


    def sum_by(adata: ad.AnnData, col: str) -> ad.AnnData:
        adata.obs[col] = adata.obs[col].astype('category')
        assert pd.api.types.is_categorical_dtype(adata.obs[col])

        cat = adata.obs[col].values
        indicator = sparse.coo_matrix(
            (np.broadcast_to(True, adata.n_obs), (cat.codes, np.arange(adata.n_obs))),
            shape=(len(cat.categories), adata.n_obs),
        )

        return ad.AnnData(
            indicator @ adata.X, var=adata.var, obs=pd.DataFrame(index=cat.categories)
        )


    # shape before merging:
    # Now merge. Note that all var columns will be dropped after merging (as we don’t specify how to merge). As the merging is done based on .obs indices in the function above, we transpose our anndata object and re-transpose it after merging.

    adata_query_unprep = sum_by(adata_query_unprep.transpose(), col="gene_ids").transpose()

    # add back gene ids:
    adata_query_unprep.var = adata_query_unprep.var.join(gene_id_to_gene_name_df).rename(columns={"gene_symbol":"gene_names"})


    # #### Surgery

    # We pad missing query genes with zeros and reorder the available ones to ensure data corectness and smooth running of the scArches reference mapping.
    adata_query = sca.models.SCANVI.prepare_query_anndata(
        adata = adata_query_unprep,
        # return_reference_var_names=True,
        reference_model = ref_model_dir,
        inplace=False)


    # This line should be kept unchanged due to the structure of the pre-trained reference model.
    adata_query.obs['scanvi_label'] = 'unlabeled'


    # Now we perform scArches “surgery”.
    # Thanks to exploring the model we can deduce what .obs we need and the name we assigned. In this case we need to define a 'dataset'
    # column in obs to define the batches.
    # We can load the model with this command or train it if we havent got it
    adata_query.obs[batch_key] = adata_query.obs[batch_key].astype('string')
    adata_query.obs['dataset'] = adata_query.obs[batch_key]

    # We convert all cetogrical types to strings to be abl to save into an h5ad file.
    # Sometimes we have to run this command after trying to save strangely.
    adata_query.var.mito = False
    for col in adata_query.obs.columns:
        if pd.api.types.is_categorical_dtype(adata_query.obs[col]):
            adata_query.obs[col] = adata_query.obs[col].astype('str')


    # TODO possible checkpoint
    # Convert all mito falgs into false because we are using the latent genes of the HCA_lung model
    # adata_query = sc.read_h5ad(f'{backup_dir}/surgeries/query_{name}.h5ad')
    # os.makedirs(f'{backup_dir}/surgeries', exist_ok=True)
    # adata_query.write_h5ad(f'{backup_dir}/surgeries/query_{name}.h5ad')
    
    
    surgery_model = sca.models.SCANVI.load_query_data(
            adata_query,
            ref_model_dir,
            freeze_dropout = True,
        )

    # Run the neural network surgery
    surgery_model.train(max_epochs=surgery_epochs,
                        **early_stopping_kwargs_surgery)


    if surgery_save:
        print('saving the model')
        surgery_model.save(f'{backup_dir}/surgeries/{name}', overwrite=True)


    # #### Get latent representation
    # Here we will calculate the “latent representation”, or “low-dimensional embedding” of your dataset. This embedding is in the same space as the HLCA core/reference embedding that you loaded in the beginning of the script. Hence, we can combine the two embeddings afterwards (HLCA + your new data), and do joint clustering, UMAP embedding, label transfer etc.!
    adata_query_latent = sc.AnnData(surgery_model.get_latent_representation(adata_query))
    adata_query_latent.obs = adata_query.obs.loc[adata_query.obs.index,:]


    # #### Combine embeddings
    # We add “reference or query” metadata to acquire more information and better analyse the integration level.
    adata_query_latent.obs['ref_or_query'] = "query"
    adata_ref.obs['ref_or_query'] = "ref"


    # We will now combine the two embeddings to enable joing clustering etc. If you expect non-unique barcodes (.obs index), set index_unique to e.g. “_” and batch_key to the obs column that you want to use as barcode suffix (e.g. “dataset”).
    combined_emb = adata_ref.concatenate(adata_query_latent, index_unique=None) # index_unique="_", batch_key="dataset") # alternative

    # #### Label transfer
    # Next, we use a knn classifier to transfer the lables from the reference to the query. We do this for every level of the annotation (i.e. level 1-5). Note that some cell types don’t have annotations for higher levels, e.g. mast cells do not have level 4 or 5 annotations. For those cell types, we “propagate” to the higher levels, i.e. you will see “3_Mast cells” in level 4 and 5 annotations. (Most cell types don’t have a level 5 annotation!) Therefore, all highest level annotations can be found under level 5.
    celltypes = f'{ref_model_dir_prefix}/HLCA_celltypes_ordered.csv'

    # url = 'https://github.com/LungCellAtlas/mapping_data_to_the_HLCA/raw/main/supporting_files/HLCA_celltypes_ordered.csv'
    # gdown.download(url, celltypes, quiet=False)
    # TODO check if exists
    cts_ordered = pd.read_csv(celltypes,index_col=0)
    
    #Now run the label transfer commands. Note that this might take quite a while if you have a large query dataset! For our small test dataset, it should not take long.
    print("Starting kkn label transfer")
    knn_transformer = sca.utils.knn.weighted_knn_trainer(
        train_adata=adata_ref,
        train_adata_emb="X",
        n_neighbors=50,
    )

    labels, uncert = sca.utils.knn.weighted_knn_transfer(query_adata=adata_query_latent, 
                                                         query_adata_emb="X", # location of our joint embedding
                                                         label_keys="Level",
                                                         knn_model=knn_transformer,
                                                         ref_adata_obs = adata_ref.obs.join(cts_ordered, on='ann_finest_level'))


    # With the commands above, we labeled every cell from the query. However, some cells might have high label transfer uncertainty. It is useful to set those to “unknown” instead of giving them a cell type label. This will help highlight cell types/states that are new (i.e. not present in the reference) and possible interesting, they’re worth taking a careful look at!
    labels.rename(columns={f'Level_{lev}':f'Level_{lev}_transfered_label_unfiltered' for lev in range(1,6)},inplace=True)
    uncert.rename(columns={f'Level_{lev}':f'Level_{lev}_transfer_uncert' for lev in range(1,6)},inplace=True)

    combined_emb.obs = combined_emb.obs.join(labels)
    combined_emb.obs = combined_emb.obs.join(uncert)

    t_labels = [f'Level_{lev}_transfered_label_unfiltered' for lev in range(1,6)]
    t_uncert = [f'Level_{lev}_transfer_uncert' for lev in range(1,6)]


    # Convert uncertainties to arrays
    combined_emb.obs[t_uncert] = list(np.array(combined_emb.obs[t_uncert]))
    # Convert cell type labels to categoricals, and set “nan” to NaN
    def remove_uncert_types(combined_emb):
        import pandas as pd
        t_labels = [f'Level_{lev}_transfered_label_unfiltered' for lev in range(1,6)]
        t_uncert = [f'Level_{lev}_transfer_uncert' for lev in range(1,6)]

        combined_emb.obs[t_uncert] = list(np.array(combined_emb.obs[t_uncert]))

        for col, uncert in zip(t_labels,t_uncert):
            filtered_colname = col.replace('_unfiltered','')
            # too high uncertainty levels => set to "Unknown"
            combined_emb.obs[filtered_colname] = combined_emb.obs[col]
            combined_emb.obs[filtered_colname] = combined_emb.obs[filtered_colname].astype('str')
            combined_emb.obs[filtered_colname].mask(
                combined_emb.obs[uncert] > uncertainty_threshold,
                'Unknown',
                inplace = True)

            # convert to categorical:
            combined_emb.obs[col] = pd.Categorical(combined_emb.obs[col])
            combined_emb.obs[filtered_colname] = pd.Categorical(combined_emb.obs[filtered_colname])
            # then replace "nan" with NaN (that makes colors better in umap)
            combined_emb.obs[col].replace('nan',np.nan,inplace=True)
            combined_emb.obs[filtered_colname].replace('nan',np.nan,inplace=True)
        return combined_emb
    
    remove_uncert_types(combined_emb)


    # Let’s take a look at the percentage of cells set to “unknown” after our filtering:
    print(f'Percentage of unknown per level, with uncertainty_threshold={uncertainty_threshold}:')
    for level in range(1,6):
        print(f"Level {level}: {np.round(sum(combined_emb.obs[f'Level_{level}_transfered_label'] =='Unknown')/adata_query.n_obs*100,2)}%")

    # ### UMAP
    # 
    # The UMAP plots help us perform downstream analysis, like clustering, label transfer, integration and more.
    # 
    # #### UMAP Query vs. Reference
    print("Starting neighbours of combined_emb")
    sc.pp.neighbors(combined_emb, n_neighbors=n_neighbors)
    print("Starting dimension reduction tu umap")
    sc.tl.umap(combined_emb)

    return combined_emb