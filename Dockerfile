FROM python:3.11

SHELL ["/bin/bash", "-c"]

# For editing in shell and pre-reqs for scArches and R
RUN apt update && apt install -y --no-install-recommends nano nodejs gfortran
 
# For cell annotation
RUN pip install scarches==0.5.9

# For notebooks
RUN pip install jupyterlab==3.2.4
RUN echo "alias jl='jupyter-lab --no-browser --ip=0.0.0.0 --allow-root /root/host_home'" >> ~/.bashrc

# For visualizing and builidng networks
RUN apt install -y --no-install-recommends libgraphviz-dev graphviz
RUN pip install pygraphviz==1.11 networkx==3.1

# For R
RUN wget https://cran.rstudio.com/src/base/R-4/R-4.3.1.tar.gz
RUN tar xvfz R-4.3.1.tar.gz && rm R-4.3.1.tar.gz
WORKDIR R-4.3.1
RUN ./configure --enable-R-shlib --with-cairo --with-libpng --prefix=/opt/R/
RUN make && make install
WORKDIR /opt/R
RUN rm -rf /opt/R/R-4.3.1
ENV PATH="/opt/R/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/R/lib/R/lib:${LD_LIBRARY_PATH}"
RUN echo 'PATH="/opt/R/bin:${PATH}"' >> ~/.bashrc
RUN echo 'LD_LIBRARY_PATH="/opt/R/lib/R/lib:${LD_LIBRARY_PATH}"' >> ~/.bashrc

RUN echo "/opt/R/lib/R/lib" >> /etc/ld.so.conf.d/RLibs.conf && ldconfig

# R libraries
RUN apt install -y --no-install-recommends libharfbuzz-dev libfribidi-dev
RUN cat ~/.bashrc

RUN Rscript -e "update.packages(ask=FALSE, repos='https://cran.itam.mx/')"
RUN Rscript -e "install.packages(c('devtools', 'gam', 'RColorBrewer', 'BiocManager', 'IRkernel','png'), repos='https://cran.itam.mx/')"
RUN Rscript -e "IRkernel::installspec(user = FALSE)"
RUN Rscript -e "BiocManager::install(c('sparseMatrixStats', 'SparseArray', 'DelayedMatrixStats','scuttle', 'scry', 'edgeR'))"
RUN Rscript -e "devtools::install_github('MatteoBlla/PsiNorm')"

# Other pip packages
RUN pip install triku==2.1.6 rpy2==3.5.14 anndata2ri==1.3.1 ikarus==0.0.3 Cython==0.29.33 tables==3.9.1 infercnvpy==0.4.3
WORKDIR ~
RUN git clone https://github.com/cvanelteren/forceatlas2
RUN cd forceatlas2 && echo "from forceatlas2 import *" >> fa2/__init__.py && pip install .

# JAVA for networks
RUN apt install default-jre


RUN apt-get clean -y && apt-get autoremove -y