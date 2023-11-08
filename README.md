# Comp-oncology
Example notebook for the recommended analyses for single cell cancer data.
We use the lung cancer data from [Laughney et. al](https://www.nature.com/articles/s41591-019-0750-6) for immune development
in lung cancer.

## Pre-reqs
We highly recommend using a computer or server with cuda enabled.
[Docker](https://docs.docker.com/engine/install/ubuntu/) must be a given to run the containers.
Another advantage would be to have docker installed with the cuda package
for transfering your cuda install to the containers.

Nevertheless, we added a list of the package versions for R and python in `requirements.txt`
and `r_packages.txt`. Or you can base yourself on the commands in the Dockerfile if you 
don't want to use docker.

The notebook has commands to download the `.h5ad` file for scanpy. It contains the raw count matrix
in `Anndata` form.

## Running

We have designed a docker image that has all the necessary libs. It is however very big, because it has 
all the R and python packages including the ones for ML.

To be able to run the notebook, you should run the container with this notebook repo
inserted as a volume. In the next command, `$HOME/2021-SC-HCA-LATAM/CONTAINER` is the
path of the repo, and the other one is where the large data files would be stored.

```
docker run --interactive --runtime=nvidia --gpus all --tty --name comp_onco --publish 8888-8892:8888-8892 --volume $HOME/2021-SC-HCA-LATAM/CONTAINER:/root/host_home --volume /datos:/root/datos --workdir /root/host_home/ netopaas/comp-onco:paga /bin/bash
```

We publish some ports to use the jupyter server.


After that just run the command `jl` inside the container and a jupyter lab server will be launched.

## Running in servers

It is recommended to run the notebook on a server with lots of RAM.
 Another advantage of having a container is that you can run everything without giving sudo access. 
If you are doing so, you can launch a reverse proxy to access this server from 
your computer with a command like this:

```bash
ssh -p 5265 -N -f -L 5432:localhost:8888 sefirot.inmegen.gob.mx
```

This forwards the port 8888 in the server to 5432 in your localhost. The -p is for when the ssh port on your server is different from the default.


## Troubleshooting

If you keep your repo inside the container, the user might change. We recommend configuring the ssh keys inside the container. These are the steps:

```
git config --global user.email ernesto.paas@ciencias.unam.mx
git config --global user.name "Ernesto Paas"
```

We suggest saving a key pair that can be generated with the command `ssh-keygen -t ed25519 -C "your_email@example.com".`
in the folder that contains the repo, and then copying them to `~/.ssh/id_ed25519` etc.
