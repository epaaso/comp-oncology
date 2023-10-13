# Comp-oncology
Example notebook for the recommended analyses for single cell cancer data

## Pre-reqs
We highly recommend using computer or server with cuda enabled.
AAnother adavnatage would be to have dokcer installed with the cuda package
for translating your cuda install to the contianers.

Nevertheless there is a list of the libraries and python and R packages at the end.


## Running

We have designed a docker image that has all the necessary libs. It is however very big, because it has 
all the R and python pacages including the ones for ML.

To be able to run the notebook, you should run the container with this notebook repo
inserted as a volume. In the next command, `$HOME/2021-SC-HCA-LATAM/CONTAINER` is the
path of the repo, and the other one is where the large data files would be stored.

```
docker run --interactive --runtime=nvidia --gpus all --tty --name comp_onco --publish 8888-8892:8888-8892 --volume $HOME/2021-SC-HCA-LATAM/CONTAINER:/root/host_home --volume /datos:/root/datos --workdir /root/host_home/ scarches /bin/bash
```

We publish some ports to use the jupyter server.


After that just run the command `jl` and a jupyter lab server will be launched.

## Running in servers

It is recommended to run the notebook in a server with much RAM. Another advantage of having a container is that
you can run everything witout having sudo acces.
If you are doing so you can do a reverse proxy to acces this sever from your computer with 
a command like this:

```bash
ssh -p 5265 -N -f -L 5432:localhost:8888 sefirot.inmegen.gob.mx
```

This forwards the port 8888 in the server to 5432 in your localhost. The -p is for when the ssh port on your server is different than the default.


## Troubleshooting

If you keep your repo inside the container, the user might change. We recommend configuring the ssh keys inside the container. These are the steps:

```
git config --global user.email ernesto.paas@ciencias.unam.mx
git config --global user.name "Ernesto Paas"
```

We suggest saving a key pair, that can be generated with the command `ssh-keygen -t ed25519 -C "your_email@example.com".`
in the folder that contains the repo