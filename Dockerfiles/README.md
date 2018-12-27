# GEM Maker Dockerfiles

Trying to implement best practices...

https://github.com/docker-library/official-images

## Build a Docker Image

```bash
docker build -t {TAG NAME} {PATH TO DOCKERFILE}

# example
cd Dockerfiles/sratoolkit/2.9.2
docker build -t systemsgenetics/sratoolkit:2.9.2 .

docker build -t systemsgenetics/base:1.0 .
docker build -t systemsgenetics/ascp:3.9.1 .

```

## Get Shell into Docker Image

```bash
docker run --rm -it systemsgenetics/sratoolkit:2.9.2 /bin/bash
docker run --rm -it systemsgenetics/ascp:3.9.1 /bin/bash

```
