# GEM Maker Dockerfiles

Trying to implement best practices...

https://github.com/docker-library/official-images

## Build a Docker Image

```bash
docker build -t {TAG NAME} {PATH TO DOCKERFILE}

# example
docker build -t systemsgenetics/gemmaker     Dockerfiles/base/1.0/
docker build -t systemsgenetics/aspera:3.8.1 Dockerfiles/aspera/3.8.1/
```

## Get Shell into Docker Image

```bash
docker run --rm -it systemsgenetics/aspera:3.8.1 /bin/bash
```
