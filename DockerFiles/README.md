# GEM Maker Dockerfiles

Trying to implemented best practices...

https://github.com/docker-library/official-images



**To build a DockerFile:**

```bash
docker build -t {TAG NAME} {PATH TO DOCKERFILE}
```

*Example. From the directory of the Dockerfile:*

```bash
docker build -t systemsgenetics/sratoolkit:2.8.2 .
```

**Enter Image Shell:**

```bash
docker run -it systemsgenetics/sratoolkit:2.9.2 /bin/bash
```
