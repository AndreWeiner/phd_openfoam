# OpenFOAM utilities and solvers

This repository contains [OpenFOAM](https://openfoam.com/) test cases, boundary conditions, utilities and solvers related to the PhD thesis entitled

> Modeling and simulation of convection-dominated species transfer at rising bubbles

to be completed by the end of 2019. Note that this repository is **work in progress**.

## Dependencies

All boundary conditions, utilities, and solvers are compiled using a special Docker image containing:

- OpenFOAM-v1906
- [PyTorch](https://pytorch.org/) 1.2.0

The Dockerfile and additional information on the build process can be found [here](https://github.com/AndreWeiner/of_pytorch_docker). The docker image is hosted on [Dockerhub](https://cloud.docker.com/u/andreweiner/repository/docker/andreweiner/of_pytorch). To pull the image containing OpenFOAM-v1906 and PyTorch 1.2, run

```
docker pull andreweiner/of_pytorch:of1906-py1.2-cpu
```

### Docker

Any installed version of [Docker](https://docs.docker.com/install/) larger than **1.10** will be able to pull and execute the Docker image hosted on [Dockerhub](https://hub.docker.com/r/andreweiner/jupyter-environment). There are convenience scripts to create and start a Docker container which require root privileges. To avoid running the scripts with *sudo*, follow the [post-installation steps](https://docs.docker.com/install/linux/linux-postinstall/).

### PyTorch models

### Geometry files

## Compiling solvers etc.

## Running test cases

```
source /opt/OpenFOAM/OpenFOAM-v1906/etc/bashrc
```

## How to reference
