# Building Physics - Applications in Python

This repository contains the code examples quoted in the book. A pdf copy of this book can be downloaded [here](link).

## Environment setup

The Python environment of this repository can be set up via Pipenv. 
To install Pipenv on your machine, we refer you to the [Pipenv documentation](https://pipenv-fork.readthedocs.io/en/latest/install.html).

```shell script
pipenv install
pipenv shell
``` 
The first command installs the project dependencies.
The second command starts a pipenv shell with the virtual env activated.

If the dependencies locking takes a long time, use :

```
PIPENV_SKIP_LOCK=true pipenv install 
```

## Starting a notebook

Examples are provided as python scripts and interactive Notebooks.
To start a notebook server :

``` shell script
pipenv shell
jupyter notebook --ip=0.0.0.0
```

A browser should start with a Jupyter server showing the root folder of the repository.

## Structure of the repository

The files are organised as in the book :

### Chapter 1 : basics
- Explicit Euler Vector
- Explicit Eurler Matrix
- Crank Nicolson example for air interface

### Chapter 2 : PCM, HVAC
- PCM ()
- HVAC 

### Chapter 3 : coupled problems and minimization
