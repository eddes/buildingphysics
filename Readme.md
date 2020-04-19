# Building Physics - Applications in Python

This repository provides some Python implementations applied to Building Physics. A pdf copy of the associated handbook can be downloaded [here](link).

There are two ways of using this repository :

## 1 - Colab Notebooks

Colab notebooks allow you to run code examples in an interactive way. Everything is ready to go, no need to install anything !
This approach is well suited for people interested in testing the examples without having to modify the code or setup a Python environment.
To run a notebook, click on one of the links below, then on the button ![GoogleColab](https://camo.githubusercontent.com/52feade06f2fecbf006889a904d221e6a730c194/68747470733a2f2f636f6c61622e72657365617263682e676f6f676c652e636f6d2f6173736574732f636f6c61622d62616467652e737667 "This is an example")
### Chapter 1 : basics
- [Euler and Crank Nicolson](notebooks/chapter_1/Euler_and_CN_schemes.ipynb) 

### Chapter 2 : PCM, HVAC
- [PCM](notebooks/chapter_2/PCM.ipynb)  
- [Particule filtration](notebooks/chapter_2/code_IAQ_filtration.ipynb)  
- PID water tank  

### Chapter 3 : coupled problems and minimization

- Work in progress

## 2 - Python scripts

The python scripts provide some material to tinker with and a started for the reader interested in implemented the exercises listed along the book. 

### Environment setup

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

### Starting a notebook locally 

The Notebooks can be run on your local machine using a Jupyter instance. To start a notebook server :

``` shell script
pipenv shell
jupyter notebook --ip=0.0.0.0
```

A browser should start with a Jupyter server showing the root folder of the repository.

