# SimulFlux

A comprehensive simulator for MINFLUX experiments.

## MATLAB

Clone this repository. Open the repository in MATLAB. 

### Setup

Install MATLAB 2023b or newer with the curve fitting toolbox.

### Usage

Open MATLAB. In MATLAB, navigate to the `MATLAB/examples` folder.

## Python

### Setup

Clone this repository. Make sure `conda` is installed, ideally through 
[Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file).
Navigate to the folder containing this repository. Then run, in this folder,

```
cd python
conda create -n simulfux python=3.11
conda activate simulflux
pip install -r requirements.txt
```

### Usage

Launch a Jupyter Lab instance in VSCode, another IDE, or through the Miniforge prompt:

```
jupyter lab
```

In Jupyter Lab, open and run the notebooks in the `python/examples` folder.
