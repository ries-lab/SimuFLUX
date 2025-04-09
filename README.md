# SimuFlux

A comprehensive simulator for MINFLUX experiments.

**Reference:**
[Marin, Z, and J Ries. Evaluating MINFLUX Experimental Performance *in silico*. bioRxiv, 2025.04.08.647786 (2025).](https://doi.org/10.1101/2025.04.08.647786)

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
conda create -n simuflux python=3.11
conda activate simuflux
pip install -r requirements.txt
```

### Usage

Launch a Jupyter Lab instance in VSCode, another IDE, or through the Miniforge prompt:

```
jupyter lab
```

In Jupyter Lab, open and run the notebooks in the `python/examples` folder.

## Google Colab

The Python notebooks can be run from Google Colab without installing any software.

Please note that you may occasionally experience an error in one of the notebook cells. 
If this happens, the notebook will stop running. In the event of an error, please navigate 
to Runtime > Restart session and start the notebook run from the beginning, following 
instructions. If the same error occurs twice, please contact the authors via the "Issues" 
tab at https://github.com/ries-lab/SimuFLUX/.

Please note that if you see a "Runtime Warning" message, this is not an error. The notebook 
will not stop running and this is not a problem.


| Notebook Name | Link | Description |
|---------------|------|-------------|
| example1_simple_MINFLUX.ipynb | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1lAS3fMO2pxzLy0JPsQ5AC4pMNMK8Ar95) | Scan a static fluorophore with a donut PSF. |
| example2_vectorial_PSF.ipynb    | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1lgm1Kpv_XJT8oJnX_mBtB1QmPyjOgoew) | Use a vectorial PSF to examine the influence of misalignment, background, and fluorescent beads on measurements. Try PhaseFLUX. |
| example2b_compare_PSFs.ipynb | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1IeVmwU5EaQTCq2FOnU_bcnPgwNQDf8bL) | Compare different excitation PSFs. Examine the effects of multiple fluorophores. |
| example2c_phaseplate_misalignment.ipynb  | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1u_FvSvzCX96bSlv1O7dAh2V8rlAZbo8b) | Look at phase plate misalignment. |
| example2d_pinhole_misalignment.ipynb  | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1CqLoiKFxxCd5wf4iGHN29yTSV5xtj7-e) | Look at pinhole misalignment. |
| example3_Abberior_sequence.ipynb  | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/19Qcq7qlGfqDA2Ti0_7Q8SxPwoc2XLDtk) | Use an Abberior sequence file to run an experiment. |
| example4_blinking_fluorophore.ipynb  | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1OQ0fQ_enYYOtbP44-mrVKRQ0CTl7UKwo) | Simulate measurement with a blinking fluorophore. Investigate averaging of flickering signal. |
| example5_moving_fluorophore.ipynb | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1pqMiLu17PSFbDDta0xxJ2tV6g0RHnFu4) | Simulate measurement with a moving fluorophore. Investigate diffusion and system vibrations. |
| example5b_max_diffusion.ipynb | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1-5EYJim_54Ebv2-uYzoMQIew9J5-rPJT) | Simulate diffusion and investigate system properties and root mean square error as a function of diffusion coefficient. |
| example6_Fluorophore_Collections.ipynb | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1AJv-4grNjnR9jf6jzghaCDeOe19fkL4d) | Simulate MINFLUX with multiple fluorophores. Image a simulated nuclear pore complex. |
| example7_estimators.ipynb | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1mLMCshwX_8bQWWLvQ7DytabVPJmn90dc) | Investigate the performance of different MINFLUX estimators. |
| example8_background.ipynb | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1xMTNvXxg_iqg_K3vMmkSQCgATMIRwgXj) | Investigate the influence of background from constant offsets, autofluorescence, and nearby fluorophores. |
| example10_imaging.ipynb | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1lbB7ysgNL9Zuzqe_6V4CssM6SAFrZg2W) | Simulate MINFLUX imaging with DNA-PAINT and dSTORM. Tune fluorophore densities to optimal levels. |
