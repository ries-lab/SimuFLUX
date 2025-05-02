# SimuFlux

A comprehensive simulator for MINFLUX experiments.

**Reference:**
[Marin, Z, and J Ries. Evaluating MINFLUX Experimental Performance *in silico*. bioRxiv, 2025.04.08.647786 (2025).](https://doi.org/10.1101/2025.04.08.647786)

## MATLAB

### Requirements

* Windows 10+, Mac 13.7+ or Linux (Ubuntu 20.04+). 8 GB RAM.
* [MATLAB](https://www.mathworks.com/products/matlab.html) 2023b+ with the curve fitting toolbox.

### Setup

Clone this repository. Open the repository in MATLAB. 

### Usage

Open MATLAB. In MATLAB, navigate to the `MATLAB/examples` folder.

### Technical Details

Tested on Windows 11 and Mac OS 14.5. Installation time for MATLAB is ~1 hour.
Installation time for SimuFLUX is < 5 min.

## Python


### Requirements

* Windows 10+, Mac 13.7+ or Linux (Ubuntu 20.04+). 8 GB RAM.
* [Miniforge](https://github.com/conda-forge/miniforge?tab=readme-ov-file) or an equivalent `conda` environment manager.

### Setup

Clone this repository. Navigate to the folder containing this repository. Then run, in this folder,

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

### Technical Details

Tested on Windows 11 and Mac OS 14.5. Installation time for Miniforge is ~20 minutes.
Installation time for SimuFLUX is < 5 min.

## Google Colab

The Python notebooks can be run from Google Colab without installing any software.

Please note that you may occasionally experience an error in one of the notebook cells. 
If this happens, the notebook will stop running. In the event of an error, please navigate 
to Runtime > Restart session and start the notebook run from the beginning, following 
instructions. If the same error occurs twice, please contact the authors via the ["Issues" 
tab](https://github.com/ries-lab/SimuFLUX/issues).

Please note that if you see a "Runtime Warning" message, this is not an error. The notebook 
will not stop running and this is not a problem.


| Notebook Name | Description | Link |
|---------------|-------------|------|
| example1_simple_MINFLUX.ipynb | Scan a static fluorophore with a donut PSF. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1lAS3fMO2pxzLy0JPsQ5AC4pMNMK8Ar95) |
| example2_vectorial_PSF.ipynb | Use a vectorial PSF to examine the influence of misalignment, background, and fluorescent beads on measurements. Try PhaseFLUX. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1lgm1Kpv_XJT8oJnX_mBtB1QmPyjOgoew) |
| example2b_compare_PSFs.ipynb | Compare different excitation PSFs. Examine the effects of multiple fluorophores. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1IeVmwU5EaQTCq2FOnU_bcnPgwNQDf8bL) |
| example2c_phaseplate_misalignment.ipynb  | Look at phase plate misalignment. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1u_FvSvzCX96bSlv1O7dAh2V8rlAZbo8b) |
| example2d_pinhole_misalignment.ipynb  | Look at pinhole misalignment. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1CqLoiKFxxCd5wf4iGHN29yTSV5xtj7-e) |
| example3_Abberior_sequence.ipynb  | Use an Abberior sequence file to run an experiment. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/19Qcq7qlGfqDA2Ti0_7Q8SxPwoc2XLDtk) |
| example4_blinking_fluorophore.ipynb  | Simulate measurement with a blinking fluorophore. Investigate averaging of flickering signal. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1OQ0fQ_enYYOtbP44-mrVKRQ0CTl7UKwo) |
| example5_moving_fluorophore.ipynb | Simulate measurement with a moving fluorophore. Investigate diffusion and system vibrations. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1pqMiLu17PSFbDDta0xxJ2tV6g0RHnFu4) |
| example5b_max_diffusion.ipynb | Simulate diffusion and investigate system properties and root mean square error as a function of diffusion coefficient. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1-5EYJim_54Ebv2-uYzoMQIew9J5-rPJT) |
| example6_Fluorophore_Collections.ipynb | Simulate MINFLUX with multiple fluorophores. Image a simulated nuclear pore complex. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1AJv-4grNjnR9jf6jzghaCDeOe19fkL4d) |
| example7_estimators.ipynb | Investigate the performance of different MINFLUX estimators. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1mLMCshwX_8bQWWLvQ7DytabVPJmn90dc) |
| example8_background.ipynb | Investigate the influence of background from constant offsets, autofluorescence, and nearby fluorophores. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1xMTNvXxg_iqg_K3vMmkSQCgATMIRwgXj) |
| example10_imaging.ipynb | Simulate MINFLUX imaging with DNA-PAINT and dSTORM. Tune fluorophore densities to optimal levels. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1lbB7ysgNL9Zuzqe_6V4CssM6SAFrZg2W) |
| example11_max_diffusion.ipynb | Simulate tracking of diffusion fluorophore under different conditions. Optimize MINFLUX parameters for tracking. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1_awyZfaorz7JCf5moyi-jgTdIWwZjXRs) |
