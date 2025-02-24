# import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path

SCRIPT_DIR = Path(os.getcwd()).parent
sys.path.append(os.path.dirname(SCRIPT_DIR))

import numpy as np
from python.psfs import PsfVectorial

psf_vec = PsfVectorial() 

I_model = psf_vec.calculatePSFs('vortex', 0)