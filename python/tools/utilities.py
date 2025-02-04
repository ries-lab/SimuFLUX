
"""
Copyright (c) 2022      Ries Lab, EMBL, Heidelberg, Germany
All rights reserved     Heintzmann Lab, Friedrich-Schiller-University Jena, Germany

@author: Rainer Heintzmann, Sheng Liu, Jonas Hellgoth
"""

import numpy as np
import scipy.signal as signal

def fft3d(tfin):
    return np.fft.fftshift(np.fft.fftn(np.fft.fftshift(tfin, axes=(-1, -2, -3))), axes=(-1, -2, -3))

def ifft3d(tfin):
    return np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(tfin, axes=(-1, -2, -3))), axes=(-1, -2, -3))

def cztfunc1(datain, param):
    A = param[0]
    Bh = param[1]
    C = param[2]
    N = A.shape[0]
    L = Bh.shape[0]
    M = C.shape[0]

    # Pad the input arrays
    Apad = np.pad(A * datain / N, ((0, 0), (0, L - N)), mode='constant')
    Apad = np.pad(Apad, ((0, L - N), (0, 0)), mode='constant')
    
    # Perform 2D FFT
    Ah = np.fft.fft2(Apad)
    
    # Perform 2D IFFT
    cztout = np.fft.ifft2(Ah * Bh / L)
    
    # Extract the relevant portion of the output
    dataout = C * cztout[..., -M:, -M:]

    return dataout
