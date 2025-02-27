## This function introduces a phase shift at the center of the input image.
# Code written by Takahiro DEGUCHI, European Molecular Biology Laboratory, October 2023.
# This is meant to be used with MATLAB library, "Electromagnetic field in the focus region of a microscope objective", 
# based on the original publication, Leutenegger et al., "Fast focus field calculations", published in Optics Express 14 (2006) 

import numpy as np

from .math import pol2cart, cart2pol

RADSHIFT = 1/np.sqrt(2) + 0.03

def pishift(sys, E, r, t, p):
    xx, yy = pol2cart(p, r)
    xx = xx + sys['maskshift'][0]
    yy = yy + sys['maskshift'][1]
    p, r = cart2pol(xx, yy)
    
    idx = r <= RADSHIFT # The area for pi shift.
    idx = 2*idx-1
    phase = idx*np.exp((np.pi)*1j) 
    if "delshift" in sys.keys():
        phase[r<=RADSHIFT] *= np.exp(sys['delshift']*1j)

    E *= phase

    return E
