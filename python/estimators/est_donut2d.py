import numpy as np

def est_donut2d(photonsi, patternpos, L, fwhm, background=0):
    photonsi = photonsi - background / len(photonsi)
    
    pi = photonsi / np.sum(photonsi, axis=0)
    
    # Equation 2.63
    xest = -1 / (1 - (L**2 * np.log(2) / fwhm**2)) * np.sum(pi * patternpos,axis=0)
    
    return xest