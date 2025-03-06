import numpy as np

def est_donut2d(photonsi, patternpos, L, fwhm, background=0):
    # print(f"photonsi: {photonsi.shape}")
    # print(f"patternpos: {patternpos.shape}")
    # print(f"L: {L}")
    # print(f"fwhm: {fwhm}")
    # print(f"background: {background}")
    photonsi = photonsi - background / len(photonsi)
    
    pi = photonsi / np.sum(photonsi, axis=0)
    
    # Equation 2.63
    xest = -1 / (1 - (L**2 * np.log(2) / fwhm**2)) * np.sum(pi * patternpos,axis=0)
    
    return xest