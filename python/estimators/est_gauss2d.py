import numpy as np

def est_gauss2d(photonsi, patternpos, L, sigma, iscenter):
    numrim = patternpos.shape[0] - iscenter
    pi = photonsi / np.sum(photonsi)
    Ls2 = L**2 / 8 / sigma**2
    xest = (numrim + np.exp(Ls2)) / Ls2 / numrim * np.sum(pi * patternpos, axis=0)
    return xest
