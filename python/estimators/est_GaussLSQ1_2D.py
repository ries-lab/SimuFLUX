import numpy as np

def est_GaussLSQ1_2D(photonsi, patternpos, L, sigma, iscenter):
    pi = photonsi / np.sum(photonsi)
    Ls2 = L**2 / 8 / sigma**2
    if iscenter:
        numrim = patternpos.shape[0] - 1
        xest = (numrim + np.exp(Ls2)) / Ls2 / numrim * np.sum(pi * patternpos, axis=0)
    else:
        xest = np.sum(pi*patternpos, axis=0) / Ls2
    return xest
