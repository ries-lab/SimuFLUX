import numpy as np

from .est_qLSQiter1D import est_qLSQiter1D

def est_zline(photonsi, L, iter=15, eps=0.1):
    print("zline", L)
    phot1 = photonsi[0:2]
    phot2 = photonsi[2:4]
    if len(photonsi) == 5:
        phot1 = np.vstack([phot1, photonsi[4]])
        phot2 = np.vstack([phot2, photonsi[4]])

    c = np.array([2.1415, -3.2122, 2.0062, -0.4677])*1e3
    fr = phot2[1]/(phot2[0]+phot2[1])
    xest = c[0]*fr**3+c[1]*fr**2+c[2]*fr+c[3]

    if np.abs(xest) < (L[0]/2*0.75):
        xest = est_qLSQiter1D(phot1, L[0], iter, eps, xest)

    return xest