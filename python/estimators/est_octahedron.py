import numpy as np

from .est_donutLSQ1_2D import est_donutLSQ1_2D
from .est_qLSQiter1D import est_qLSQiter1D

def est_octahedron(photonsi, patternposi, L, fwhm, background=0):
    photonsxy = photonsi[0:4]
    photonsz = photonsi[4:6]
    patternposxy = patternposi[0:4,:]

    if len(photonsi) == 7:
        # center probed
        photonsxy = np.vstack([photonsxy, photonsi[6]])
        photonsz = np.vstack([photonsz, photonsi[6]])
        patternposxy = np.vstack([patternposxy,
                                  np.zeros(patternposxy.shape[1], 
                                           dtype=patternposxy.dtype)])

    xest = est_donutLSQ1_2D(photonsxy, patternposxy, L, fwhm, background)
    xest[2] = est_qLSQiter1D(photonsz, L)

    return xest
