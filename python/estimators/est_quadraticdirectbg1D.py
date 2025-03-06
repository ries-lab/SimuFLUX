import numpy as np

def est_quadraticdirectbg1D(photonsi, L):
    if len(photonsi) > 4:
        i1, i2, i0 = 2, 0, 4
    elif len(photonsi) == 3:
        i1, i2, i0 = 0, 1, 2
    elif len(photonsi) == 2:
        return q2point(photonsi, L)
    
    p = photonsi[np.array([i1, i2, i0])]  # for testing
    xest = (-(L*p[0]) + L*p[1])/(4*(2*p[2] - p[0] - p[1]))

    relmax=0.75
    xest[np.abs(xest)>(relmax*L)] = relmax*L  # extrapolation does not work well

    return xest

def q2point(photons, L):
    raise NotImplementedError("I don't know how to handle two points.")