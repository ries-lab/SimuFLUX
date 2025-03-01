import numpy as np

def est_quadraticdirect1D(photonsi, L):
    if len(photonsi) > 4:
        i1, i2, i0 = 2, 0, 4
    else:
        i1, i2, i0 = 0, 1, 2
    
    maxrange = 150

    pi = photonsi / (photonsi[i1]+photonsi[i2])

    # eq 2.63
    xest1 = (L-np.sqrt(L**2*(1-(pi[i1]-pi[i2])**2)))/(2*(pi[i1]-pi[i2]))
    xest2 = (L+np.sqrt(L**2*(1-(pi[i1]-pi[i2])**2)))/(2*(pi[i1]-pi[i2]))

    if (xest2 > maxrange) or np.sum((pnorm(L,xest1)-pi[[i1,i2,i0]])**2) < np.sum((pnorm(L,xest2)-pi[[i1,i2,i0]])**2):
        xest = xest1
    else:
        xest = xest2

    xest[np.abs(xest)>maxrange] = np.nan  # deviation of Gaussian from Donut

    return xest    

def pnorm(L, x0):
    po = np.array([
        1/2 + (2 * L * x0) / (L**2 + 4 * x0**2),
        1/2 - (2 * L * x0) / (L**2 + 4 * x0**2),
        2 * x0**2 / (L**2 + 4 * x0**2)
    ])
    return po
