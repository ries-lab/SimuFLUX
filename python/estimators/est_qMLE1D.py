import numpy as np

def est_qMLE1D(photonsi, L):
    if len(photonsi)>3:
        i1, i2, i0 = 2, 0, 4
    else:
        i1, i2, i0 = 0, 1, 2

    if isinstance(photonsi, list):
        photonsi = np.array(photonsi)

    maxrange=L
    # pi = photonsi/(photonsi[i1]+photonsi[i2])
    # eq 2.63
    xest1=-L/2+L/(1+np.sqrt(photonsi[i2]/photonsi[i1]))
    xest2=-L/2+L/(1-np.sqrt(photonsi[i2]/photonsi[i1]))
    
    p1 = pnorm(L, xest1)
    p2 = pnorm(L, xest2)
    LL1 = loglikelihood(p1/np.sum(p1),photonsi[[i1, i2, i0]])
    LL2 = loglikelihood(p2/np.sum(p2),photonsi[[i1, i2, i0]])

    if LL1 > LL2:
        xest=xest1
    else:
        xest=xest2
    xest[np.abs(xest)>maxrange] =np.nan  # deviation of Gaussian from Donut

    return xest

def pnorm(L, x0):
    return np.array([0.5+(2*L*x0)/(L**2+4*x0**2),
                     0.5-(2*L*x0)/(L**2+4*x0**2),
                     2*x0**2/(L**2+4*x0**2)])


def loglikelihood(pi, ni):
    return np.sum(ni*np.log(pi))