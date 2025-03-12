import numpy as np

def est_quad1Diter(photonsi, L, iter=15, eps=0.1):
    """
    Eilers 2.64, k=1
    """

    if len(photonsi) == 2:
        itfun = iteration
    else:
        itfun = iterationcenter

    pi = photonsi / np.sum(photonsi)
    xest = np.array([0.0])
    for _ in range(iter):
        xo = xest.copy()
        xest = itfun(pi, xest, L)
        if np.sum((xest - xo)**2) < eps**2:
            break
    # xest = np.clip(xest, -L, L) # avoid crazy high numbers
    xest[np.abs(xest)>1.3*L] = np.nan
    return xest

def iteration(p, x0, L):
    d12 = p[0] - p[1]
    dr = ((L**2 + 4 * x0**2) * (-4 * L * x0 + d12 * (L**2 + 4 * x0**2))) / (4 * (L**3 - 4 * L * x0**2))
    dro = dr + x0
    return dro

def iterationcenter(p, x0, L):
    # p[0] at -L/2, p[1] at L/2, p[2] at 0
    d12 = p[0] - p[1]
    pt = p[0] + p[1] - 2 * p[2]
    dr = ((L**2 + 6 * x0**2) * (-L**3 * (3 + pt) * x0 - 6 * L * (-4 + pt) * x0**3 + 
           d12 * (L**4 - 36 * x0**4))) / (4 * (L**5 - 9 * L**3 * x0**2 + 36 * L * x0**4))
    dro = dr + x0
    return dro
