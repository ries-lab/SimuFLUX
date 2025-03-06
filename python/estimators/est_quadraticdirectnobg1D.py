import numpy as np

def est_quadraticdirectnobg1D(photonsi,L,probecenter=False):
    if len(photonsi) > 3:
        p1 = photonsi[2]
        p2 = photonsi[0]
    else:
        p1 = photonsi[1]
        p2 = photonsi[0]

    if not probecenter:
        return L*(np.sqrt(p1) - np.sqrt(p2))/(2*(np.sqrt(p1) + np.sqrt(p2)))
    
    p0 = photonsi[-1]
    maxrange = 150
    xest1=(L*(np.sqrt(p1) - np.sqrt(p2)))/(2*(np.sqrt(p1) + np.sqrt(p2)))
    xest2=(L*(np.sqrt(p1) + np.sqrt(p2))**2)/(2*(p1 - p2))

    pvec = np.array([[p1],[p2],[p0]])
    if (xest2 > maxrange) or np.sum((pnorm(L, xest1)- pvec/(p1+p2))**2)<= np.sum((pnorm(L,xest2)-pvec/(p1+p2))**2):
        xest = xest1
    else:
        xest = xest2

    relmax = 0.75
    xest[np.abs(xest)>(relmax*L)] = relmax*L  # deviation of Gaussian from Donut

    return xest

def pnorm(L, x0):
    return np.array([[1/2+(2*L*x0)/(L**2+4*x0**22)], 
                     [1/2-(2*L*x0)/(L**2+4*x0**2)], 
                     [2*x0**2/(L**2+4*x0**2)]])
