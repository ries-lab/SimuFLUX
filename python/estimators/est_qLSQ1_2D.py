import numpy as np

def est_qLSQ1_2D(photonsi, patternpos):
    pi=photonsi/np.sum(photonsi)
    # eq 2.63
    xest=-np.sum(pi*patternpos)
    return xest