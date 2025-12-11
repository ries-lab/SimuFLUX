from .est_GaussLSQ1_2D import est_GaussLSQ1_2D

def est_pinholeorbit(photonsi, patternpos, L, sigma, iscenter):
    xesth = est_GaussLSQ1_2D(photonsi, patternpos, L, sigma, iscenter)
    xest = -xesth

    return xest