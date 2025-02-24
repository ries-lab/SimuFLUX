#Phase ramp.
#
# sys    System data
#
#
# E      Electric field [V/m]
# r      Radial position [Da/2]
# t      Incidence angle [rad]
# p      Polar angle [rad]
#

import numpy as np

from .math import pol2cart, cart2pol

def phaseramp(sys, E, r, t, p):
    xx, yy = pol2cart(p, r)
    xx = xx + sys['maskshift'][1]
    yy = yy + sys['maskshift'][0]
    p, r = cart2pol(xx, yy)
    
    e = np.exp(1j * p)
    try:
        E = E * np.tile(e, (E.shape[0], 1))
    except AttributeError:
        E *= e

    return E
