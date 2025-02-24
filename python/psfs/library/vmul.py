#o=vmul(s,t)
#-----------
#
#Vector multiplication.
#
#  Leutenegger Marcel � 3.3.2005
#
#Input:
# s      vector
# t      scalar
#
#Output:
# o      vector: s*t = s.*repmat(t,3,1)
#

import numpy as np

def vmul(s=None, t=None):
    if s is None and t is None:
        print("\nVector multiplication.\n\n\tLeutenegger Marcel © 3.3.2005\n")
        return None
    elif s is not None and t is not None:
        if s.size == 0 or t.size == 0:
            return np.array([])
        elif s.shape[0] != 3 or t.shape[0] != 1:
            raise ValueError("Incompatible dimensions.")
        elif t.size == 1:
            return s * t
        elif s.size == 3:
            return np.tile(s, (t.size, 1)).T * np.tile(t, (3, 1))
        else:
            return s * np.tile(t, (3, 1))
    else:
        raise ValueError("Incorrect number of arguments.")
