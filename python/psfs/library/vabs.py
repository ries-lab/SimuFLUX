#o=vabs(s)
#---------
#
#Vector length.
#
#  Leutenegger Marcel � 4.12.2005
#
#Input:
# s      vector
#
#Output:
# o      scalar: |s| = sqrt(dot(s,s))
#

import numpy as np

def vabs(s=None):
    if s is None:
        print("\nVector length.\n\n\tLeutenegger Marcel © 4.12.2005\n")
        return None
    elif s is not None:
        if s.size == 0:
            return np.array([])
        elif s.shape[0] != 3:
            raise ValueError("Incompatible dimensions.")
        else:
            return np.sqrt(np.sum(np.real(np.conj(s) * s)))
    else:
        raise ValueError("Incorrect number of arguments.")
