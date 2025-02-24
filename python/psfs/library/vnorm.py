#o=vnorm(s)
#----------
#
#Vector norm.
#
#  Leutenegger Marcel � 3.3.2005
#
#Input:
# s      vector
#
#Output:
# o      scalar: |s|^2 = dot(s,s)
#

import numpy as np

def vnorm(s=None):
    if s is None:
        print("\nVector norm.\n\tLeutenegger Marcel © 3.3.2005\n")
        return None
    elif s is not None:
        if s.size == 0:
            return np.array([])
        elif s.shape[0] != 3:
            raise ValueError("Incompatible dimensions.")
        else:
            return np.sum(np.real(np.conj(s) * s), axis=0)
    else:
        raise ValueError("Incorrect number of arguments.")
