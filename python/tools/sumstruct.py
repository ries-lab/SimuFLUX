import copy

import numpy as np

def sumstruct(sin=None, sadd=None):
    if sin is None or getattr(sin, 'phot') is None:
        return sadd
    if sadd is None or getattr(sadd, 'phot') is None:
        return sin
    
    sout = sin
    fn = sin.__dict__.keys()
    for k in fn:
        saddk = getattr(sadd, k)
        if isinstance(saddk, (int, float)) or (isinstance(saddk, np.ndarray) and hasattr(saddk, '__add__') and hasattr(saddk, '__mul__')):
            # print(f"key: {k} sink: {getattr(sin, k)} saddk: {getattr(sadd, k)}")
            setattr(sout, k, getattr(sin, k) + getattr(sadd, k))
        # else:
        #     setattr(sout, k, getattr(sin, k))

    return sout