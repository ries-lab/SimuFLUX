from types import SimpleNamespace

import numpy as np

def appendstruct(sin, sadd):
    if not sin:
        return sadd
    
    sout = sin.copy()
    try:
        fn = list(sin.keys())  # Get field names (keys) of the dictionary
    except AttributeError:
        fn = list(sin.__dict__.keys())
    # ln = len(sadd[fn[0]])  # Get the length of the first field in sadd

    for key in fn:
        if key in sadd.keys():
            sout[key] = np.concatenate((sout[key], sadd[key]))

    return SimpleNamespace(**sout)