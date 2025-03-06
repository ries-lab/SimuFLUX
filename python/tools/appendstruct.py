from types import SimpleNamespace

import numpy as np

def appendstruct(sin, sadd):
    if not sin:
        return sadd
    
    is_dict = True
    try:
        sout = sin.copy()
    except AttributeError:
        # SimpleNamespace
        is_dict = False  # keep types consistent
        sout = sin.__dict__

    fn = sout.keys()
    try:
        sadd.keys()
    except AttributeError:
        sadd = sadd.__dict__

    for key in fn:
        if key in sadd.keys():
            sout[key] = np.concatenate((sout[key], sadd[key]))                

    if is_dict:
        return sout
    else:
        return SimpleNamespace(**sout)