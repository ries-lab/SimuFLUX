import numpy as np

def pol2cart(p, r):
    x = r*np.cos(p)
    y = r*np.sin(p)

    return x, y

def cart2pol(x, y):
    r = np.sqrt(x*x+y*y)
    p = np.arctan2(y,x)

    return p, r