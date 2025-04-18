from collections import namedtuple
import numpy as np

Properties = namedtuple("Properties", ["brightness", "pos"])

class Fluorophore:
    def __init__(self, pos=[0,0,0], brightness=1000):
        self.pos = pos  # nm
        self.brightness = brightness  # kHz
        self.remainingphotons = np.inf
        # Dummy properties for compatibility with collection
        self.numberOfFluorophores = 1
        self.allbleached = False

    @property
    def pos(self):
        return self._pos
    
    @pos.setter
    def pos(self, val):
        if not isinstance(val, np.ndarray):
            self._pos = np.array(val, dtype=float)
        else:
            self._pos = val
    
    def intensity(self, I0, dwelltime, time, phfac):
        # if props is not None:  # For performance: avoid property access
        #     brightness = props.brightness
        # else:
        brightness = self.brightness
        
        return brightness * I0 * dwelltime #* phfac
    
    def photons(self, I0, *args):
        return np.random.poisson(self.intensity(I0, *args))
    
    def reset(self, time=None):
        pass
    
    def position(self, time):
        pos = self.pos
        
        return pos, True
    
    # Dummy functions for compatibility with collection
    def updateonoff(self, time):
        pass
    
    def getproperties(self):
        return Properties(brightness=self.brightness, pos=self.pos)
