import numpy as np

from .Fluorophore import Fluorophore

class FlBleach(Fluorophore):
    def __init__(self, pos=..., brightness=1000):
        super().__init__(pos, brightness)
        self._photonbudget = np.inf

    def intensity(self, I0, dwelltime, time, phfac, props=None):
        remainingphotons = self.remainingphotons
        if remainingphotons == 0:
            return 0
        intensity = self.brightness*I0*dwelltime
        
        Io = np.minimum(intensity, remainingphotons)
        self.remainingphotons=remainingphotons-Io/phfac

        return Io
    
    @property
    def photonbudget(self):
        return self._photonbudget
    
    @photonbudget.setter
    def photonbudget(self, val):
        self._photonbudget = val
        self.reset(0)

    def reset(self, time):
        # resets the counter for bleaching, i.e. switches the fluorophore on again
        self.remainingphotons = np.random.exponential(self.photonbudget)
