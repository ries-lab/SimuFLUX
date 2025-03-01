import numpy as np

class FlCollection:
    def __init__(self):
        self.flall = None
        self.numberOfFluorophores = 0
    
    def set(self, fllist):
        if not isinstance(fllist, list):
            fllist = [fllist]
        self.flall = fllist
        self.numberOfFluorophores = len(fllist)
    
    def add(self, fllist):
        if not isinstance(fllist, list):
            fllist = [fllist]

        if self.flall is None:
            self.flall = fllist
        else:
            self.flall.extend(fllist)
        self.numberOfFluorophores = len(self.flall)

    def position(self, time, props=None):
        pos = np.zeros((self.numberOfFluorophores, 3))

        for k in range(self.numberOfFluorophores):
            pos[k,:], _ = self.flall[k].position(time)
        
        isactive = np.ones(self.numberOfFluorophores).astype(bool)

        return pos, isactive
    
    def intensity(self, intin, dt, time, phfac, props=None):
        ih = np.zeros((self.numberOfFluorophores,1))
        for k in range(self.numberOfFluorophores):
            ih[k] = self.flall[k].intensity(intin[k], dt, time, phfac[k])
        
        return ih

    def reset(self, time):
        pass

    def updateonoff(self, time):
        pass

    @property
    def allbleached(self):
        return False
    
    def getproperties(self):
        return []
    
    @property
    def remainingphotons(self):
        out = 0
        for k in range(self.numberOfFluorophores):
            out += self.flall[k].remainingphotons
        return out
    
    @property
    def brightness(self):
        out = np.zeros((self.numberOfFluorophores, 1))
        for k in range(self.numberOfFluorophores,-1,-1):
            out[k] = self.flall[k].brightness
        return out    
