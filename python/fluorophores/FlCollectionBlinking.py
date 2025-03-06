import numpy as np

from ..tools import appendstruct
from ..tools import copyfields
from .FlCollection import FlCollection
from .FlBleach import FlBleach

class FlCollectionBlinking(FlCollection):
    def __init__(self):
        super().__init__()
        self.flprop = None # moving, isactive, nexton, nextoff, remaining_activations, startpos
        self.switchpar = {'starton': -1,
                          'tonsmlm': 1e3,
                          'toffsmlm': 1e5,
                          'photonbudget': 1000,
                          'activations': 1e6,
                          'brightness': 1000}  # starton, tonsmlm, toffsmlm, microseconds
        
    def position(self, time, props=None):
        flall = self.flall
        flprop = self.flprop
        pos = np.zeros((len(flall), 3))
        flmoving = np.where(flprop['moving'])[0]

        for fl in flmoving:
            pos[fl, :] = flall[fl].position(time)

        pos[~flprop['moving'],:] = flprop['startpos'][~flprop['moving'],:]

        return pos, flprop['isactive']
    
    def intensity(self, intin, dt, time, phfac, props=None):
        fact = np.where(self.flprop['isactive'])[0]
        ih = np.zeros((len(fact), 1))
        for k in range(len(fact)):
            ih[k] = self.flall[fact[k]].intensity(intin[k], dt, time, phfac[k])
        
        return ih
    
    def add(self, fllist):
        if isinstance(fllist, np.ndarray) and hasattr(fllist, '__add__') and hasattr(fllist, '__mul__'):  # call addstatic if poslist
            self.addstatic(fllist)
            return 
        if not isinstance(fllist, list):
            fllist = [fllist]

        try:
            self.flall.extend(fllist)
        except AttributeError:
            self.flall = fllist
        numfl = len(fllist)
        flproph = {}
        flproph['isactive'] = np.ones(numfl, dtype=bool)
        flproph['nexton'] = np.zeros(numfl)
        flproph['nextoff'] = np.zeros(numfl)
        flproph['moving'] = np.ones(numfl, dtype=bool)
        flproph['activations'] = np.zeros(numfl)
        flproph['remaining_activations'] = np.ones(numfl)
        self.flprop = appendstruct(self.flprop, flproph)
        self.numberOfFluorophores = len(self.flall)
        self.updatepos()

    def addstatic(self, poslist):
        if poslist.shape[1] == 2:  # 2D data
            poslistp = np.zeros((poslist.shape[0], poslist.shape[1]+1))
            poslistp[:,:-1] = poslist
            poslist = poslistp  # 3D data
        numfl = poslist.shape[0]
        flproph = {}
        flproph['moving'] = np.ones(numfl, dtype=bool)
        fllist = [None]*numfl
        for k in reversed(range(numfl)):
            flproph['moving'][k] = False

            addfl = FlBleach()
            addfl.pos = poslist[k,:]
            addfl.reset()
            fllist[k] = addfl
        self.flprop = appendstruct(self.flprop, flproph)
        try:
            self.flall.extend(fllist)
        except AttributeError:
            self.flall = fllist

        self.reset()  # reset fluorophores

    def updateonoff(self, time):
        flprop = self.flprop
        switchpar = self.switchpar
        isactive = flprop['isactive']

        # switch off
        switchoff = (time >= flprop['nextoff']) & flprop['isactive'] & (~flprop['moving'])
        # add remaining photons = 0?
        activef = np.where(isactive & (~flprop['moving']))[0]
        for k in activef:
            if self.flall[k].remainingphotons < 1:
                switchoff[k] = True
        ssoff = np.sum(switchoff)

        if ssoff > 0:
            isactive[switchoff] = False
            self.flprop['nexton'][switchoff] = np.random.exponential(switchpar['toffsmlm'], ssoff) + time

        # switch on
        switchon = np.where((time >= flprop['nexton']) \
                            & (~flprop['isactive'] \
                               & (~flprop['moving']) \
                                & (flprop['remaining_activations'] > 0)))[0]
        sson = len(switchon)
        if sson > 0:
            isactive[switchon] = True
            self.flprop['activations'][switchon] += 1
            self.flprop['nextoff'][switchon] = np.random.exponential(switchpar['tonsmlm'], sson) + time
            self.flprop['remaining_activations'][switchon] -= 1
            for k in range(sson):
                self.flall[switchon[k]].reset()

        if (sson > 0) or (ssoff > 0):
            self.flprop['isactive'] = isactive
    
    def setpar(self, **kwargs):
        self.switchpar, _ = copyfields(self.switchpar, kwargs)

    def reset(self, time=None):
        if time is None:
            time = 0
        flproph = self.flprop
        flall = self.flall
        numfl = len(flall)
        switchpar = self.switchpar

        # check for each element
        try:
            if numfl != flproph['moving'].shape[0]:
                moving = np.ones(numfl, dtype=bool)
                moving[:flproph['moving'].shape[0]] = flproph['moving']
                flproph['moving'] = moving
        except KeyError:
            flproph['moving'] = np.ones(numfl, dtype=bool)
        try: 
            if numfl != flproph['isactive'].shape[0]:
                moving = np.ones(numfl, dtype=bool)
                moving[:flproph['isactive'].shape[0]] = flproph['isactive']
                flproph['isactive'] = moving
        except KeyError:
            flproph['isactive'] = np.ones(numfl, dtype=bool)
        try:
            if numfl != flproph['nexton'].shape[0]:
                moving = np.zeros(numfl)
                moving[:flproph['nexton'].shape[0]] = flproph['nexton']
                flproph['nexton'] = moving
        except KeyError:
            flproph['nexton'] = np.zeros(numfl)
        try:
            if numfl != flproph['nextoff'].shape[0]:
                moving = np.zeros(numfl)
                moving[:flproph['nextoff'].shape[0]] = flproph['nextoff']
                flproph['nextoff'] = moving
        except KeyError:
            flproph['nextoff'] = np.zeros(numfl)
        

        for k in reversed(range(numfl)):
            flall[k].reset()
            if np.any(flproph['moving'][k]):
                continue
            if switchpar['starton'] > -1:
                flproph['isactive'][k] = switchpar['starton'] == 1
            else:
                flproph['isactive'][k] = np.random.rand() < (switchpar['tonsmlm']/(switchpar['toffsmlm'] + switchpar['tonsmlm']))
    
            flall[k].brightness = switchpar['brightness']
            flall[k].photonbudget = switchpar['photonbudget']

            flproph['nexton'][k] = np.random.exponential(switchpar['toffsmlm']) + time
            flproph['nextoff'][k] = np.random.exponential(switchpar['tonsmlm']) + time
        
        flproph['activations'] = flproph['isactive'].astype(float)
        flproph['remaining_activations'] = np.maximum(1,np.random.poisson(switchpar['activations'], numfl))
        self.flprop = flproph
        self.numberOfFluorophores = len(self.flall)
        self.updatepos()

    @property
    def allbleached(self):
        return np.sum(self.flprop['remaining_activations']) == 0

    def updatepos(self):
        fl = self.flall
        self.flprop['startpos'] = np.zeros((len(fl), 3))
        for k in reversed(range(len(fl))):
            pos, _ = fl[k].position(0)
            self.flprop['startpos'][k,:] = pos