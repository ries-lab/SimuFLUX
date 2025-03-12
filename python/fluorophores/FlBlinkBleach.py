import numpy as np

from .Fluorophore import Fluorophore

class FlBlinkBleach(Fluorophore):
    def __init__(self, pos=[0,0,0], brightness=1000):
        super().__init__(pos, brightness)

        self.fast_ton = 0.1 # ms
        self.fast_toff = 0
        self.starton = False
        self.photonbudget = np.inf  # 1000 photons
        self.time = 0
        self.blinkingtrace = np.array([0, 0, 0])
        self.tind = 0

    def intensity(self, I0, dwelltime, time, phfac, props=None):
        fraction = self.measure(dwelltime, time)
        intensity = (self.brightness)*I0*dwelltime*fraction
        remainingphotons = self.remainingphotons
        Io = np.minimum(intensity, remainingphotons)
        self.remainingphotons = remainingphotons-Io/phfac

        return Io

    def makeblinkingtrace(self):
        self.blinkingtrace = calculatetrace(self.fast_ton, self.fast_toff)
        self.tind = 0
        self.time = 0

    def reset(self, time=None):
        self.makeblinkingtrace()
        if not self.starton:
            self.blinkingtrace[:,2] -= np.random.rand()*self.blinkingtrace[-1,2]/2
        
        self.blinkingtrace[:,2] += time
        # recalculate photon budget
        self.remainingphotons = np.random.exponential(self.photonbudget)

    def measure(self, dwelltime, t):
        trace = self.blinkingtrace
        
        if t>(trace[-1,2]-dwelltime*2):  #outside range
            self.extendblinkingtrace()  # extend blinkingtrace
    
        ind1 = self.tind
        if trace[ind1,2] > t:  #time asked that was before
            result = np.where(trace[ind1,2] <= t)[0]
            if len(result) > 0:
                ind1 = result[-1]
            else:
                ind1 = 0

        while trace[ind1,2] <= t:  # find last index below
            ind1 += 1
        ind1 = np.maximum(0, ind1-1)
        ind2 = ind1
        while trace[ind2,2] < t+dwelltime:
            ind2 += 1
        ind2 = np.maximum(0, ind2-1)
        ontime = np.sum(trace[ind1:ind2,0])
        ontime -= np.minimum(t-trace[ind1,2], trace[ind1,0])  # first block
        ontime += np.minimum(t+dwelltime-trace[ind2,2], trace[ind2,0])  # last block
        fraction = np.maximum(ontime/dwelltime, 0)

        self.tind = ind2

        return fraction
    
    def extendblinkingtrace(self):
        blt2 = calculatetrace(self.fast_ton, self.fast_toff)
        blt2[:,2] += np.sum(self.blinkingtrace[-1,:])
        self.blinkingtrace = np.vstack([self.blinkingtrace, blt2])

def calculatetrace(ton,toff):
    maxsimultime = 100*(ton+toff)
    numpoints = np.maximum(10, np.ceil(maxsimultime/(ton+toff)*2)).astype(int)
    tont = np.random.exponential(ton, (numpoints,1))
    tofft = np.random.exponential(toff, (numpoints,1))
    tsum = tont+tofft
    timestart = np.cumsum(np.vstack([[0], tsum[:-1]]), axis=0)
    # print(f"tont: {tont.shape} tofft: {tofft.shape} timestart: {timestart.shape} tsum: {tsum.shape}")
    blinkingtrace = np.hstack([tont,tofft,timestart])
    # print("blinkingtrace ", blinkingtrace.shape)

    return blinkingtrace
