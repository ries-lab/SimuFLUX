import numpy as np

from .psf import Psf

class PsfGauss2D(Psf):
    def __init__(self):
        super().__init__()

        self.sigma = 310/2.35

    def intensity(self, flpos, patternpos, phasepattern=None, L=None):
        flposrel=np.atleast_2d(flpos-patternpos)
        r2=flposrel[:,:2]**2
        sigma=self.sigma
        rs2=np.sum(r2,1)
        
        io = np.exp(-rs2/2/sigma**2)
        if self.sigmaz>0:
            phfac=self.pinholezfac(flposrel)
            io=io*phfac
        else:
            phfac=np.ones(io.shape)

        return io, phfac

    def calculatePSFs(self, phasepattern=None, Lxs=None):
        pass