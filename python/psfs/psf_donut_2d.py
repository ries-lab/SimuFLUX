import numpy as np

from .psf import Psf

class PsfDonut2D(Psf):
    def __init__(self):
        super().__init__()

        self.fwhm = 310  # comparison with calculated PSF

    def intensity(self, flpos, patternpos, phasepattern, L):
        flposrel=np.atleast_2d(flpos-patternpos)
        r2=flposrel[:,:2]**2
        fwhm=self.fwhm
        rs2=np.sum(r2,1)/fwhm**2
        zerooffset=self.zerooffset
        # io=2.2610*rs2*np.exp(-2.7726*rs2)+zerooffset
        io = 2.7726*rs2*np.exp(-2.7726*rs2) + zerooffset
        if self.sigmaz>0:
            phfac=self.pinholezfac(flposrel)
            io=io*phfac
        else:
            phfac=np.ones(io.shape)

        return io, phfac

    def calculatePSFs(self, phasepattern, Lxs):
        pass