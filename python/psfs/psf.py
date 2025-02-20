"""
Adapation of original code by @kiwibogo:
https://github.com/ries-lab/uiPSF
"""

from abc import ABCMeta, abstractmethod
import pickle
from dataclasses import dataclass, field

import numpy as np
import scipy.special as spf

from ..tools import utilities as im
from ..tools.math import pol2cart, cart2pol
@dataclass
class ROI:
    roi_size: list = field(default_factory=lambda: [160,140,140])  # vector of 2 or 3 element, roi size in [y, x] or [z, y, x]
@dataclass
class Model:
    zernike_nl: list = None           # set the Zernike terms for PSF learning, e.g. [(2,2),(2,-2)], if empty, all zernike terms defined by n_max will be used 
    pupilsize: int = 35               # unit: pixel
    n_max: int =  8                   # maximum zernike order
    blur_sigma: float = 0             # unit: pixel
    with_apoid: bool = True           # with theoretical apoidization term
    bin: int = 1                      # upsamplling pixel size equal to camera pixel size divided by bin

@dataclass
class Ri:
    imm: float = 1.406
    med: float = 1.406
    cov: float = 1.512

@dataclass
class Imaging:
    RI: Ri = field(default_factory=lambda: Ri())
    emission_wavelength: float = 0.7  # unit: micron
    NA: float = 1.35

@dataclass
class Option:
    model: Model = field(default_factory=lambda: Model())
    imaging: Imaging = field(default_factory=lambda: Imaging())

@dataclass
class Parameters:
    roi: ROI = field(default_factory=lambda: ROI())
    option: Option = field(default_factory=lambda: Option())
    pixelsize_x: float = 0.01        # unit: micron
    pixelsize_y: float = 0.01
    pixelsize_z: float = 0.01
    bead_radius: float = 0.0          # unit: micron
    maskshift: list = None

class Psf():
    """
    Interface that ensures consistency and compatability between all old and new implementations of data classes, fitters and psfs.
    Classes implementing this interafce define a psf model/parametrization. They describe how the parameters of the psf are used to calculate a forward image
    at a specific position. They also provide initial values and postprocessing of the variables for the fitter,
    since they depend on the nature of the psf model/parametrization.
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        super().__init__()

        self.sigmaz = 0
        self.zerooffset = 0
        self._bead_kernel = None
        self.pupil = None
        self.phi = None
        self.kr = None

    @abstractmethod
    def intensity(self, flpos ,patternpos, phasepattern, L):
        """
        Calculates PSF intensity and pinhole factor.

        Parameters
        ----------
        flpos : np.typing.ArrayLike
            with respect to optical axis ([0,0] of PSF coordinate system)
        patternpos: np.typing.ArrayLike
            position by EOD only of excitation pattern, no descanning. Pinhole is with
            respect to [0,0] (but can be shifted in pinhole definition).
        """
        raise NotImplementedError("You need to implement a 'calc_forward_images' method in your psf class.")

    @abstractmethod
    def calc_forward_images(self, variables: list) -> np.typing.ArrayLike:
        """
        Calculates the forward images.
        """
        raise NotImplementedError("You need to implement a 'calc_forward_images' method in your psf class.")

    @abstractmethod
    def postprocess(self, variables: list) -> list:
        """
        Postprocesses the optimized variables. For example, normalizes the psf or calculates global positions.
        """
        raise NotImplementedError("You need to implement a 'postprocess' method in your psf class.")

    def save(self, filename: str) -> None:
        """
        Save object to file.
        """
        with open(filename, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load(filename: str):
        """
        Load object from file.
        """
        with open(filename, "rb") as f:
            self = pickle.load(f)
        return self
    
    @property
    def bead_kernel(self):
        if self._bead_kernel is None:
            self.gen_bead_kernel()
        return self._bead_kernel
    
    def gen_bead_kernel(self,isVolume = False):
        pixelsize_z = self.parameters.pixelsize_z
        bead_radius = self.parameters.bead_radius
        if isVolume:
            Nz = self.parameters.roi.roi_size[-3]
            bin = 1
        else:
            Nz = self.parameters.roi.roi_size[-3]+np.int32(bead_radius//pixelsize_z)*2+4
            bin = self.options.model.bin
        
        Lx = self.parameters.roi.roi_size[-1]*bin
        pixelsize_x = self.parameters.pixelsize_x/bin
        pixelsize_y = self.parameters.pixelsize_y/bin

        xrange = np.linspace(-Lx/2+0.5,Lx/2-0.5,Lx)+1e-6
        zrange = np.linspace(-Nz/2+0.5,Nz/2-0.5,Nz)
        [xx,yy,zz] = np.meshgrid(xrange,xrange,zrange)
        xx = np.swapaxes(xx,0,2)
        yy = np.swapaxes(yy,0,2)
        zz = np.swapaxes(zz,0,2)

        pkx = 1/Lx/pixelsize_x
        pky = 1/Lx/pixelsize_y
        pkz = 1/Nz/pixelsize_z
        if bead_radius>0:
            Zk0 = np.sqrt((xx*pkx)**2+(yy*pky)**2+(zz*pkz)**2)*bead_radius
            mu = 1.5
            kernel = spf.jv(mu,2*np.pi*Zk0  )/(Zk0**mu)*bead_radius**3
            kernel = kernel/np.max(kernel)
            kernel = np.float32(kernel)
        else:
            kernel = np.ones((Nz,Lx,Lx),dtype=np.float32)
        self._bead_kernel = kernel + 0j

        return 
    
    def calpupilfield(self,fieldtype='vector',Nz=None,datatype='bead',polarization=None, mask=None):
        if Nz is None:
            # Nz = self.bead_kernel.shape[0]
            pixelsize_z = self.parameters.pixelsize_z
            bead_radius = self.parameters.bead_radius
            Nz = Nz = self.parameters.roi.roi_size[-3]+np.int32(bead_radius//pixelsize_z)*2+4
        bin = self.options.model.bin
        Lx = self.parameters.roi.roi_size[-1]*bin 
        # Ly = self.parameters.roi.roi_size[-2]*bin
        # Lz = self.parameters.roi.roi_size[-3]
        # xsz = self.options.model.pupilsize
        xsz = self.parameters.roi.roi_size[-1]*bin 

        xrange = np.linspace(-Lx/2+0.5,Lx/2-0.5,Lx)
        [xx,yy] = np.meshgrid(xrange,xrange)
        pkx = xx/Lx
        pky = yy/Lx     
        self.kspace = np.float32(pkx*pkx+pky*pky)
        self.kspace_x = np.float32(pkx*pkx)
        self.kspace_y = np.float32(pky*pky)

        pixelsize_x = self.parameters.pixelsize_x/bin
        pixelsize_y = self.parameters.pixelsize_y/bin
        NA = self.options.imaging.NA
        emission_wavelength = self.options.imaging.emission_wavelength
        nimm = self.options.imaging.RI.imm
        nmed = self.options.imaging.RI.med
        ncov = self.options.imaging.RI.cov
        n_max = self.options.model.n_max
        # Zk = im.genZern1(n_max,xsz)

        n1 = np.array(range(-1,n_max,2))
        self.spherical_terms = (n1+1)*(n1+2)//2

        pupilradius = 1
        krange = np.linspace(-pupilradius+pupilradius/xsz,pupilradius-pupilradius/xsz,xsz)
        [xx,yy] = np.meshgrid(krange,krange)
        self.kr = np.lib.scimath.sqrt(xx**2+yy**2)
        kz = np.lib.scimath.sqrt((nimm/emission_wavelength)**2-(self.kr*NA/emission_wavelength)**2)

        cos_imm = np.lib.scimath.sqrt(1-(self.kr*NA/nimm)**2)
        cos_med = np.lib.scimath.sqrt(1-(self.kr*NA/nmed)**2)
        cos_cov = np.lib.scimath.sqrt(1-(self.kr*NA/ncov)**2)
        # kz_med = nmed/emission_wavelength*cos_med
        FresnelPmedcov = 2*nmed*cos_med/(nmed*cos_cov+ncov*cos_med)
        FresnelSmedcov = 2*nmed*cos_med/(nmed*cos_med+ncov*cos_cov)
        FresnelPcovimm = 2*ncov*cos_cov/(ncov*cos_imm+nimm*cos_cov)
        FresnelScovimm = 2*ncov*cos_cov/(ncov*cos_cov+nimm*cos_imm)
        Tp = FresnelPmedcov*FresnelPcovimm
        Ts = FresnelSmedcov*FresnelScovimm
        Tavg = (Tp+Ts)/2

        self.phi = np.arctan2(yy,xx)

        if mask == "phaseramp" and self.parameters.maskshift is not None:
            print(f"self.parameters.maskshift is {self.parameters.maskshift}")
            xshift, yshift = self.parameters.maskshift
            xx,yy = pol2cart(self.phi, self.kr)
            xx = xx + xshift
            yy = yy + yshift
            self.phi, self.kr = cart2pol(xx, yy)

        cos_phi = np.cos(self.phi).astype('complex64')
        sin_phi = np.sin(self.phi).astype('complex64')
        sin_med = (self.kr*NA/nmed).astype('complex64')

        px, py, pz = cos_med*cos_phi,cos_med*sin_phi,-sin_med
        sx, sy, sz = -sin_phi,cos_phi,np.zeros(cos_phi.shape).astype('complex64')

        if polarization == "circular":
            px = (px/np.sqrt(2) - py*1j/np.sqrt(2)).astype('complex64')
            py = (px*1j/np.sqrt(2) + py/np.sqrt(2)).astype('complex64')
            sx = (sx/np.sqrt(2) - sy*1j/np.sqrt(2)).astype('complex64')
            sy = (sx*1j/np.sqrt(2) + sy/np.sqrt(2)).astype('complex64')

        # if mask == "phaseramp":
        #     pr = np.exp(1j*self.phi).astype('complex64')
        #     px *= pr
        #     py *= pr
        #     pz *= pr
        #     sx *= pr
        #     sy *= pr
        #     sz *= pr

        pvec = Tp*np.stack([px,py,pz])
        svec = Ts*np.stack([sx,sy,sz])

        hx = cos_phi*pvec-sin_phi*svec
        hy = sin_phi*pvec+cos_phi*svec

        # if polarization == "circular":
        #     # hx = hx/np.sqrt(2) - hy*1j/np.sqrt(2)
        #     # hy = hx*1j/np.sqrt(2) + hy/np.sqrt(2)
        #     hx[0,...] = hx[0,...]/np.sqrt(2) - hx[1,...]*1j/np.sqrt(2)
        #     hx[1,...] = hx[0,...]*1j/np.sqrt(2) + hx[1,...]/np.sqrt(2)
        #     hy[0,...] = hy[0,...]/np.sqrt(2) - hy[1,...]*1j/np.sqrt(2)
        #     hy[1,...] = hy[0,...]*1j/np.sqrt(2) + hy[1,...]/np.sqrt(2)

        self.dipole_field = np.concatenate((hx,hy),axis=0,dtype='complex64')

        if self.options.model.with_apoid:
            #apoid = 1/np.lib.scimath.sqrt(cos_med)
            apoid = np.lib.scimath.sqrt(cos_imm)/cos_med
            #apoid = np.lib.scimath.sqrt(cos_med)/cos_imm
            if fieldtype=='scalar':
                apoid=apoid*Tavg
        else:
            apoid = 1

        kpixelsize = 2.0*NA/emission_wavelength/xsz
        self.paramxy = im.prechirpz1(kpixelsize,pixelsize_x,pixelsize_y,xsz,Lx)

        self.aperture = np.complex64(self.kr<1)
        pupil = self.aperture*apoid

        self.Zrange = np.linspace(-Nz/2+0.5,Nz/2-0.5,Nz,dtype=np.complex64).reshape((Nz,1,1))
        # self.kx = np.complex64(xx*NA/emission_wavelength)*pixelsize_x
        # self.ky = np.complex64(yy*NA/emission_wavelength)*pixelsize_y
        self.kz = np.complex64(kz)*self.parameters.pixelsize_z

        self.pupil = np.asarray(pupil, dtype=np.complex64)

        if mask == "phaseramp":
            self.pupil *= np.exp(1j*self.phi)

        self.normf = self.calnorm(self.pupil)
        
        # if fieldtype=='scalar':
        #     psfA = im.cztfunc1(self.pupil,self.paramxy)
        #     # self.normf = 1/np.sum((psfA*np.conj(psfA)).real)
        #     self.norm = 1/np.max((psfA*np.conj(psfA)).real)
        # else:
        #     I_res = np.zeros(self.paramxy[2].shape, dtype=np.float32)
        #     for h in self.dipole_field:
        #         PupilFunction = self.pupil*h
        #         psfA = im.cztfunc1(PupilFunction,self.paramxy)     
        #         I_res += (psfA*np.conj(psfA)).real
        #     # self.normf = 1/np.sum(I_res)
        #     self.normf = 1/np.max(I_res)
        # #if datatype == 'bead':
        # #    self.Zrange = -1*np.linspace(-Nz/2+0.5,Nz/2-0.5,Nz,dtype=np.complex64).reshape((Nz,1,1))
        # #elif datatype == 'insitu':
        # self.kz_med = np.complex64(kz_med)*self.parameters.pixelsize_z
        # self.k = np.complex64(nmed/emission_wavelength)*self.parameters.pixelsize_z
        # self.apoid = np.complex64(apoid)
        # self.nimm = nimm
        # self.nmed = nmed
        # self.Zk = np.float32(Zk)

        # # only for bead data, precompute phase ramp
        # Lx = self.parameters.roi.roi_size[-1]      
        # Ly = self.parameters.roi.roi_size[-2]
        # Lz = self.parameters.roi.roi_size[-3]

        # self.zv = np.linspace(0,Lz-1,Lz,dtype=np.float32).reshape(Lz,1,1)-Lz/2
        # self.kxv = np.linspace(-Lx/2+0.5,Lx/2-0.5,Lx,dtype=np.float32)/Lx
        # self.kyv = (np.linspace(-Ly/2+0.5,Ly/2-0.5,Ly,dtype=np.float32).reshape(Ly,1))/Ly
        # self.kzv = (np.linspace(-Lz/2+0.5,Lz/2-0.5,Lz,dtype=np.float32).reshape(Lz,1,1))/Lz

    def calnorm(self,pupil):
        psfA = im.cztfunc1(pupil,self.paramxy)   
        # normf = np.sum((psfA * np.conj(psfA)).real)
        normf = np.max((psfA * np.conj(psfA)).real)
        return normf

    def setpinholesimple(self, lambda_nm=600, AU=1, diameter=None, NA=1.5, refractive_index=1.5):
        """
        Sets the pinhole size based on given parameters.

        Parameters
        ----------
        lambda_nm (float): Wavelength in nm (default: 600 nm)
        AU (float): Airy unit scaling factor (default: 1)
        diameter (float or None): Pinhole diameter in nm (default: computed if None)
        NA (float): Numerical aperture (default: 1.5)
        refractive_index (float): Medium refractive index (default: 1.5)
        """
        # Compute diameter if not provided
        if diameter is None:
            diameter = AU * 1.22 * lambda_nm / NA

        n = refractive_index
        fwhm2 = (0.88 * lambda_nm / (n - np.sqrt(n**2 - NA**2)))**2 + (np.sqrt(2) * n * diameter / NA)**2
        self.sigmaz = np.sqrt(fwhm2) / 2.35

    def pinholezfac(self, flposrel):
        """
        Computes the pinhole z-factor based on fluorophore positions.

        Parameters
        ----------
        flposrel (numpy.ndarray): Array of fluorophore positions (N x 3).

        Returns
        -------
        numpy.ndarray: Computed pinhole z-factor.
        """
        sigmaz = self.sigmaz
        if sigmaz > 0:
            phfac = np.exp(-flposrel[:, 2]**2 / (2 * sigmaz**2))
        else:
            phfac = np.ones((flposrel.shape[0], 1))

        return phfac
    
    def normfactor(self, *args, **kwargs):
        return 1
    