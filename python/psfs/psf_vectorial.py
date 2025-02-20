"""
Copyright (c) 2022      Ries Lab, EMBL, Heidelberg, Germany
All rights reserved     

@author: Sheng Liu, Jonas Hellgoth
"""
# import time
from dataclasses import dataclass, field

# import tensorflow as tf
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .psf import Psf, Parameters
from ..tools import utilities as im

@dataclass
class System:
    Zr: np.typing.ArrayLike = field(default_factory=lambda: np.zeros(2,3))
    maskshift: list = None
class PsfVectorial(Psf):
    """
    PSF class that uses a 3D volume to describe the PSF.
    Should only be used with single-channel data.
    """
    def __init__(self, parameters=None) -> None:
        super().__init__()
        if parameters is None:
            self.parameters = Parameters()
        else:
            self.parameters = parameters
        self.Zphase = None
        self.zT = None
        self.options = self.parameters.option
        self.initpupil = None
        self.defocus = np.float32(0)
        self.psftype = 'vector'
        self.psfs = {}
        self.psfph = None
    
    def intensity(self, flpos, patternpos, phasepattern, L=None):
        key = str(phasepattern)
        if L is not None:
            key = f"{phasepattern}{L}"

        try:
            psfint = self.psfs[key]
        except KeyError:
            psfint = self.calculate_psf(key, phasepattern=phasepattern, Lxs=L)

        phkey=self.psfph
        flposrel=flpos-patternpos
        # iexc=psfint.interp(flposrel) + self.zerooffset
        iexc=psfint(flposrel) + self.zerooffset
        iexc=np.maximum(iexc,0)
        try:
            psfph=self.psfs[phkey]
            # phfac=psfph.interp(flpos)
            phfac=psfph(flpos)
        except KeyError:
            phfac=np.ones(iexc.shape)
        idet = iexc * phfac

        return idet, phfac

    def genpsfmodel(self,sigma,pupil,addbead=False):
        phiz = np.exp(-1j*2*np.pi*self.kz*(self.Zrange+self.defocus))
        if self.psftype == 'vector':
            I_res = np.zeros((phiz.shape[0], self.paramxy[2].shape[0], self.paramxy[2].shape[1]), dtype=np.float32)
            PupilFunction = phiz[None,...]*pupil[None,None,...]*self.dipole_field[:,None,...]
            for h in range(self.dipole_field.shape[0]):
                # start = time.time()
                psfA = im.cztfunc1(PupilFunction[h,...],self.paramxy)
                # stop = time.time()
                # duration = stop-start
                # print(f"psfA calculation took {duration:.2f} s. Dtype is {psfA.dtype}")
                I_res += (psfA*np.conj(psfA)).real*self.normf
            # print(f"I_res shape {I_res.shape} dtype: {I_res.dtype}")
        else:
            PupilFunction = pupil*phiz
            I_res = im.cztfunc1(PupilFunction,self.paramxy)      
            I_res = (I_res*np.conj(I_res)).real*self.normf

        # I_res = I_res/np.max(I_res)#*self.normf

        # filter2 = np.exp(-2*sigma[1]*sigma[1]*self.kspace_x-2*sigma[0]*sigma[0]*self.kspace_y)
        # filter2 = filter2/np.max(filter2) + 0.0j
        # if addbead:
        #     I_blur = np.real(im.ifft3d(im.fft3d(I_res)*filter2*self.bead_kernel))
        # else:
        #     I_blur = np.real(im.ifft3d(im.fft3d(I_res)*filter2))
        # I_blur = np.expand_dims(I_blur,axis=-1)
        
        # bin = self.options.model.bin
        # kernel = np.ones((bin,bin,1,1),dtype=np.float32)
        # # kernel = np.ones((1,bin,bin,1),dtype=np.float32)
        # # I_modeltf = tf.nn.convolution(I_blur,kernel,strides=(1,bin,bin,1),padding='SAME',data_format='NHWC')
        # I_model = im.convolve(I_blur, kernel, bin)

        # # np.testing.assert_allclose(I_modeltf, I_model)

        # return I_model[...,0]
        return I_res
    
    def calculate_psf(self, key, phasepattern=None, Lxs=None, force=False):
        if not force:
            try:
                return self.psfs[key]
            except KeyError:
                pass
        
        # Apply phasepattern to pupil
        if phasepattern == 'vortex':
            self.calpupilfield(fieldtype=self.psftype,mask="phaseramp")
                            #    polarization="circular",
            # self.phase_ramp()
            # self.circular_polarization()
        else:
            self.calpupilfield(fieldtype=self.psftype)

        # Calculate PSF
        sigma = np.ones((len(self.parameters.roi.roi_size),))*self.options.model.blur_sigma*np.pi
        I_model = self.genpsfmodel(sigma, self.pupil)
        I_model = I_model/np.max(I_model[I_model.shape[0]//2,...])  # Normalize
        xy = (np.arange(I_model.shape[-1])-I_model.shape[-1]//2)*self.parameters.pixelsize_x*1000  # nm
        z = (np.arange(I_model.shape[0])-I_model.shape[0]//2)*self.parameters.pixelsize_z*1000     # nm
        # X, Y, Z = np.meshgrid(xy,xy,z)
        # self.psfs[key] = GriddedInterpolant(xy,xy,z,I_model)
        # print("Interp shapes", xy.shape, z.shape, I_model.shape)
        self.psfs[key] = RegularGridInterpolator((xy,xy,z), I_model.T, method='cubic')

        return self.psfs[key]

    def imagestack(self, key):
        try:
            psf = self.psfs[key]
        except KeyError:
            psf = self.calculate_psf(key)

        return psf.values

    def normfactor(self, phasemask, number):
        key =f"{phasemask}{number}"
        psf = self.calculate_psf(key, phasepattern=phasemask, Lxs=number)
        return psf.values[...,psf.values.shape[-1]//2].max()
