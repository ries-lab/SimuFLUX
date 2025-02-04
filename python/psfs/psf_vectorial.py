"""
Copyright (c) 2022      Ries Lab, EMBL, Heidelberg, Germany
All rights reserved     

@author: Sheng Liu, Jonas Hellgoth
"""

import numpy as np
from .psf import Psf
from ..tools import utilities as im

class PsfVectorial(Psf):
    """
    PSF class that uses a 3D volume to describe the PSF.
    Should only be used with single-channel data.
    """
    def __init__(self,options=None) -> None:
        self.parameters = None
        self.data = None
        self.Zphase = None
        self.zT = None
        self.bead_kernel = None
        self.options = options
        self.initpupil = None
        self.defocus = np.float32(0)
        self.psftype = 'vector'
    
    def intensity(self, flpos ,patternpos, phasepattern, L):
        pass

    def calc_forward_images(self, variables):
        """
        Calculate forward images from the current guess of the variables.
        Shifting is done by Fourier transform and applying a phase ramp.
        """

        pos, backgrounds, intensities, pupilR, pupilI, sigma, gxy = variables

        if self.const_mag:
            pupil_mag = np.complex128(1.0 + 0.0j)
        else:
            pupil_mag = np.complex128(pupilR*self.weight[4] + 0.0j)

        if self.initpupil is not None:
            pupil = self.initpupil
            normp = np.complex128(1.0 + 0.0j)
        else:
            pupil_phase = (np.cos(pupilI * self.weight[3]) + 1j * np.sin(pupilI * self.weight[3])) * self.aperture
            pupil_phase0 = (np.cos(pupilI * 0.0) + 1j * np.sin(pupilI * 0.0)) * self.aperture
            normp = self.calnorm(pupil_phase)/self.calnorm(pupil_phase0)
            pupil = pupil_phase*pupil_mag*self.apoid
        
        self.psfnorm = normp

        Nz = self.Zrange.shape[0]
        pos = np.array(np.reshape(pos, pos.shape + (1,1,1)), dtype=np.complex128)
        phiz = -1j*2*np.pi*self.kz*(pos[:,0]+self.Zrange+self.defocus)
        if pos.shape[1]>3:
            phixy = 1j*2*np.pi*self.ky*pos[:,2]+1j*2*np.pi*self.kx*pos[:,3]
            phiz = 1j*2*np.pi*(self.kz_med*pos[:,1]-self.kz*(pos[:,0]+self.Zrange))
        else:
            phixy = 1j*2*np.pi*self.ky*pos[:,1]+1j*2*np.pi*self.kx*pos[:,2]

        
        if self.psftype == 'vector':
            I_res = 0.0
            for h in self.dipole_field:
                PupilFunction = pupil * np.exp(phiz + phixy) * h
                psfA = im.cztfunc1(PupilFunction,self.paramxy)     
                I_res += psfA*np.conj(psfA)*self.normf
        else:
            PupilFunction = pupil*np.exp(phiz+phixy)
            I_res = im.cztfunc1(PupilFunction,self.paramxy)
            I_res = I_res*np.conj(I_res)*self.normf

        bin = self.options.model.bin
        if not self.options.model.var_blur:
            sigma = np.ones((2,))*self.options.model.blur_sigma*np.pi
        filter2 = np.exp(-2*sigma[1]*sigma[1]*self.kspace_x-2*sigma[0]*sigma[0]*self.kspace_y)
        filter2 = filter2/np.max(filter2) + 0.0j
        I_blur = im.ifft3d(im.fft3d(I_res)*self.bead_kernel*filter2)
        #I_blur = im.ifft3d(im.fft3d(I_res)*filter2)
        I_blur = np.expand_dims(np.real(I_blur),axis=-1)
        
        kernel = np.ones((1,bin,bin,1,1),dtype=np.float32)
        I_blur_bin = tf.nn.convolution(I_blur,kernel,strides=(1,1,bin,bin,1),padding='SAME',data_format='NDHWC')

        psf_fit = I_blur_bin[...,0]
        
        st = (self.bead_kernel.shape[0]-self.data.rois[0].shape[-3])//2
        psf_fit = psf_fit[:,st:Nz-st]

        if self.options.model.estimate_drift:
            gxy = gxy*self.weight[2]
            psf_shift = self.applyDrfit(psf_fit,gxy)
            forward_images = psf_shift*intensities*self.weight[0] + backgrounds*self.weight[1]
        else:
            forward_images = psf_fit*intensities*self.weight[0] + backgrounds*self.weight[1]

        return forward_images

    def genpsfmodel(self,sigma,pupil,addbead=False):
        phiz = -1j*2*np.pi*self.kz*(self.Zrange+self.defocus)
        if self.psftype == 'vector':
            I_res = 0.0
            for h in self.dipole_field:
                PupilFunction = pupil*np.exp(phiz)*h
                psfA = im.cztfunc1(PupilFunction,self.paramxy)      
                I_res += psfA*np.conj(psfA)*self.normf
        else:
            PupilFunction = pupil*np.exp(phiz)
            I_res = im.cztfunc1(PupilFunction,self.paramxy)      
            I_res = I_res*np.conj(I_res)*self.normf
    
        bin = self.options.model.bin

        filter2 = np.exp(-2*sigma[1]*sigma[1]*self.kspace_x-2*sigma[0]*sigma[0]*self.kspace_y)
        filter2 = filter2/np.max(filter2) + 0.0j
        if addbead:
            I_blur = np.real(im.ifft3d(im.fft3d(I_res)*filter2*self.bead_kernel))
        else:
            I_blur = np.real(im.ifft3d(im.fft3d(I_res)*filter2))
        I_blur = np.expand_dims(np.real(I_blur),axis=-1)
        
        kernel = np.ones((bin,bin,1,1),dtype=np.float32)
        I_model = tf.nn.convolution(I_blur,kernel,strides=(1,bin,bin,1),padding='SAME',data_format='NHWC')
        I_model = I_model[...,0]

        return I_model

    def postprocess(self, variables):
        """
        Applies postprocessing to the optimized variables. In this case calculates
        real positions in the image from the positions in the roi. Also, normalizes
        psf and adapts intensities and background accordingly.
        """
        positions, backgrounds, intensities, pupilR,pupilI,sigma,gxy = variables
        z_center = (self.Zrange.shape[-3] - 1) // 2
        bin = self.options.model.bin
        positions[:,1:] = positions[:,1:]/bin

        pupil_mag = pupilR*self.weight[4] + 0.0j
        if self.initpupil is not None:
            pupil = self.initpupil
        else:
            pupil = (np.cos(pupilI*self.weight[3]) + 1j * np.sin(pupilI*self.weight[3]))*pupil_mag*self.aperture*self.apoid

        I_model = self.genpsfmodel(sigma,pupil)
        I_model_bead = self.genpsfmodel(sigma,pupil,addbead=True)
        #I_model_bead = np.real(im.ifft3d(im.fft3d(I_res)*self.bead_kernel*filter2))

        images, _, centers, _ = self.data.get_image_data()
        if positions.shape[1]>3:
            global_positions = np.swapaxes(np.vstack((positions[:,0]+z_center,positions[:,1],centers[:,-2]-positions[:,-2],centers[:,-1]-positions[:,-1])),1,0)
        else:
            global_positions = np.swapaxes(np.vstack((positions[:,0]+z_center,centers[:,-2]-positions[:,-2],centers[:,-1]-positions[:,-1])),1,0)

        return [global_positions.astype(np.float32),
                backgrounds*self.weight[1], # already correct
                intensities*self.weight[0], # already correct
                I_model_bead,
                I_model,
                np.complex64(pupil),
                sigma,
                gxy*self.weight[2],
                np.flip(I_model,axis=-3),
                variables] # already correct


    def res2dict(self,res):
        res_dict = dict(pos=res[0],
                        bg=np.squeeze(res[1]),
                        intensity=np.squeeze(res[2]),
                        I_model_bead = res[3],
                        I_model = res[4],
                        pupil = res[5],
                        sigma = res[6]/np.pi,
                        drift_rate=res[7],
                        I_model_reverse = res[8],
                        offset=np.min(res[4]),
                        apodization = self.apoid,
                        cor_all = self.data.centers_all,
                        cor = self.data.centers)    
        return res_dict