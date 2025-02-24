
# import time
import os
import pathlib
from dataclasses import dataclass

# import tensorflow as tf
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import convolve

from .psf import Psf
from .library import effField, effIntensity

@dataclass
class PSFStruct:
    normalization : float = None
    interp : callable = None

def addzernikeaberrations(sys, addpar):
    if 'Zr' in sys and sys['Zr']:
        sys['Ei'].append('zernike')
    elif 'Zr' in addpar and addpar['Zr']:
        sys['Ei'].append('zernike')
        sys['Zr'] = addpar['Zr']
    return sys

def normpsf(psf):
    return psf, np.max(psf)

def range4PSF(psf,dr,dz):
    nx = (np.arange(psf.shape[0]) - (psf.shape[0]+1)/2)*dr*1e9
    nz = (np.arange(psf.shape[2]) - (psf.shape[2]+1)/2)*dz*1e9
    
    return nx, nx, nz

def meshgrid4PSF(psf,dr,dz):
    nx, nz = range4PSF(psf,dr,dz)
    X, Y, Z = np.meshgrid(nx, nx, nz)
    
    return X, Y, Z

def genkey(phasepattern, Lxs=None):
    if Lxs is None:
        return f"{phasepattern}"
    
    try:
        if isinstance(Lxs[0], (int, float)):
            Lx = Lxs
            Lxs = Lx.astype(str)
        else:
            Lx = Lxs.astype(float)
    except TypeError:
        if isinstance(Lxs, (int, float)):
            Lx = [Lxs]
            Lxs = [str(x) for x in Lx]
        else:
            Lx = [float(Lxs)]
    return f"{phasepattern}" + "".join([str(x) for x in Lxs])

class PsfVectorial(Psf):
    """
    PSF class that uses a 3D volume to describe the PSF.
    Should only be used with single-channel data.
    """
    def __init__(self, parameters=None) -> None:
        super().__init__()
        
        if parameters is None:
            import yaml

            base_path = pathlib.Path(__file__).parent.parent
            param_path = os.path.join(base_path, 'settings/default_microscope.yaml')

            with open(param_path, 'r') as fn:
                self.parameters = yaml.load(fn, Loader=yaml.FullLoader)
        else:
            self.parameters = parameters
        self.PSFs = {}
        self.PSFph = None
        self.normfactgauss = None
        self.pinholepar = None
        self.sigma = 310/2.355
        self.fwhm = 310

    def setpar(self, **kwargs):
        # overwrite any of the parameter fields if they exist
        intersection_keys = list(self.parameters.keys() & kwargs.keys())

        for key in intersection_keys:
            self.parameters[key].update(kwargs[key])

        # otherwise add to addpar
        exclusive_keys = list(set(kwargs.keys()) - set(self.parameters.keys()))
        if len(exclusive_keys) > 0:
            self.parameters['addpar'] = {}
            for key in exclusive_keys:
                self.parameters['addpar'][key] = kwargs[key]

        if len(intersection_keys) > 0:
            # recalculate
            self.PSFs = {}
            self.normfactgauss = None
            if self.pinholepar is not None:
                self.setpinhole(**self.pinholepar)

    def intensity(self, flpos, patternpos, phasepattern, L=None):
        key = genkey(phasepattern, L)

        try:
            psfint = self.PSFs[key]
        except KeyError:
            key = self.calculatePSFs(phasepattern, L)
            psfint = self.PSFs[key]

        flposrel = flpos - patternpos
        iexc = psfint.interp(flposrel) + self.zerooffset
        iexc = np.maximum(iexc,0)
        try:
            phfac = self.PSFs[self.PSFph].interp(flpos)
        except KeyError:
            phfac = np.ones(iexc.shape)
        idet = iexc * phfac

        return idet, phfac
    
    def calculatePSFs(self, phasepattern, Lxs=None, forcecalculation=False):
        # print(phasepattern)
        key = genkey(phasepattern, Lxs)

        if key in self.PSFs and not forcecalculation:
            return key
        
        intmethod = 'cubic' # linear,cubic?
        bounds_error = False
        extraolation_method = None  # If None, extrapolate
        
        print(key, ", ", end="")

        opt = self.parameters['opt']
        out = self.parameters['addpar']
        sys = self.parameters['sys']

        fwhm = 0.51*sys['loem']/sys['NA']*1e9
        self.sigma = fwhm/2.35

        if self.normfactgauss is None:
            sys['Ei'] = ['circular']
            outc = effField(sys, out, opt)
            outc = effIntensity(sys, outc)
            zmid=int(np.ceil(outc['I'].shape[2]/2))
            self.normfactgauss=np.max(np.max(outc['I'][...,zmid]))  # normalized to max =1

        if phasepattern == "flat":
            sys['Ei'] = ['circular']
            sys = addzernikeaberrations(sys,out)
            out = effField(sys, out, opt)
            out = effIntensity(sys,out)
            PSF = out['I'].squeeze()/self.normfactgauss
            PSF = self.beadsize(PSF, sys['beadradius'])
            PSFdonut = PSFStruct()
            PSF, PSFdonut.normalization = normpsf(PSF)
            xx, yy, zz = range4PSF(PSF,out['dr'],out['dz'])
            PSFdonut.interp = RegularGridInterpolator((xx, yy, zz), 
                                                       PSF, 
                                                       method = intmethod, 
                                                       bounds_error = bounds_error,
                                                       fill_value = extraolation_method)
            self.PSFs[key] = PSFdonut
        elif phasepattern == "vortex":
            sys['Ei'] = ['phaseramp', 'circular']
            sys = addzernikeaberrations(sys, out)
            out = effField(sys, out, opt)       
            out = effIntensity(sys, out)
            PSF = out['I'].squeeze()/self.normfactgauss
            PSF = self.beadsize(PSF, sys['beadradius'])
            PSFdonut = PSFStruct()
            PSF, PSFdonut.normalization = normpsf(PSF)
            xx, yy, zz = range4PSF(PSF,out['dr'],out['dz'])
            # print(xx.shape, yy.shape, zz.shape, PSF.shape, PSF.dtype, np.sum(np.isnan(PSF)))
            PSFdonut.interp = RegularGridInterpolator((xx, yy, zz), 
                                                      PSF, 
                                                      method = intmethod, 
                                                      bounds_error = bounds_error,
                                                      fill_value = extraolation_method)
            self.PSFs[key] = PSFdonut
        else:
            raise UserWarning(f"{phasepattern} PSF name not defined.")
        
        return key


    def beadsize(self, psf, R):
        if R == 0:
            return psf
        
        dr = self.parameters['addpar']['dr']
        dz = self.parameters['addpar']['dz']
        
        nx = np.arange(1, psf.shape[0] + 1)
        nx = (nx - np.mean(nx)) * dr
        
        ny = np.arange(1, psf.shape[1] + 1)
        ny = (ny - np.mean(ny)) * dr
        
        nz = np.arange(1, psf.shape[2] + 1)
        nz = (nz - np.mean(nz)) * dz
        
        X, Y, Z = np.meshgrid(nx, ny, nz, indexing='ij')
        circpsf = (X**2 + Y**2 + Z**2 <= R**2).astype(float)
        circpsf /= np.sum(circpsf)
        
        psfo = convolve(psf, circpsf, mode='constant')
        
        return psfo
    
    def setpinhole(self, lambda_nm=600, AU=1, diameter=None, offset=[0,0]):
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
        if len(offset) != 2:
            raise ValueError('PSF_vectorial.setpinhole: offset needs to have two entries [x,y]')
        # Compute diameter if not provided
        if diameter is None:
            diameter = np.round(AU*1.22*self.parameters['sys']['loem']*1e9/self.parameters['sys']['NA'])

        phdefined = self.PSFph is not None and self.pinholepar is not None
        notchanged = phdefined and (diameter == self.pinholepar['diameter']) and all(offset == self.pinholepar['offset'])
        calculated = notchanged and f"pinhole{diameter}{offset}" in self.PSFs.keys()

        if not calculated:
            self.PSFph = self.calculatePSFs('pinhole',[diameter,offset],1)

        self.pinholepar = {}
        self.pinholepar['offset'] = offset
        self.pinholepar['diameter'] = diameter

    def imagestack(self, key, show=False):
        key = self.calculatePSFs(key)
        psf = self.PSFs[key]
        vout = psf.interp.values

        if self.PSFph is not None:
            # Apply pinhole
            phpsf = self.PSFs[self.PSFph]
            vout = vout*phpsf.interp.values

        return vout.transpose(1,0,2)
    
    def normfactor(self, phasemask, number):
        return self.PSFs[f"{phasemask}{number}"].normalization