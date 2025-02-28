
# import time
import os
import pathlib
from dataclasses import dataclass

# import tensorflow as tf
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .psf import Psf
from .library import effField
from .library import effIntensity
from .library import conv2fft
from .library import convnfft

@dataclass
class PSFStruct:
    normalization : float = 1
    interp : callable = None

def addzernikeaberrations(sys, addpar):
    if 'Zr' in sys:
        sys['Ei'].append('zernike')
    elif 'Zr' in addpar:
        sys['Ei'].append('zernike')
        sys['Zr'] = addpar['Zr']
    return sys

def normpsf(psf):
    return psf, np.max(psf)

def range4PSF(psf,dr,dz):
    nx = (np.arange(psf.shape[0]) - (psf.shape[0]-1)/2)*dr*1e9
    nz = (np.arange(psf.shape[2]) - (psf.shape[2]-1)/2)*dz*1e9
    
    return nx, nx, nz

def meshgrid4PSF(psf,dr,dz):
    nx, nz = range4PSF(psf,dr,dz)
    X, Y, Z = np.meshgrid(nx, nx, nz)
    
    return X, Y, Z

def genkey(phasepattern, Lxs=None):
    if Lxs is None:
        return f"{phasepattern}"
    elif isinstance(Lxs, (int, float)):
        Lx = [[Lxs]]
    else:
        # list or array (we hope)
        Lx = list(Lxs.copy())
        for i, x in enumerate(Lxs):
            if not isinstance(x , list):
                Lx[i] = [x]

    # unlist
    return f"{phasepattern}" + "".join([str(x) for xs in Lx for x in xs])

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
        # # overwrite any of the parameter fields if they exist
        # intersection_keys = list(self.parameters.keys() & kwargs.keys())

        # for key in intersection_keys:
        #     self.parameters[key].update(kwargs[key])

        # # otherwise add to addpar
        # exclusive_keys = list(set(kwargs.keys()) - set(self.parameters.keys()))
        # if len(exclusive_keys) > 0:
        #     if 'addpar' not in self.parameters.keys():
        #         self.parameters['addpar'] = {}
        #     for key in exclusive_keys:
        #         self.parameters['addpar'][key] = kwargs[key]

        # if len(intersection_keys) > 0:
        #     # recalculate
        #     self.PSFs = {}
        #     self.normfactgauss = None
        #     if self.pinholepar is not None:
        #         self.setpinhole(**self.pinholepar)
        fn = kwargs.keys()
        fnpar = self.parameters.keys()

        recalculate = False
        for k in fn:
            assigned = False
            for j in fnpar:
                if k in self.parameters[j].keys():
                    try:
                        if self.parameters[j][k] != kwargs[k]:
                            self.parameters[j][k] = kwargs[k]
                            recalculate = True
                    except ValueError:
                        # probably an array
                        if np.any(self.parameters[j][k] != kwargs[k]):
                            self.parameters[j][k] = kwargs[k]
                            recalculate = True
                    assigned = True
            if not assigned:
                self.parameters['addpar'][k] = kwargs[k]

        if recalculate:
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
            phfac = np.ones_like(iexc)
        idet = iexc * phfac

        return idet, phfac
    
    def calculatePSFs(self, phasepattern, Lxs=None, forcecalculation=False):
        key = genkey(phasepattern, Lxs)

        if key in self.PSFs and not forcecalculation:
            return key
        
        if Lxs is None:
            # safety, but not so safe so alert the user
            Lxs = [0]
            print("Warning! Unknown Lxs. Using a default of [0].")
        
        intmethod = 'cubic' # linear,cubic?
        bounds_error = False
        extraolation_method = None  # If None, extrapolate
        
        print(f"{key}, ", end="")

        opt = self.parameters['opt']
        out = self.parameters['addpar']
        sys = self.parameters['sys']

        fwhm = 0.51*sys['lo']/sys['NA']*1e9
        self.sigma = fwhm/2.35

        if self.normfactgauss is None:
            sys['Ei'] = ['circular']
            outc = effField(sys, out, opt)
            outc = effIntensity(sys, outc)
            zmid=int(np.ceil(outc['I'].shape[2]/2))
            self.normfactgauss=np.max(outc['I'][:,:,zmid])  # normalized to max =1

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
        elif "halfmoon" in phasepattern:
            # convert L into phase
            dxdphi = 1.4  # nm/degree, from LSA paper Deguchi
            # dzdphi = -3.6  #nm/degree
            sys['delshift'] = np.deg2rad(float(Lxs)/dxdphi)
            sys['Ei'] = ['halfmoon', 'linear']
            sys = addzernikeaberrations(sys, out)
            out = effField(sys, out, opt)
            out = effIntensity(sys, out)
            PSF = out['I'].squeeze()/self.normfactgauss
            PSF = self.beadsize(PSF, sys['beadradius'])
            PSFx = PSFStruct()
            PSFy = PSFStruct()
            PSF, PSFx.normalization = normpsf(PSF)
            xx, yy, zz = range4PSF(PSF,out['dr'],out['dz'])
            PSFx.interp = RegularGridInterpolator((xx, yy, zz), 
                                                   PSF.transpose(1,0,2), 
                                                   method = intmethod, 
                                                   bounds_error = bounds_error,
                                                   fill_value = extraolation_method)
            PSFy.interp = RegularGridInterpolator((xx, yy, zz), 
                                                   PSF, 
                                                   method = intmethod, 
                                                   bounds_error = bounds_error,
                                                   fill_value = extraolation_method)   
            PSFy.normalization = PSFx.normalization
            self.PSFs[genkey("halfmoonx",Lxs)] = PSFx
            self.PSFs[genkey("halfmoony",Lxs)] = PSFy
        elif phasepattern == "pinhole":
            sys['Ei'] = ['circular']
            sys = addzernikeaberrations(sys, out)
            phdiameter = float(Lxs[0])
            if len(Lxs) == 3:
                phpos = [float(Lxs[1]), float(Lxs[2])]
            else:
                phpos = [0, 0]
            out = effField(sys, out, opt)
            out = effIntensity(sys, out)
            psfg = out['I'].squeeze()
            norm = psfg[:,:,int(round(psfg.shape[2]/2))].sum()
            psfg /= norm
            pixelsize = opt['pixSize']*1e9
            nx = np.arange(psfg.shape[0])
            nx = (nx-np.mean(nx))*pixelsize
            Xk, Yk = np.meshgrid(nx, nx) 
            kernel = np.double(((Xk - phpos[1])**2 + (Yk - phpos[0])**2) < (phdiameter / 2)**2)
            # print(kernel.shape, psfg.shape)
            psfph = conv2fft(psfg, kernel)
            xx, yy, zz = range4PSF(psfph, out['dr'], out['dz'])
            PSFdonut = PSFStruct()
            PSFdonut.interp = RegularGridInterpolator((xx, yy, zz),
                                                      psfph, 
                                                      method = intmethod, 
                                                      bounds_error = bounds_error,
                                                      fill_value = extraolation_method)
            self.PSFs[key] = PSFdonut
        elif phasepattern == "tophat":
            dzdphi = -3.6  #nm/degree
            sys['delshift'] = np.deg2rad(float(Lxs)/dzdphi)
            sys['Ei'] = ['pishift', 'circular']
            sys = addzernikeaberrations(sys, out)
            out = effField(sys, out, opt)
            out = effIntensity(sys, out)
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
            xx, yy, zz = range4PSF(PSF, out['dr'], out['dz'])
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
        
        nx = np.arange(psf.shape[0])
        nx = (nx - np.mean(nx)) * dr
        
        ny = np.arange(psf.shape[1])
        ny = (ny - np.mean(ny)) * dr
        
        nz = np.arange(psf.shape[2])
        nz = (nz - np.mean(nz)) * dz
        
        X, Y, Z = np.meshgrid(nx, ny, nz)
        circpsf = np.double((X**2 + Y**2 + Z**2) <= R**2)
        circpsf /= np.sum(circpsf)
        
        psfo = convnfft(psf, circpsf)
        
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
            # print(f"calculated diameter is: {diameter}")

        phdefined = self.PSFph is not None and self.pinholepar is not None
        notchanged = phdefined and (diameter == self.pinholepar['diameter']) and offset == self.pinholepar['offset']
        calculated = notchanged and genkey('pinhole', [diameter,offset]) in self.PSFs.keys()

        if not calculated:
            self.PSFph = self.calculatePSFs('pinhole',
                                            Lxs = [diameter,offset],
                                            forcecalculation=True)

        self.pinholepar = {}
        self.pinholepar['offset'] = offset
        self.pinholepar['diameter'] = diameter

    def imagestack(self, key, show=False):
        if key not in self.PSFs.keys():
            key = self.calculatePSFs(key, 0)
        psf = self.PSFs[key]
        vout = psf.interp.values

        if self.PSFph is not None:
            # Apply pinhole
            phpsf = self.PSFs[self.PSFph]
            vout = vout*phpsf.interp.values

        return vout.transpose(1,0,2)
    
    def normfactor(self, phasemask, number):
        return self.PSFs[f"{phasemask}{number}"].normalization