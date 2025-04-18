from abc import ABCMeta, abstractmethod
import pickle

import numpy as np

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
    def intensity(self, flpos, patternpos, phasepattern, L):
        """
        Calculates PSF intensity and pinhole factor.

        Parameters
        ----------
        flpos : np.typing.ArrayLike
            with respect to optical axis ([0,0] of PSF coordinate system)
        patternpos: np.typing.ArrayLike
            position by EOD only of excitation pattern, no descanning. Pinhole is with
            respect to [0,0] (but can be shifted in pinhole definition).
        phasepattern: str
            Type of PSF, e.g. vortex
        L : float
            Displacement of PSF from 0
        """
        raise NotImplementedError("You need to implement a 'intensity' method in your psf class.")

    @abstractmethod
    def calculatePSFs(self, phasepattern, Lxs):
        """
        Calculates the forward images.
        """
        raise NotImplementedError("You need to implement a 'calculatePSFs' method in your psf class.")

    def savePSF(self, filename: str) -> None:
        """
        Save object to file.
        """
        with open(filename, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def loadPSF(filename: str):
        """
        Load object from file.
        """
        with open(filename, "rb") as f:
            self = pickle.load(f)
        return self


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
    