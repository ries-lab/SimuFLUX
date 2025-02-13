import numpy as np
from scipy import interpolate

class GriddedInterpolant:
    """Helper class designed to mimic https://www.mathworks.com/help/matlab/ref/griddedinterpolant.html """
    def __init__(self, X, Y, Z, PSF, intmethod='cubic', extraolation_method=None):
        self.points = (Z, Y, X,)
        self.psf = PSF
        self.intmethod = intmethod

        # Kind of useless. For interpn, extrapolation and interpolation methods are the same.
        self.extraolation_method = extraolation_method

    @property
    def values(self):
        """Interp evaluated on default grid points"""
        return self.psf

    def interp(self, xi):
        """ Expect xi as (x, y, z), but psf is stored as ZYX. """
        # fill_value = None guarantees extrapolation
        return interpolate.interpn(self.points, 
                                   self.psf, 
                                   xi[::-1], 
                                   method=self.intmethod, 
                                   bounds_error=False, 
                                   fill_value=None)
