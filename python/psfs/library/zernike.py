# Jonas Ries: zernike based on Seidel function
# 
#Aberrated wavefront (Seidel coefficients).
#
# sys    System data
#  .Da      Aperture diameter [m]
#  .Zr      Zernike polynomials; [n,l,amplitude;...] n: radial (0...) l:
#  azimuthal. In units of Pi. Or better in rad?
#  .nm      Refraction index of the medium
#  .lo      Wavelength in free space [m]
#
# E      Electric field [V/m]
# r      Radial position [Da/2]
# t      Incidence angle [rad]
# p      Polar angle [rad]

#Copyright � Marcel Leutenegger, 2003-2007, �cole Polytechnique F�d�rale de Lausanne (EPFL),
#Laboratoire d'Optique Biom�dicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
#
#    This library is free software; you can redistribute it and/or modify it under
#    the terms of the GNU Lesser General Public License as published by the Free
#    Software Foundation; version 2.1 of the License.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT ANY
#    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License along
#    with this library; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

import numpy as np

from .zernfun import zernfun

def zernike(sys, E, r, t, p):
    Zr = sys['Zr']

    r = np.minimum(r, 1)  # Calculated only for r < 1, then extended
    e = np.zeros_like(r)

    for k in range(Zr.shape[0]):
        zpol = zernfun(Zr[k, 0], Zr[k, 1], r, p)
        zpolr = np.reshape(zpol, r.shape)
        e += zpolr * Zr[k, 2]

    # On axis aberrations
    # if 'defoc' in W:
    #    e += W['defoc'] * t  # r^2

    #E=E.*repmat(exp(2i*pi*e),size(E,1),1);  %B = repmat(A,m,n) returns an array containing m x n copies of A in the row and column dimensions.
    try:
        E *= np.exp(2j * np.pi * e)  # Multiply E with the phase factor
    except ValueError:
        # scale up to 3D
        E = E[:,None,None] * np.exp(2j * np.pi * e)[None,...]

    return E
