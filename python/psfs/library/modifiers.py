#Execute field modifiers.
#
# sys    System parameters
# M      Samples on aperture radius
# N      Samples on lateral wavevector
#
# E      Electric field [V/m]
# r      Radial position [Da/2]
# t      Incidence angle [rad]
# p      Polar angle [rad]
# ct     cos(t)
# st     sin(t)

#Copyright (c) Marcel Leutenegger, 2003-2007, Ecole Polytechnique Federale de Lausanne (EPFL),
#Laboratoire d'Optique Biomedicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
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

# Functions we may need in locals()
from .phaseramp import phaseramp
from .linear import linear
from .circular import circular
from .zernike import zernike
from .pishift import pishift
from .halfmoon import halfmoon

def modifiers(sys, M, N):
    x = sys['NA'] / sys['nm'] / M * np.arange(-int(np.ceil(M+2)), int(np.ceil(M+2))+1)
    y = x.reshape(1, x.size)
    p = np.arctan2(np.tile(y, (x.size, 1)), np.tile(x, (y.size, 1)).T)
    r = np.tile(x**2, (y.size, 1)) + np.tile(y**2, (x.size, 1)).T
    ct = np.sqrt(1 - r, dtype='complex64')
    st = np.sqrt(r)
    t = np.arcsin(st, dtype='complex64')
    r = sys['nm'] / sys['NA'] * st  # radial position in aperture
    E = -0.5j * (sys['Da'] / 2 / M)**2 * sys['Mt'] / (sys['ft'] * sys['lo'])

    try:
        globals()['strat'] = None
    except KeyError:
        pass
    
    for func in sys['Ei']:
        E = globals()[func](sys, E, r, t, p)
        # print(func, E.shape)


    # E = E * np.tile((1 + np.tanh(2 * len(r) * (1 - r))) / ct, (E.shape[0], 1))
    if len(E.shape) > 1:
        E *= ((1 + np.tanh(2 * len(r) * (1 - r))) / ct)
    else:
        E = E[:,None, None]*((1 + np.tanh(2 * len(r) * (1 - r))) / ct)

    return E, r, t, p, ct, st

