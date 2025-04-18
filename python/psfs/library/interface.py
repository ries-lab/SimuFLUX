#Fresnel reflection and transmission coefficients at an interface
#for the electric field (p- and s-polarisation).
#
#Input:
# e1     Relative dielectric constants
# e2
# g0     Incidence sqrt(kx^2 + ky^2)/ko
#
#Output:
# tp     Transmission coefficients
# ts
# b2     Propagation kz/ko
# rp
# rs     Reflection coefficients

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

def interface(e1, e2, g0):
    
    b1 = np.sqrt(e1 - g0 * g0, dtype='complex64')
    b2 = np.sqrt(e2 - g0 * g0, dtype='complex64')
    rp = (e2 * b1 - e1 * b2) / (e2 * b1 + e1 * b2)
    tp = np.sqrt(e1 / e2, dtype='complex64') * (1 + rp)
    rs = (b1 - b2) / (b1 + b2)
    ts = 1 + rs

    return tp, ts, b2, rp, rs
