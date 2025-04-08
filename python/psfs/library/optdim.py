#Optimal FFT dimension.
#
#Synopsis:
#  n=optdim(m)
#     Optimize the FFT performance by selecting the optimal dimension,
#     either an integer power of two or trice an integer power of two.

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

def optdim(m):
    n = np.ceil(np.log2(m) - 1)
    np15 = 1.5*2**n
    if np15 >= m:
        n = np15
    else:
        n = 2**(n+1)
    return int(n)

