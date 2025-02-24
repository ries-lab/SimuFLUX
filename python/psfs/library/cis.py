#z=cis(p,r)
#----------
#
#Sine and cosine.
#
#Input:
# p     Phase [rad]
# r     Magnitude {1}
#
#Output:
# z     r.*complex(cos(p),sin(p))
#
#Alternative:
# [c,s] To get the sine and cosine separately.

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

def cis(p=None, r=None, ret_both=False):
    if p is None:
        print("\nSine and cosine.\n\n\tMarcel Leutenegger © 1.5.2005\n")
        c, s = None, None
        return c, s
    elif r is None:
        c = np.exp(1j * p)
    elif p is not None and r is not None:
        c = r * np.exp(1j * p)
    else:
        raise ValueError("Incorrect number of arguments.")

    if ret_both:
        s = np.imag(c)
        c = np.real(c)
        return c, s
    return c
