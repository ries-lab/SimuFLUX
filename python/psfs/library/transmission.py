#Transmission through an air-lens-air-immersion objective.
#
#The input is a plane wave incident on a lens focusing the
#beam in the medium. The refraction is balanced on the three
#interfaces air-lens, lens-air and air-medium.
#
# sys    System data
#  .nl      Refraction index of the lens
#  .nm      Refraction index of the medium
#
# t      Total deflection angle [rad]
# tp     Fresnel transmission coefficients
# ts

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

def transmission(sys, t):
    tp = t * (sys['nm'] - 1) / (2 * sys['nl'] + sys['nm'] - 3)
    ts = (t - tp) / 2
    
    #
    # Propagation constants
    #
    b = np.arctan(np.sin(ts) / (sys['nl'] - np.cos(ts)))
    a = np.cos(ts + b)
    b = np.cos(b)
    
    if sys['nm'] > 1:
        d = np.arctan(np.sin(tp) / (sys['nm'] - np.cos(tp)))
        c = np.cos(tp + d)
        d = np.cos(d)
    else:
        c = 1
        d = 1
    
    #
    # Fresnel coefficients
    #
    ts = (a - sys['nl'] * b) / (a + sys['nl'] * b)
    tp = (sys['nl'] * a - b) / (sys['nl'] * a + b)
    ts = (1 - ts**2) * 2 * c / (c + sys['nm'] * d)
    tp = (1 - tp**2) * 2 * c / (sys['nm'] * c + d)
    
    return tp, ts
