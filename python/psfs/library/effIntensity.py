#out=effIntensity(sys,out)
#-------------------------
#
#Electromagnetic intensity in the focus region.
#
#Input:
# sys    System data
#  .nm      Refraction index of the immersion
#  .ns      Refraction index of the sample
# out    Result data
#  .E       Electric field [V/m]
#
#Output:
# out    Result data
#  .I       Intensity [W/m^2]

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

C0 = 299792458
U0 = 4e-7*np.pi
E0 = 1/(U0*C0*C0)
SQE0U0 = np.sqrt(E0 / U0)

def effIntensity(sys, out):
    if 'ns' in sys and sys['ns'] < sys['nm']:
        n = sys['ns']
    else:
        n = sys['nm']
    
    out['I'] = SQE0U0 * np.real(n) / 2 * np.sum(np.real(np.conj(out['E'])*out['E']),axis=0)
    
    return out