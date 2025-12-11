#2D-Gaussian amplitude.
#
# sys    System data
#  .Da      Aperture diameter [m]
#  .Po      Laser power [W]
#  .wa      Beam waist at aperture [m]
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

c0 = 299792458
u0 = 4e-7*np.pi
e0 = 1/(u0*c0**2)

def gauss(sys,E,r,t,p):
    e=2*np.sqrt(np.sqrt(u0/e0)*sys['Po']/np.pi)/sys['wa']
    e=e*np.exp(-(sys['Da']/sys['wa']/2)**2*r**2)
    try:
        E = E * np.tile(e, (E.shape[0], 1))
    except AttributeError:
        E *= e
    
    return E
