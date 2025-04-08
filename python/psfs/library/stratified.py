#Transmission through stratified media. The wavefront
#is returned with respect to the media the objective
#was designed for. The media are inserted between the
#immersion and the sample. The design focus is at z=0
#(z pointing away from the objective).
#
# sys    System data
#  .eff     Effective situation
#    .nl       Refraction indices of the layers
#    .hl       Thickness of the layers [m]
#  .ref     Reference/design situation
#    .nl       Refraction indices of the layers
#    .hl       Thickness of the layers [m]
#  .lo      Wavelength in free space [m]
#  .nm      Refraction index of the immersion
#
# t      Total deflection angle [rad]
# tp     Transmission coefficients
# ts
# b2     Propagation kz/ko

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

from .interface import interface
from .layer import layer

def stratified(sys, t):
    go = sys['nm'] * np.sin(t)
    ko = 2j * np.pi / sys['lo']
    er = np.array([sys['nm']] + sys['eff']['nl']).astype(np.float64) ** 2
    Tp, Ts = transmission(np.array([sys['nm']] + sys['ref']['nl']).astype(np.float64) ** 2, sys['ref']['hl'], go, ko)
    tp, ts = transmission(er, sys['eff']['hl'], go, ko)
    tp = tp * np.conj(np.sign(Tp))
    ts = ts * np.conj(np.sign(Ts))
    Tp, Ts, b2 = interface(er[-1], sys['ns'] ** 2, go)
    tp = tp * Tp
    ts = ts * Ts

    return tp, ts, b2

def transmission(er, hz, go, ko):
    if len(er) < 2:
        tp = 1
        ts = 1
    else:
        n = np.imag(er) == 0  # avoid undefined results in resonances
        er[n] = er[n] + 1e-9j
        tp, ts, bz, rp, rs = interface(er[-2], er[-1], go)
        bz = np.exp(1j * hz[-1] * ko * bz)
        for n in range(len(er) - 2, 0, -1):
            ta, tb, pz, ra, rb = interface(er[n-1], er[n], go)
            pz = np.exp(1j * hz[n-1] * ko * pz)
            rp, tp = layer(ra, ta, rp, tp, pz)
            rs, ts = layer(rb, tb, rs, ts, pz)
        tp = tp * bz
        ts = ts * bz

    return tp, ts
