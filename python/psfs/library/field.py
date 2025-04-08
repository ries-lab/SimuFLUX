#Electric field transmitted through an objective.
#
# sys    System data
# out    Simulation results
#
# E      Electric field [V/m]
# k      Wavevector [rad/m]
# N      Number of samples

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

import matplotlib.pyplot as plt
import numpy as np

from .samples import samples
from .modifiers import modifiers
from .stratified import stratified
from .transmission import transmission
from .interface import interface
from .cis import cis

strat = getattr(globals(), 'strat', None)

def field(sys, out):
    M, N = samples(sys, out)  # Sampling over aperture and lateral wavevector.
    E, r, t, p, ct, st = modifiers(sys, M, N)  # Execute field modifiers.

    # print(f"post modifiers: {E.shape}")

    if sys.get('phaseImage', False):
        plt.figure()
        plt.imshow(np.angle(E[0, :, :]).squeeze(), [np.min(np.angle(E[0, :, :])), np.max(np.angle(E[0, :, :]))])
        plt.figure()
        plt.imshow(np.angle(E[1, :, :]).squeeze(), [np.min(np.angle(E[1, :, :])), np.max(np.angle(E[1, :, :]))])

    if sys.get('intImage', False):
        plt.figure()
        plt.imshow(np.abs(E[0, :, :]).squeeze(), [np.min(np.abs(E[0, :, :])), np.max(np.abs(E[0, :, :]))])
        plt.figure()
        plt.plot(np.abs(E[0, :, 29]))

    if strat is not None:
        cp, sp, k = stratified(sys, t)  # Transmission through stratified media.
    else:
        cp, sp, k, _, _ = interface(sys['nm']**2, sys['ns']**2, sys['NA'] * r)  # Fresnel reflection and transmission coefficients at an interface.
    
    tp, ts = transmission(sys, t)  # Transmission through an air-lens-air-immersion objective.
    tp = tp*cp
    ts = ts*sp
    cp, sp = cis(p, ret_both=True)  # Sine and cosine.

    try:
        tp = tp * (cp * E[0, :, :] + sp * E[1, :, :])
        ts = ts * (sp * E[0, :, :] - cp * E[1, :, :])
    except (IndexError, TypeError):
        tp = tp * cp * E
        ts = ts * sp * E

    ct = tp * k / sys['ns']
    E = np.stack([ct * cp + sp * ts, ct * sp - cp * ts, tp * st / sys['ns']])
    k = 2 * np.pi / sys['lo'] * k
    E[np.isnan(E)] = 0

    return E, k, N