#[out,err]=effField(sys,out,opt)
#-------------------------------
#
#Electric field near the focus of an objective with the
#chirp z-transform implementation of the Debye integral.
#
#Input:
# sys    System data
# out    Result data
# opt    Options
#  .Et      Display transmitted field           {0}
#  .Ef      Display focus field                 {1}
#  .mem     System memory (RAM) [MB]        {512MB}
#
#Output:
# out    Result data
#  .x
#  .y       Coordinates [m]
#  .z
#  .E       Electric field [V/m]
# err    Error message

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

from .field import field
from .cis import cis
from .optdim import optdim

def effField(sys, out=None, opt=None):
    if opt is None or not isinstance(opt, dict):
        opt = {}

    if out is None:
        out = {}
        out['dr'] = opt.get('pixSize', 1)
        out['dz'] = opt.get('pixSize', 1)

    m = int(np.ceil(sys['rz'] / out['dz']))
    if 'zo' in sys:
        z = sys['zo'] + out['dz'] * np.arange(-m, m+1)
    else:
        z = out['dz'] * np.arange(m+1)
        sys['zo'] = 0

    if 'ns' not in sys:
        sys['ns'] = sys['nm']
    elif sys['ns'] < sys['nm']:
        z = z[z >= 0]
        if len(z) == 0:
            raise ValueError('Nothing to do: region is only in immersion.')

    nz = len(z)
    z = z.reshape(1, 1, nz)
    Et, kz, m = field(sys, out)
    # print("Et: ", Et[0,98:104,0])
    dk = 2 * np.pi / m
    out['Et'] = Et

    if z.ravel()[0]:
        Et = Et * np.exp(1j * z.ravel()[0] * kz)[None,...]

    m = int(np.ceil(sys['rx'] / out['dr']))
    if 'xo' in sys:
        out['x'] = sys['xo'] + out['dr'] * np.arange(-m, m+1).reshape(-1, 1)
    else:
        out['x'] = out['dr'] * np.arange(m+1).reshape(-1, 1)

    m = int(np.ceil(sys['ry'] / out['dr']))
    if 'yo' in sys:
        out['y'] = sys['yo'] + out['dr'] * np.arange(-m, m+1)
    else:
        out['y'] = out['dr'] * np.arange(m+1)

    kz = np.exp(1j * out['dz'] * kz)
    nx = len(out['x'])
    ny = len(out['y'])
    out['z'] = z

    if Et.dtype == np.float32:
        s = 24 / 1024**2
    else:
        s = 48 / 1024**2

    if s * nx * ny * nz > max(opt.get('mem', 4096) - 128, opt.get('mem', 4096) / 2):
        raise MemoryError('Out of memory.')

    s, m, n = Et.shape
    M = optdim(m + nx - 1)
    N = optdim(n + ny - 1)
    
    cout = cis(-dk / 2 * np.arange(m)**2 - dk * out['x'][0] / out['dr'] * np.arange(m))
    Et = Et * cout[None, :, None]
    cisa0 = cis(dk / 2 * ((m-1) * np.arange(nx) - np.arange(nx)**2))
    cisa1 = cis(-dk / 2 * np.arange(n)**2 - dk * out['y'][0] / out['dr'] * np.arange(n))
    # print("cisa0", cisa0[:3])
    # print("cisa1", cisa1[:3])
    a = cisa0[:,None]*cisa1[None,:]
    b = cis(dk / 2 * ((n-1) * np.arange(ny) - np.arange(ny)**2))
    v = np.fft.fft(cis(dk / 2 * (np.arange(1 - m, nx)**2)), M)
    w = np.fft.fft(cis(dk / 2 * (np.arange(1 - n, ny)**2)), N)
    # print(Et.shape, a.shape, b.shape, v.shape, w.shape) # (3, 207, 207) (141, 207) (141,) (384,) (384,)
    # print(f"Et: {Et[0,0,:3]}, a: {a[0,:3]}, a T: {a[:3,0]}, b: {b[:3]}, v: {v[:3]}, w: {w[:3]}")

    # print(f"m: {m}, n: {n}, nx: {nx}, ny: {ny}")
    m = np.arange(m-1, m + nx -1)
    n = np.arange(n-1, n + ny -1)

    out['E'] = np.zeros((3, nx, ny, nz), dtype=Et.dtype)

    for z in range(nz):
        E = np.fft.ifft(np.fft.fft(Et, M, axis=1) * v[None, :, None], M, axis=1)
        # print("first, ", E[0,0,:3])
        # print(f"shape: {E[:, m, :].shape} {a.shape}")
        # print("intermediate", (E[:, m, :] * a[None,...])[0,0,:3])
        # print("intermediate T", (E[:, m, :] * a[None,...])[0,:3,0])
        # print("intermediate2", (np.fft.fft(E[:, m, :] * a[None,...], N, axis=2) * w[None, None, :])[0,0,:3])
        E = np.fft.ifft(np.fft.fft(E[:, m, :] * a[None,...], N, axis=2) * w[None, None, :], N, axis=2)
        # print("second, ", E[0,0,:3])
        E = E[:, :, n] * b[None, None, :]
        # print("third, ", E[0,0,:3])
        out['E'][:, :, :, z] = E

        # raise UserWarning("test")

        Et = Et * kz[None,...]

    return out
