import numpy as np

def makeNPC(pos=[0,0,0], R=50, copynumber=32, dz=50, rotation=0, twistangle=np.pi/32, shiftangle=np.pi/16, dR=3):
    if not isinstance(pos, np.ndarray):
        pos = np.array(pos, dtype=float)

    dphi = np.pi/4
    maxphi = 2 * np.pi - dphi - rotation
    nphi = int(maxphi / dphi)
    phi = np.arange(nphi+1)*dphi
    v0 = np.zeros_like(phi)

    if copynumber == 8:
        phiall = phi
        zall = v0
        Rall = v0 + R
    elif copynumber == 16:
        phiall = np.hstack([phi, -phi+twistangle])
        zall = np.hstack([v0+dz/2, v0-dz/2])
        Rall = np.hstack([v0,v0]) + R
    elif copynumber == 32:
        phiall = np.hstack([phi, phi+shiftangle, -phi+twistangle, -phi + twistangle - shiftangle])
        zall = np.hstack([v0+dz/2, v0+dz/2, v0-dz/2, v0-dz/2])
        Rall = np.hstack([v0+R+dR, v0+R, v0+R+dR, v0+R])
    else:
        raise ValueError("NPC copy number: 8, 16 or 32")

    posnpc = np.vstack([Rall*np.cos(phiall), Rall*np.sin(phiall), zall]).T
    posnpc += pos

    return posnpc
