import numpy as np

def convnfft(a, b):
    afft = np.fft.fftn(a, axes=(0,1,2))
    bfft = np.fft.fftn(b, axes=(0,1,2), s=a.shape)
    cfft = afft*bfft
    out = np.real(np.fft.ifftn(cfft, axes=(0,1,2)))
    for k in range(len(a.shape)):
        out = np.fft.ifftshift(out, k)

    return out