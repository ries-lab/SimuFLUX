import numpy as np

def conv2fft(a, b):
    afft = np.fft.fft2(a, axes=(0,1))
    bfft = np.fft.fft2(b, s=(a.shape[0], a.shape[1]))
    cfft = afft*bfft[...,None]
    c = np.real(np.fft.ifft2(cfft, axes=(0,1)))
    for k in [0,1]:
        c = np.fft.ifftshift(c, k)

    return c
