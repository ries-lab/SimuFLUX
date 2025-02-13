
"""
Copyright (c) 2022      Ries Lab, EMBL, Heidelberg, Germany
All rights reserved     Heintzmann Lab, Friedrich-Schiller-University Jena, Germany

@author: Rainer Heintzmann, Sheng Liu, Jonas Hellgoth
"""

from math import factorial 

import numpy as np
from scipy import signal
from scipy import fft

def fft3d_numpy(tfin):
    return np.fft.fftshift(np.fft.fftn(np.fft.fftshift(tfin, axes=(-1, -2, -3))), axes=(-1, -2, -3))

def ifft3d_numpy(tfin):
    return np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(tfin, axes=(-1, -2, -3))), axes=(-1, -2, -3))

def fft3d(tfin):
    return fft.fftshift(fft.fftn(fft.fftshift(tfin, axes=(-1, -2, -3))), axes=(-1, -2, -3))

def ifft3d(tfin):
    return fft.ifftshift(fft.ifftn(fft.ifftshift(tfin, axes=(-1, -2, -3))), axes=(-1, -2, -3))


def convolve(I_blur, kernel, bin):
    # replacement for tf.nn.convolution(I_blur,kernel,strides=(1,bin,bin,1),padding='SAME',data_format='NHWC')
    
    # Calculate padding
    pad_height = (kernel.shape[0] - 1) // 2
    pad_width = (kernel.shape[1] - 1) // 2

    # Apply padding
    I_blur_padded = np.pad(I_blur, ((0, 0), (pad_height, pad_height), (pad_width, pad_width), (0, 0)), mode='constant')

    # Perform convolution
    output = signal.convolve(I_blur_padded, kernel, mode='same', method='fft')

    # Apply strides
    output_strided = output[:, ::bin, ::bin, :]

    return output_strided

def prechirpz1(kpixelsize,pixelsize_x,pixelsize_y,N,M):
    krange = np.linspace(-N/2+0.5,N/2-0.5,N,dtype=np.float32)
    [xxK,yyK] = np.meshgrid(krange,krange)
    xrange = np.linspace(-M/2+0.5,M/2-0.5,M,dtype=np.float32)
    [xxR,yyR] = np.meshgrid(xrange,xrange)
    a = 1j*np.pi*kpixelsize
    A = np.exp(a*(pixelsize_x*xxK*xxK+pixelsize_y*yyK*yyK))
    A = A / A.shape[0]
    C = np.exp(a*(pixelsize_x*xxR*xxR+pixelsize_y*yyR*yyR))

    brange = np.linspace(-(N+M)/2+1,(N+M)/2-1,N+M-1,dtype=np.float32)
    [xxB,yyB] = np.meshgrid(brange,brange)
    B = np.exp(-a*(pixelsize_x*xxB*xxB+pixelsize_y*yyB*yyB))
    Bp = np.zeros((fft.next_fast_len(B.shape[0]), fft.next_fast_len(B.shape[1])), dtype=B.dtype)
    Bp[:B.shape[0], :B.shape[1]] = B
    Bh = fft.fft2(Bp)
    Bh = Bh / Bh.shape[0]

    return A,Bh,C

def cztfunc1(datain, param):
    A = param[0]
    Bh = param[1]
    C = param[2]
    N = A.shape[0]
    L = Bh.shape[0]
    M = C.shape[0]

    sh = datain.shape
    lsh = len(sh)
    ln = L-N

    if lsh == 2:
        Apad = np.pad(A*datain/N, ((0,ln),(0,ln)))
    elif lsh == 3:
        Apad = np.pad(A*datain/N, ((0,fft.next_fast_len(A.shape[0])-A.shape[0]),(0,ln),(0,ln)))
    elif lsh > 3:
        # Pad the input arrays
        Apad = np.concat((A*datain/N,np.zeros(datain.shape[:-1]+(ln,),dtype=np.complex64)),axis=-1)
        Apad = np.concat((Apad,np.zeros(Apad.shape[:-2]+(ln,Apad.shape[-1],),np.complex64)),axis=-2)

    # Perform 2D FFT
    Ah = fft.fft2(Apad)

    cztout = fft.ifft2(Ah * Bh[None,...])

    # Extract the relevant portion of the output
    dataout = C * cztout[..., -M:, -M:]

    return dataout

def noll2nl(j):
    n = np.ceil((-3 + np.sqrt(1 + 8*j)) / 2)
    l = j - n * (n + 1) / 2 - 1
    if np.mod(n, 2) != np.mod(l, 2):
       l = l + 1
    
    if np.mod(j, 2) == 1:
       l= -l
    
    return np.int32(n),np.int32(l)

def radialpoly(n,m,rho):
    if m==0:
        g = np.sqrt(n+1)
    else:
        g = np.sqrt(2*n+2)
    r = np.zeros(rho.shape)
    for k in range(0,(n-m)//2+1):
        coeff = g*((-1)**k)*factorial(n-k)/factorial(k)/factorial((n+m)//2-k)/factorial((n-m)//2-k)
        p = rho**(n-2*k)
        r += coeff*p

    return r

def genZern1(n_max,xsz):
    Nk = (n_max+1)*(n_max+2)//2
    Z = np.ones((Nk,xsz,xsz))
    pkx = 2/xsz
    xrange = np.linspace(-xsz/2+0.5,xsz/2-0.5,xsz)
    [xx,yy] = np.meshgrid(xrange,xrange)
    rho = np.lib.scimath.sqrt((xx*pkx)**2+(yy*pkx)**2)
    phi = np.arctan2(yy,xx)

    for j in range(Nk):
        n, l = noll2nl(j+1)
        m = np.abs(l)
        r = radialpoly(n,m,rho)
        if l<0:
            Z[j] = r*np.sin(phi*m)
        else:
            Z[j] = r*np.cos(phi*m)
    return Z
