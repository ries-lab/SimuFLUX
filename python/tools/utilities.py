
"""
Copyright (c) 2022      Ries Lab, EMBL, Heidelberg, Germany
All rights reserved     Heintzmann Lab, Friedrich-Schiller-University Jena, Germany

@author: Rainer Heintzmann, Sheng Liu, Jonas Hellgoth
"""

import os
import sys
import ctypes
import pickle
import multiprocessing
from math import factorial

import numpy as np
from scipy import signal
from scipy import ndimage
from scipy.fft import next_fast_len

import pyfftw

WISDOMFILE = os.path.join(os.path.split(__file__)[0], 'fftw_wisdom.pkl')


NUM_THREADS = multiprocessing.cpu_count()
PLANNER_FLAGS = ('FFTW_MEASURE',)

pyfftw.config.NUM_THREADS = NUM_THREADS
pyfftw.config.PLANNER_EFFORT = PLANNER_FLAGS[0]

pyfftw.interfaces.cache.enable()

PYFFTW_ARRAYS = {}
PYFFTW_BUILDERS = {}

if sys.platform == 'win32':
    memcpy = ctypes.cdll.msvcrt.memcpy
elif sys.platform == 'darwin':
    memcpy = ctypes.CDLL('libSystem.dylib').memcpy
else: #linux
    memcpy = ctypes.CDLL('libc.so.6').memcpy

def load_wisdom(wisdomfile):
    if os.path.exists(wisdomfile):
        try:
            with open(wisdomfile, 'rb') as f:
                pyfftw.import_wisdom(pickle.load(f))
        except Exception:
            pass

load_wisdom(WISDOMFILE)

def save_wisdom(wisdomfile):
    with open(wisdomfile, 'wb') as f:
            pickle.dump(pyfftw.export_wisdom(), f)

def fft3d(tfin):
    a = pyfftw.empty_aligned(tfin.shape, dtype='complex64')
    b = pyfftw.builders.fftn(a, axes=(-3,-2,-1))
    a[:] = pyfftw.interfaces.numpy_fft.fftshift(tfin, axes=(-3, -2, -1))
    return pyfftw.interfaces.numpy_fft.fftshift(b(), axes=(-3, -2, -1))

def ifft3d(tfin):
    a = pyfftw.empty_aligned(tfin.shape, dtype='complex64')
    b = pyfftw.builders.ifftn(a, axes=(-3,-2,-1))
    a[:] = pyfftw.interfaces.numpy_fft.ifftshift(tfin, axes=(-3, -2, -1))
    return pyfftw.interfaces.numpy_fft.ifftshift(b(), axes=(-3, -2, -1))

# def convolve(I_blur, kernel, bin):
#     # replacement for tf.nn.convolution(I_blur,kernel,strides=(1,bin,bin,1),padding='SAME',data_format='NHWC')
    
#     # Calculate padding
#     pad_height = (kernel.shape[0] - 1) // 2
#     pad_width = (kernel.shape[1] - 1) // 2

#     # Apply padding
#     I_blur_padded = np.pad(I_blur, ((0, 0), (pad_height, pad_height), (pad_width, pad_width), (0, 0)), mode='constant')

#     # Perform convolution
#     output = signal.convolve(I_blur_padded, kernel, mode='same', method='fft')

#     # Apply strides
#     output_strided = output[:, ::bin, ::bin, :]

#     return output_strided

def convolve(I_blur, kernel, bin):
    # Output shape calculation
    output_shape = (I_blur.shape[0], I_blur.shape[1] // bin, I_blur.shape[2] // bin, I_blur.shape[3])
    output = np.zeros(output_shape, dtype=I_blur.dtype)
    
    # Perform convolution on each slice
    for i in range(I_blur.shape[0]):
        slice_I_blur = I_blur[i, :, :, :]
        convolved_slice = ndimage.convolve(slice_I_blur, kernel[:, :, :, 0], mode='constant', cval=0.0)
        output[i, :, :, :] = convolved_slice[::bin, ::bin]

    return output

# def prechirpz1(kpixelsize,pixelsize_x,pixelsize_y,N,M):
#     krange = np.linspace(-N/2+0.5,N/2-0.5,N,dtype=np.float32)
#     [xxK,yyK] = np.meshgrid(krange,krange)
#     xrange = np.linspace(-M/2+0.5,M/2-0.5,M,dtype=np.float32)
#     [xxR,yyR] = np.meshgrid(xrange,xrange)
#     a = 1j*np.pi*kpixelsize
#     A = np.exp(a*(pixelsize_x*xxK*xxK+pixelsize_y*yyK*yyK))
#     A = A / A.shape[0]
#     C = np.exp(a*(pixelsize_x*xxR*xxR+pixelsize_y*yyR*yyR))

#     brange = np.linspace(-(N+M)/2+1,(N+M)/2-1,N+M-1,dtype=np.float32)
#     [xxB,yyB] = np.meshgrid(brange,brange)
#     B = np.exp(-a*(pixelsize_x*xxB*xxB+pixelsize_y*yyB*yyB))
#     Bp = pyfftw.zeros_aligned((next_fast_len(B.shape[0]),next_fast_len(B.shape[1])), dtype=B.dtype)
#     b = pyfftw.builders.fftn(Bp, axes=(-1,-2))
#     Bp[:B.shape[0],:B.shape[1]] = B
#     Bh = b()
#     Bh = Bh / Bh.shape[0]
#     # Bh = np.fft.fft2(B)

#     return A,Bh,C

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
    Bp = pyfftw.zeros_aligned((next_fast_len(B.shape[0]),next_fast_len(B.shape[1])), dtype=B.dtype)
    b = pyfftw.builders.fftn(Bp, axes=(-1,-2))
    Bp[:B.shape[0],:B.shape[1]] = B
    Bh = b()
    Bh = Bh / Bh.shape[0]
    # Bh = np.fft.fft2(B)

    return A,Bh,C

# def cztfunc1(datain, param):
#     A = param[0]
#     Bh = param[1]
#     C = param[2]
#     N = A.shape[0]
#     L = Bh.shape[0]
#     M = C.shape[0]

#     sh = datain.shape
#     lsh = len(sh)
#     ln = L-N

#     if lsh == 2:
#         newsh = (sh[0]+ln, sh[1]+ln)
#     elif lsh == 3:
#         newsh = (next_fast_len(sh[0]), sh[1]+ln, sh[2]+ln)

#     key = str(newsh)
#     fkey, bkey = key+'forward', key+'backward'
#     try:
#         Apad, Ah = PYFFTW_ARRAYS['Apad'], PYFFTW_ARRAYS['Ah']
#         cztout = PYFFTW_ARRAYS['cztout']
#         b, d = PYFFTW_BUILDERS[fkey], PYFFTW_BUILDERS[bkey]
#     except KeyError:
#         PYFFTW_ARRAYS['Apad'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
#         PYFFTW_ARRAYS['Ah'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
#         PYFFTW_ARRAYS['cztout'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
#         PYFFTW_BUILDERS[fkey] = pyfftw.FFTW(PYFFTW_ARRAYS['Apad'], 
#                                             PYFFTW_ARRAYS['Ah'], 
#                                             axes=(-2,-1), 
#                                             direction='FFTW_FORWARD',
#                                             threads=NUM_THREADS,
#                                             flags=PLANNER_FLAGS)
#         PYFFTW_BUILDERS[bkey] = pyfftw.FFTW(PYFFTW_ARRAYS['Ah'], 
#                                             PYFFTW_ARRAYS['cztout'], 
#                                             axes=(-2,-1), 
#                                             direction='FFTW_BACKWARD',
#                                             threads=NUM_THREADS,
#                                             flags=PLANNER_FLAGS)
#         Apad, Ah = PYFFTW_ARRAYS['Apad'], PYFFTW_ARRAYS['Ah']
#         cztout = PYFFTW_ARRAYS['cztout']
#         b, d = PYFFTW_BUILDERS[fkey], PYFFTW_BUILDERS[bkey]

#     if lsh == 2:
#         Apad[:sh[0],:sh[1]] = A*datain
#     elif lsh == 3:
#         Apad[:sh[0],:sh[1],:sh[2]] = A*datain

#     # Perform 2D FFT
#     b.execute()

#     Ah *= Bh

#     # Perform 2D IFFT
#     d.execute()

#     # Extract the relevant portion of the output
#     dataout = C * cztout[..., -M:, -M:]

#     return dataout

# def cztfunc1(datain, param):
#     A = param[0]
#     Bh = param[1]
#     C = param[2]
#     N = A.shape[0]
#     L = Bh.shape[0]
#     M = C.shape[0]

#     sh = datain.shape
#     lsh = len(sh)
#     ln = L-N

#     if lsh == 2:
#         # newsh = (sh[0]+ln, sh[1]+ln)
#         newsh = (sh[0]+ln, sh[1]+ln)
#     elif lsh == 3:
#         # newsh = (next_fast_len(sh[0]), sh[1]+ln, sh[2]+ln)
#         newsh = (sh[0]+ln, sh[1]+ln, next_fast_len(sh[2]))

#     key = str(newsh)
#     fkey, bkey = key+'forward', key+'backward'
#     try:
#         Apad, Ah = PYFFTW_ARRAYS['Apad'], PYFFTW_ARRAYS['Ah']
#         cztout = PYFFTW_ARRAYS['cztout']
#         b, d = PYFFTW_BUILDERS[fkey], PYFFTW_BUILDERS[bkey]
#     except KeyError:
#         PYFFTW_ARRAYS['Apad'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
#         PYFFTW_ARRAYS['Ah'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
#         PYFFTW_ARRAYS['cztout'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
#         PYFFTW_BUILDERS[fkey] = pyfftw.FFTW(PYFFTW_ARRAYS['Apad'], 
#                                             PYFFTW_ARRAYS['Ah'], 
#                                             axes=(0,1), 
#                                             direction='FFTW_FORWARD',
#                                             threads=NUM_THREADS,
#                                             flags=PLANNER_FLAGS)
#         PYFFTW_BUILDERS[bkey] = pyfftw.FFTW(PYFFTW_ARRAYS['Ah'], 
#                                             PYFFTW_ARRAYS['cztout'], 
#                                             axes=(0,1), 
#                                             direction='FFTW_BACKWARD',
#                                             threads=NUM_THREADS,
#                                             flags=PLANNER_FLAGS)
#         Apad, Ah = PYFFTW_ARRAYS['Apad'], PYFFTW_ARRAYS['Ah']
#         cztout = PYFFTW_ARRAYS['cztout']
#         b, d = PYFFTW_BUILDERS[fkey], PYFFTW_BUILDERS[bkey]

#         save_wisdom(WISDOMFILE)

#     if lsh == 2:
#         Apad[:sh[0],:sh[1]] = A*datain

#         b.execute()

#         Ah *= Bh

#         d.execute()

#         cztout[-M:,-M:] *= C
#         return cztout[-M:,-M:].T

#     elif lsh == 3:
#         # Apad[:sh[0],:sh[1],:sh[2]] = datain * A[...,None]

#         # Apad[:sh[0],:sh[1],:sh[2]] = datain
#         # Apad[:sh[0],:sh[1],:sh[2]] *= A[...,None]

#         # memcpy(Apad[:sh[0],:sh[1],:sh[2]].ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)), 
#         #        datain.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
#         #        datain.nbytes)

#         # print(f"A.dtype: {A.dtype} datain.dtype: {datain.dtype} Apad.dtype: {Apad.dtype}")

#         Apadp = datain*A[...,None]

#         Apad[:sh[0],:sh[1],:sh[2]] = Apadp[:]

#         # np.multiply(A[...,None], datain, out=Apad[:sh[0],:sh[1],:sh[2]])

#         b.execute()

#         Ah *= Bh[...,None]

#         d.execute()

#         cztout[-M:,-M:,:] *= C[...,None]
#         return cztout[-M:,-M:,:].T

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
        newsh = (sh[0]+ln, sh[1]+ln)
    elif lsh == 3:
        # newsh = (next_fast_len(sh[0]), sh[1]+ln, sh[2]+ln)
        newsh = (sh[0], sh[1]+ln, sh[2]+ln)

    # print("sh ", sh)
    # print("newsh ", newsh)

    key = str(newsh)
    fkey, bkey = key+'forward', key+'backward'
    try:
        Apad, Ah = PYFFTW_ARRAYS[key+'Apad'], PYFFTW_ARRAYS[key+'Ah']
        cztout = PYFFTW_ARRAYS[key+'cztout']
        b, d = PYFFTW_BUILDERS[fkey], PYFFTW_BUILDERS[bkey]
    except KeyError:
        PYFFTW_ARRAYS[key+'Apad'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
        PYFFTW_ARRAYS[key+'Ah'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
        PYFFTW_ARRAYS[key+'cztout'] = pyfftw.zeros_aligned(newsh, dtype='complex64')
        PYFFTW_BUILDERS[fkey] = pyfftw.FFTW(PYFFTW_ARRAYS[key+'Apad'], 
                                            PYFFTW_ARRAYS[key+'Ah'], 
                                            axes=(-2,-1), 
                                            direction='FFTW_FORWARD',
                                            threads=NUM_THREADS,
                                            flags=PLANNER_FLAGS)
        PYFFTW_BUILDERS[bkey] = pyfftw.FFTW(PYFFTW_ARRAYS[key+'Ah'], 
                                            PYFFTW_ARRAYS[key+'cztout'], 
                                            axes=(-2,-1), 
                                            direction='FFTW_BACKWARD',
                                            threads=NUM_THREADS,
                                            flags=PLANNER_FLAGS)
        
        save_wisdom(WISDOMFILE)
        Apad, Ah = PYFFTW_ARRAYS[key+'Apad'], PYFFTW_ARRAYS[key+'Ah']
        cztout = PYFFTW_ARRAYS[key+'cztout']
        b, d = PYFFTW_BUILDERS[fkey], PYFFTW_BUILDERS[bkey]

    if lsh == 2:
        Apad[:sh[0],:sh[1]] = A*datain

        b.execute()

        Ah *= Bh

        d.execute()

        cztout[-M:,-M:] *= C
        return cztout[-M:,-M:]

    elif lsh == 3:
        # Apad[:sh[0],:sh[1],:sh[2]] = datain * A[...,None]

        # Apad[:sh[0],:sh[1],:sh[2]] = datain
        # Apad[:sh[0],:sh[1],:sh[2]] *= A[...,None]

        # memcpy(Apad[:sh[0],:sh[1],:sh[2]].ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)), 
        #        datain.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8)),
        #        datain.nbytes)


        for i in range(sh[0]):
            Apad[i, :sh[1], :sh[2]] = A*datain[i,...]

        # np.multiply(A[...,None], datain, out=Apad[:sh[0],:sh[1],:sh[2]])

        b.execute()

        Ah *= Bh

        d.execute()

        cztout[...,-M:,-M:] *= C
        return cztout[...,-M:,-M:]

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
