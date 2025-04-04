# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

factor=1
L=100
N=256*factor
axis = np.linspace(-factor/2*L,factor/2*L,N,endpoint=False)
charge = np.exp(-(axis+10)**2)
qs = scp.fft.fftfreq(N,d=factor*L/N)*2*np.pi
qs2 = scp.fft.fftfreq(2*N,d=factor*L/N)*2*np.pi
kernel = np.concatenate(([0],1/(2*np.abs(qs[1:]))))


charge = np.zeros(N)
charge[int(N/4)] = 1
stencil = np.concatenate((np.linspace(0,1,int(N/2),endpoint=False),np.linspace(1,0,int(N/2),endpoint=False)))

axis = np.linspace(0,L,N,endpoint = False)
axis2 = np.linspace(0,2*L,2*N,endpoint = False)

res = scp.signal.fftconvolve(charge, stencil, mode = "same")
plt.plot(axis, res)
plt.plot(axis, charge)
plt.savefig("./img/ffttest_a.png")
plt.close()

kernel = scp.fft.fft(scp.fft.fftshift(stencil))
res2 = scp.fft.ifft(scp.fft.fft(charge)*kernel)
plt.plot(axis, res2)
plt.plot(axis, charge)
plt.savefig("./img/ffttest_b.png")
plt.close()

kernel2 = scp.fft.fft(scp.fft.fftshift(np.concatenate((scp.fft.fftshift(scp.fft.ifft(kernel)),np.zeros(N)))))
charge2 = np.concatenate((charge,np.zeros(N)))
res4 = scp.fft.fftshift(scp.fft.ifft(scp.fft.fft(charge2)*kernel2))
plt.plot(axis, res4[int(N/2):3*int(N/2)])
plt.plot(axis, charge)
plt.savefig("./img/ffttest_c.png")
plt.close()

plt.plot(axis,scp.fft.ifft(kernel))
plt.plot(axis2,scp.fft.ifft(kernel2))
#plt.xlim((-1,1))
plt.savefig("./img/ffttest_d.png")
plt.close()