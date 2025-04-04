# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

L=100
N=256
axis = np.linspace(0,L,N,endpoint=False)
axis2 = np.linspace(0,2*L,2*N,endpoint=False)
qs = scp.fft.fftfreq(N,d=L/N)*2*np.pi
qs2 = scp.fft.fftfreq(2*N,d=L/N)*2*np.pi


charge = np.exp(-(axis-10)**2/30**2)
kernel = np.concatenate(([0],1/(2*np.abs(qs[1:]))))
potential = scp.fft.ifft(kernel * scp.fft.fft(charge))
potential = potential-np.max(potential)+1

plt.plot(axis,charge, c='y')
plt.plot(axis,potential, c='r')
plt.savefig("./img/ffttest0.png")
plt.close()

charge2 = np.concatenate((charge,np.zeros(N)))
kernel2 = np.concatenate(([0],1/(2*np.abs(qs2[1:]))))
potential2 = scp.fft.ifft(kernel2 * scp.fft.fft(charge2))
potential2 = potential2-np.max(potential2)+1

plt.plot(axis,charge, c='y')
plt.plot(axis,potential, c='r')
plt.plot(axis,potential2[:N], c='r')
plt.savefig("./img/ffttest1.png")
plt.close()

charge3 = np.concatenate((charge2,np.zeros(2*N)))
kernel3 = scp.fft.fft(scp.fft.fftshift(np.concatenate((scp.fft.fftshift(scp.fft.ifft(kernel2)),np.zeros(2*N)))))
potential3 = scp.fft.fftshift(scp.fft.ifft(kernel3 * scp.fft.fft(charge3)))
potential3 = potential3-np.max(potential3)+1

plt.plot(axis,charge, c='y')
plt.plot(axis,potential3[N:2*N], c='r')
plt.savefig("./img/ffttest2.png")
plt.close()