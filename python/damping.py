# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
import scipy.optimize as opt
import mpmath as mpm

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

0.2088,0.2088,0.2094,0.2066,0.2233,0.4444,0.7448

0.1764,0.1764,0.1826,0.2319,0.2394,0.3249,0.4084

A = np.array([0.2088,0.2088,0.2094,0.2230,0.2380,0.4444,1.5184])
B = np.array([0.1764,0.1764,0.1826,0.2006,0.2115,0.3249,0.4906])

T=0.05
u = np.array([0.0001,0.005,0.015,0.035,0.05,0.15,0.25])

up=np.linspace(0.0001,0.25,100)

alpha = T * np.log(2+2*np.cosh(up/T))
beta = (T**2) * ((np.pi**2)/6 + (up/T)**2 / 2)

plt.plot(up,np.pi*alpha, c='b')
plt.plot(up,3*beta/alpha, c='r')
plt.scatter(u,A, c='purple')
plt.scatter(u,B, c='orange')
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.savefig("./img/QuadFermi.png")
plt.close()

print(np.pi*T * np.log(2+2*np.cosh(0.25/T)))