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

A = np.array([0.0825,0.0825,0.0801,0.0890,0.0934,0.1621,0.2580])
B = np.array([0.0043,0.0043,0.0049,0.0053,0.0060,0.0173,0.0383])
W = np.array([1.9756,1.9879,1.9972,1.9990,1.9927])
M = np.array([0.0165,0.0343,0.0594,0.1478,0.3036])

C = np.array([0.2234,0.2244,0.2321])
D = np.array([0.6969,0.6974,0.7042])
u2=np.array([0.0001,0.01,0.05])

T=0.05
u = np.array([0.0001,0.005,0.015,0.035,0.05,0.15,0.25])
up=np.linspace(0.0001,0.25,100)

kT = np.array([0.005,0.01,0.02,0.05,0.1])
kTp = np.linspace(0.005,0.1,100)

alpha = T * np.log(2+2*np.cosh(up/T))
beta = (T**2) * ((np.pi**2)/6 + (up/T)**2 / 2)

omega = 2+0*kTp
emme = 3*kTp

plt.plot(up/T,alpha, c='b', label = r"simulated $\alpha(\mu,T)$")
plt.plot(up/T,beta, c='r', label = r"simulated $\beta(\mu,T)$")
plt.scatter(u/T,A, c='purple', label= r"theoretical $\alpha(\mu,T)$")
plt.scatter(u/T,B, c='orange', label= r"theoretical $\beta(\mu,T)$")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$\mu/k_B T$")
plt.legend()
plt.tight_layout()
plt.savefig("./img/QuadFermi.png")
plt.close()

plt.plot(kTp,omega, c='b', label=r"simulated $\omega_p$")
plt.plot(kTp,emme, c='r',  label=r"simulated $\eta$")
plt.scatter(kT,W, c='purple', label=r"simulated $\omega_p$")
plt.scatter(kT,M, c='orange',  label=r"simulated $\eta$")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$k_B T$")
plt.legend()
plt.tight_layout()
plt.savefig("./img/QuadMax.png")
plt.close()

plt.plot(up/T,alpha, c='b', label = r"simulated $\alpha(\mu,T)$")
plt.plot(up/T,1+0*beta, c='r', label = r"simulated $\beta(\mu,T)$")
plt.scatter(u2/T,C, c='purple', label= r"theoretical $\alpha(\mu,T)$")
plt.scatter(u2/T,D, c='orange', label= r"theoretical $\beta(\mu,T)$")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$\mu/k_B T$")
plt.legend()
plt.tight_layout()
plt.savefig("./img/LinFermi.png")
plt.close()