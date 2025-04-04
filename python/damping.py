# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
import scipy.optimize as opt
import mpmath as mpm

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 16})

0.2088,0.2088,0.2094,0.2066,0.2233,0.4444,0.7448

0.1764,0.1764,0.1826,0.2319,0.2394,0.3249,0.4084

vF=0.5

A = np.array([0.9888,0.9888,0.988,0.985,0.9823,0.959])
B = np.array([0.0115,0.0115,0.0115,0.01548,0.03345,0.1054])
W = np.array([1.9756,1.9879,1.9972,1.9990,1.9927])
M = np.array([0.0165,0.0343,0.0594,0.1478,0.3036])

C = np.array([0.963,0.912,0.81,0.574,0.43,0.234])
D = np.array([0.1837,0.1876,0.21,0.2,0.184,0.192])

T=0.01
T2=0.1
u = np.array([0.0001,0.001,0.01,0.03,0.1,0.3])
u2 = np.array([0.01,0.05,0.1,0.3,0.5,1])
up=np.linspace(0.00009,0.4,100)
up2=np.linspace(0.009,1.2,100)


kT = np.array([0.005,0.01,0.02,0.05,0.1])
kTp = np.linspace(0.005,0.1,100)

beta = np.zeros(100)
for i in range(100):
    beta[i]=3/2*T*(mpm.polylog(2,-np.exp(up[i]/T))/mpm.polylog(1,-np.exp(up[i]/T))+mpm.polylog(2,-np.exp(-up[i]/T))/mpm.polylog(1,-np.exp(-up[i]/T)))

alpha = 1+0*beta

alpha2 = np.copy(beta)
for i in range(100):
    alpha2[i]=mpm.polylog(1,-np.exp(up2[i]/T2))/mpm.polylog(2,-np.exp(up2[i]/T2))/(4*T2)*vF

beta2 = 0.75*vF**2+0*beta

omega = 2+0*kTp
emme = 3*kTp

plt.scatter(u/T,A, c='purple', label= r"theoretical $\alpha(\mu,T)$",zorder=10)
plt.scatter(u/T,3*B, c='orange', label= r"theoretical $\beta(\mu,T)$",zorder=10)
plt.plot(up/T,alpha, c='b', label = r"simulated $\alpha(\mu,T)$")
plt.plot(up/T,beta, c='r', label = r"simulated $\beta(\mu,T)$")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$\mu/k_B T$")
plt.legend()
plt.tight_layout()
plt.savefig("./img/QuadFermi.png")
plt.close()

plt.scatter(kT,W, c='purple', label=r"theoretical $\omega_p$",zorder=10)
plt.scatter(kT,M, c='orange',  label=r"theoretical $\eta$",zorder=10)
plt.plot(kTp,omega, c='b', label=r"simulated $\omega_p$")
plt.plot(kTp,emme, c='r',  label=r"simulated $\eta$")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$k_B T$")
plt.legend()
plt.tight_layout()
plt.savefig("./img/QuadMax.png")
plt.close()

plt.scatter(u2/T2,C, c='purple', label= r"theoretical $\alpha(\mu,T)$",zorder=10)
plt.scatter(u2/T2,D, c='orange', label= r"theoretical $\beta(\mu,T)$",zorder=10)
plt.plot(up2/T2,alpha2, c='b', label = r"simulated $\alpha(\mu,T)$")
plt.plot(up2/T2,beta2, c='r', label = r"simulated $\beta(\mu,T)$")
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r"$\mu/k_B T$")
plt.legend()
plt.tight_layout()
plt.savefig("./img/LinFermi.png")
plt.close()