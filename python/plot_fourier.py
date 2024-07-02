# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
import scipy.optimize as opt

#############################
name = 'dispersion_2DFermiQuadratic6'
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

# Creates movie
fig = plt.figure()

# Reads Test_Data.h5 file in ../data/ directory
file=h5py.File('./output/'+name+'.h5','r')

pos_min = file.attrs['Position Min.']
pos_max = file.attrs['Position Max.']
mom_min = file.attrs['Momentum Min.']
mom_max = file.attrs['Momentum Max.']
Nt = file.attrs['Number of Timesteps']
diag_f = file.attrs['Diagnostic Frequency']
pos_num = file.attrs['Number of Position Points'][0]
mom_num = file.attrs['Number of Momentum Points'][0]
dx = file.attrs['Position Delta'][0]
dt = file.attrs['Time Delta']
n_figs = Nt//diag_f -1

Dt = dt * diag_f
wmax = 2*np.pi/Dt
kmax = 2*np.pi/dx
dw = 2*np.pi/(Nt*dt)
dk = 2*np.pi/(pos_max[0] - pos_min[0])

padding=2

phi = np.zeros((n_figs,pos_num-2*padding))
rho = np.zeros((n_figs,pos_num-2*padding))

print("dx = ",dx)
print("dt = ",dt)
print("Dt = ",Dt)
print("wmax = ",wmax)
print("kmax = ",kmax)
print("dw = ", dw)
print("dk = ", dk)

for i in range(n_figs):

    potential = file['/Fields'+str((i+1)*diag_f)][0]['real']
    charge1 = file['/Sources'+str((i+1)*diag_f)][0]['real']
    charge2 = file['/Sources'+str((i+1)*diag_f)][1]['real']

    phi[i,:] = potential[padding:-padding]
    rho[i,:] = charge1[padding:-padding]-charge2[padding:-padding]
    #rho[i,:] = charge1[padding:-padding]

    #print(np.average(phi[i,:]))

Frho = scp.fft.fftshift(scp.fft.fft2(rho))
Fphi = scp.fft.fftshift(scp.fft.fft2(phi))

Srho = np.abs(Frho)**2
Sphi = np.abs(Fphi)**2

wpe = 1
ve = 0.2
T = 0.05
u = 0
ge = np.pi
alpha = T * np.log(2+2*np.cosh(u/T))
beta = (T**2) * ((np.pi**2)/6 + (u/T)**2 / 2) 

qs = np.linspace(0,50,1000)

sc = plt.imshow(np.log10(Srho+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
#plt.plot(qs, np.sqrt(wpe**2 + 3/2 * ve**2 * qs**2), c = 'r') # 3D quadratic
plt.plot(qs, np.sqrt(ge * alpha * qs + 3 * beta / alpha * qs**2), c = 'r')     # 2D quadratic
#plt.plot(qs, np.sqrt(0.2/(4*np.pi) * qs + qs**2), c = 'r')    # 2D linear
plt.xlim((-10,10))
plt.ylim((0,6))
#plt.clim((-10,7))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$\omega$[$\omega_{pe}$]")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_chargeFFT.png", dpi=200)
plt.close()

sc = plt.imshow(np.log10(Sphi+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
#plt.plot(qs, np.sqrt(wpe**2 + 3/2 * ve**2 * qs**2), c = 'r') # 3D quadratic
plt.plot(qs, np.sqrt(ge * alpha * qs + 3 * beta / alpha * qs**2), c = 'r')     # 2D quadratic
#plt.plot(qs, np.sqrt(qs * u + qs**2), c = 'r')                # 2D linear
plt.xlim((-10,10))
plt.ylim((0,6))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$\omega$[$\omega_{pe}$]")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_fieldFFT.png", dpi=200)
plt.close()

kays = np.linspace(0,kmax/2,(pos_num-4)//2)

omega = np.argmax(Srho[n_figs//2: , (pos_num-4)//2:], axis = 0) * dw

def disp_rel(x, a, b):
    return a*x+b*x*x

limit = 20

params, cov = opt.curve_fit(disp_rel,kays[:limit],omega[:limit]**2)
print(params[0])
print(params[1])

plt.plot(kays, omega)
plt.plot(kays[:limit], np.sqrt(params[0]*kays[:limit]+params[1]*kays[:limit]**2))
plt.xlim((0,5))
plt.tight_layout()
plt.savefig("./img/"+name+"_omega.png", dpi=200)
plt.close()

Drho = np.abs(scp.fft.fftshift(scp.fft.fft(rho,axis=-1),axes=-1))
Dphi = np.abs(scp.fft.fftshift(scp.fft.fft(phi,axis=-1),axes=-1))

sc = plt.imshow(np.log(Drho+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,0,Nt*dt), aspect='auto', origin = 'lower')
plt.xlim((-10,10))
#plt.ylim((0,6))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$t$[$\omega_{pe}^{-1}$]")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_chargeDamp.png", dpi=200)
plt.close()

sc = plt.imshow(np.log(Dphi+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,0,Nt*dt), aspect='auto', origin = 'lower')
plt.xlim((-10,10))
#plt.ylim((0,6))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$t$[$\omega_{pe}^{-1}$]")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_fieldDamp.png", dpi=200)
plt.close()
