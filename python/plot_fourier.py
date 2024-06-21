# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

#############################
name = 'twostream2D'
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
n_figs = Nt//diag_f

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

    potential = file['/Fields'+str(i*diag_f)][0]['real']
    charge = file['/Sources'+str(i*diag_f)][0]['real']

    phi[i,:] = potential[padding:-padding]
    rho[i,:] = charge[padding:-padding]

Frho = scp.fft.fftshift(scp.fft.fft2(rho))
Fphi = scp.fft.fftshift(scp.fft.fft2(phi))

Srho = np.abs(Frho)**2
Sphi = np.abs(Fphi)**2

wpe = 1
ve = 0.2
T = 0.001
u = 0.2

qs = np.linspace(0,50,1000)

sc = plt.imshow(np.log(Srho+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
#plt.plot(qs, np.sqrt(wpe**2 + 3/2 * ve**2 * qs**2), c = 'r') # 3D quadratic
plt.plot(qs, np.sqrt(qs + 3/2 * ve**2 * qs**2), c = 'r')     # 2D quadratic
#plt.plot(qs, np.sqrt(0.2/(4*np.pi) * qs + qs**2), c = 'r')    # 2D linear
plt.xlim((-10,10))
plt.ylim((0,6))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$\omega$[$\omega_{pe}$]")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_chargeFFT.png", dpi=200)
plt.close()

sc = plt.imshow(np.log(Sphi+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
#plt.plot(qs, np.sqrt(wpe**2 + 3/2 * ve**2 * qs**2), c = 'r') # 3D quadratic
plt.plot(qs, np.sqrt(qs + 3/2 * ve**2 * qs**2), c = 'r')     # 2D quadratic
#plt.plot(qs, np.sqrt(qs * u + qs**2), c = 'r')                # 2D linear
plt.xlim((-10,10))
plt.ylim((0,6))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$\omega$[$\omega_{pe}$]")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_fieldFFT.png", dpi=200)
plt.close()


print(pos_num)
print(kmax/dk)
    
