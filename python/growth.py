# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

#############################
name = 'anomalous4'
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
pos_num = file.attrs['Number of Position Points'][0]-4
mom_num = file.attrs['Number of Momentum Points'][0]
dx = file.attrs['Position Delta'][0]
dt = file.attrs['Time Delta']

kmax = 2*np.pi/dx
dk = 2*np.pi/(pos_max[0] - pos_min[0])

qs = np.linspace(-kmax/2,kmax/2,pos_num)

n_figs = (Nt)//diag_f

Srho = np.zeros(n_figs)
Sphi = np.zeros(n_figs)

for i in range(n_figs):

    phi = file['/Fields'+str(i*diag_f)][0]['real']
    charge = file['/Sources'+str(i*diag_f)]
    rho = charge[0]['real']-charge[1]['real']

    rho = rho[2:-2]
    phi = phi[2:-2]

    Frho = scp.fft.fftshift(scp.fft.fft(rho))
    Fphi = scp.fft.fftshift(scp.fft.fft(phi))

    Srho[i] = np.max(np.abs(Frho))
    Sphi[i] = np.max(np.abs(Fphi))

plt.plot(np.linspace(0,Nt*dt,n_figs), Srho)
plt.xlabel(r"$t$ [$\omega_{pe}^{-1}$]")
plt.ylabel(r"Spectral Density") 
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_growth"+".png", dpi=100)
plt.close()

plt.plot(np.linspace(0,Nt*dt,n_figs), Sphi)
plt.xlabel(r"$t$ [$\omega_{pe}^{-1}$]")
plt.ylabel(r"Spectral Density") 
plt.tight_layout()
plt.savefig("./img/"+name+"_field_growth"+".png", dpi=100)
plt.close()
