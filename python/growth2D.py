# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
import scipy.optimize as opt

#############################
name = 'anomalous4'
iter = 2500
start = 30
end = 65
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
Nt = file.attrs['Number of Timesteps']
diag_f = file.attrs['Diagnostic Frequency']
pos_num = file.attrs['Number of Position Points']
dx = file.attrs['Position Delta']
dt = file.attrs['Time Delta']

kmax = 2*np.pi/dx
dk = 2*np.pi/(pos_max - pos_min)

n_figs = (Nt)//diag_f

Srho = np.zeros((n_figs, pos_num[0]-4, pos_num[1]-4))
Sphi = np.zeros((n_figs, pos_num[0]-4, pos_num[1]-4))
Mrho = np.zeros(n_figs)
Mphi = np.zeros(n_figs)
Kphi = np.zeros(n_figs)
Krho = np.zeros(n_figs)
Tphi = np.zeros(n_figs)
Trho = np.zeros(n_figs)

for i in range(n_figs):

    phi = file['/Fields'+str(i*diag_f)][0]['real']
    #charge = file['/Sources'+str(i*diag_f)]
    #rho = charge[0]['real']-charge[1]['real']

    phi = np.reshape(phi, (pos_num[0],pos_num[1]))
    #rho = np.reshape(rho, (pos_num[0],pos_num[1]))

    #rho = rho[2:-2,2:-2]
    phi = phi[2:-2,2:-2]

    #Frho = scp.fft.fftshift(scp.fft.fft(rho))
    Fphi = scp.fft.fftshift(scp.fft.fft2(phi))

    Sphi[i,:,:] = np.abs(Fphi)
    #Srho[i,:,:] = np.abs(Frho)

    Mphi[i] = np.max(Sphi[i,:,:])
    kx_phi = (np.argmax(np.abs(Fphi))//(pos_num[1]-4)) * dk[0] - kmax[0]/2
    ky_phi = (np.argmax(np.abs(Fphi))%(pos_num[1]-4)) * dk[1] - kmax[1]/2
    Kphi[i] = np.sqrt(kx_phi**2 + ky_phi**2)
    Tphi[i] = np.arctan(ky_phi/(kx_phi+0.0000001))*180/np.pi
    #Mrho[i] = np.max(np.abs(Frho))
    #Krho[i] = np.argmax(np.abs(Frho))
    #kx_rho = np.argmax(np.abs(Frho))[0] * dk[0] - kmax[0]/2
    #ky_rho = np.argmax(np.abs(Frho))[1] * dk[1] - kmax[1]/2
    #Krho[i] = np.sqrt(kx_rho**2 + ky_rho**2)
    #Trho[i] = np.arctan2(ky_rho,kx_rho)

    if(i*diag_f == iter):
        plt.imshow(np.abs(Fphi.T), extent = (-kmax[0]/2, kmax[0]/2, -kmax[1]/2, kmax[1]/2), aspect='auto', origin = 'lower')
        plt.colorbar()
        plt.savefig("./img/"+name+"_"+"field_fft.png", dpi=100)
        plt.close()
        print(phi)


"""
plt.plot(np.linspace(0,Nt*dt,n_figs), np.log(Srho))
plt.xlabel(r"$t$ [$\omega_{pe}^{-1}$]")
plt.ylabel(r"Log(Spectral Density)") 
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_growth"+".png", dpi=100)
plt.close()
"""
plt.plot(np.linspace(0,Nt*dt,n_figs-1), Mphi[1:])
plt.plot(np.linspace(start,end,n_figs-1), 2*np.exp(0.086*np.linspace(0,end-start,n_figs-1)))
plt.xlabel(r"$t$ [$\omega_{pe}^{-1}$]")
plt.ylabel(r"Log(Spectral Density)") 
plt.yscale('log')
plt.tight_layout()
plt.savefig("./img/"+name+"_field_growth"+".png", dpi=100)
plt.close()

plt.plot(np.linspace(0,Nt*dt,n_figs), Kphi)
plt.xlabel(r"$t$ [$\omega_{pe}^{-1}$]")
plt.ylabel(r"max($k$)") 
plt.tight_layout()
plt.savefig("./img/"+name+"_field_vector"+".png", dpi=100)
plt.close()

plt.plot(np.linspace(0,Nt*dt,n_figs), Tphi)
plt.xlabel(r"$t$ [$\omega_{pe}^{-1}$]")
plt.ylabel(r"$\theta_k$ [º]") 
plt.tight_layout()
plt.savefig("./img/"+name+"_field_angle"+".png", dpi=100)
plt.close()

def line_fit(x, a, b):
    return a*x+b

growth_rate = np.zeros((pos_num[0]-4,pos_num[1]-4))

start = int(start//(dt*diag_f))
end = int(end//(dt*diag_f))

for i in range(pos_num[0]-4):
    for j in range(pos_num[1]-4):
        spectre = np.log(Sphi[start:end,i,j])
        params, covs = opt.curve_fit(line_fit, np.linspace(start*dt*diag_f,end*dt*diag_f,end-start), spectre)
        growth_rate[i,j] = params[0]

plt.imshow(growth_rate.T, extent = (-kmax[0]/2, kmax[0]/2, -kmax[1]/2, kmax[1]/2), aspect='auto', origin = 'lower')
plt.colorbar()
plt.savefig("./img/"+name+"_"+"field_damping.png", dpi=100)
plt.close()

spectre = np.log(Mphi[start:end])
params, covs = opt.curve_fit(line_fit, np.linspace(start*dt*diag_f,end*dt*diag_f,end-start), spectre)
print(params[0])