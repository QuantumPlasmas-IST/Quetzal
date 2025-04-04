# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

#############################
name = 'weibel4'
Ktest = 3
fitC = -2
fitM = 0.2
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})

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

print("tmax =",Nt*dt)

kmax = 2*np.pi/dx
dk = 2*np.pi/(pos_max[0] - pos_min[0])

nk = int(Ktest/dk)

Dt = dt*diag_f
wmax = 2*np.pi/Dt
dw = 2*np.pi/(Nt*dt)

qs = np.linspace(-kmax/2,kmax/2,pos_num)

n_figs = (Nt)//diag_f

Srho = np.zeros(n_figs)
Sphi = np.zeros(n_figs)

RFrho = np.zeros((pos_num,n_figs))
RFphi = np.zeros((pos_num,n_figs))

fullRho = np.zeros((n_figs,pos_num))
fullPhi = np.zeros((n_figs,pos_num))

print(nk*dk)

for i in range(n_figs):

    phi = file['/Fields'+str(i*diag_f)][2]['real']
    charge = file['/Sources'+str(i*diag_f)]
    rho = charge[0]['real']#-charge[1]['real']

    Frho = scp.fft.fftshift(scp.fft.fft(rho))
    Fphi = scp.fft.fftshift(scp.fft.fft(phi))

    Srho[i] = np.abs(Frho)[pos_num//2+nk]
    Sphi[i] = np.abs(Fphi)[pos_num//2+nk]

    RFrho[:,i] = np.abs(Frho)[:]
    RFphi[:,i] = np.abs(Fphi)[:]

    fullRho[i,:] = rho
    fullPhi[i,:] = phi 

Orho = np.abs(scp.fft.fftshift(scp.fft.fft2(fullRho)))**2
Ophi = np.abs(scp.fft.fftshift(scp.fft.fft2(fullPhi)))**2

timebase = np.linspace(dt*diag_f, dt*diag_f*n_figs, n_figs-1)

print(np.sqrt(-scp.special.jv(1,nk*dk)*nk*dk))

plt.plot(timebase, Srho[1:])
plt.plot(timebase, np.min(Srho)*np.exp(fitM*(timebase-fitC)))
plt.xlabel(r"$t$ [$\omega_{pe}^{-1}$]")
plt.ylabel(r"Spectral Density") 
plt.yscale('log')
#plt.xlim((0,50))
plt.ylim((0.5*np.min(Srho), 2*np.max(Srho)))
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_growth"+".png", dpi=100)
plt.close()

def line(x,a,b):
    return a*x+b

timebase2=np.linspace(0, dt*diag_f*n_figs, n_figs)
totK = int(5//dk)
startK = 0
fit_Ms = np.zeros(totK)
fit_Bs = np.zeros(totK)
start_time=150
end_time=250
time_cond = np.logical_and(timebase2>start_time,timebase2<end_time)
for k in range(totK):
    params, covs = scp.optimize.curve_fit(line, timebase2[time_cond], np.log(np.abs(RFphi[pos_num//2+startK+k,time_cond])))
    fit_Ms[k] = params[0]
    fit_Bs[k] = params[1]
    #if(np.min(np.abs(RFphi[pos_num//2+startK+k,1:]))<np.abs(RFphi[pos_num//2+startK+k,1])):
    #    fit_Ms[k] = 0
    #if(covs[0,0]>0.005):
    #    fit_Ms[k] = 0

plt.plot(timebase2, np.log(np.abs(RFphi[pos_num//2+nk,:])))
plt.plot(timebase2[time_cond], line(timebase2[time_cond],fit_Ms[nk-startK],fit_Bs[nk-startK]))
plt.xlim((start_time,end_time))
#plt.ylim((0,2))
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_test.png", dpi=200)
plt.close()

plt.plot(np.linspace(dt*diag_f, dt*diag_f*n_figs, n_figs-1), Sphi[1:])
plt.plot(timebase, np.min(Sphi)*np.exp(fitM*(timebase-fitC)))
plt.xlabel(r"$t$ [$\omega_{pe}^{-1}$]")
plt.ylabel(r"Spectral Density") 
plt.yscale('log')
#plt.xlim((0,500))
plt.ylim((0.5*np.min(Sphi), 2*np.max(Sphi)))
plt.tight_layout()
plt.savefig("./img/"+name+"_field_growth"+".png", dpi=100)
plt.close()

T = 10
alpha = T * np.log(4)
beta = (T**2) * ((np.pi**2)/6) 
sc = plt.imshow(np.log10(Orho+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
#plt.plot(qs, np.sqrt(np.pi * 2 * 10 * qs + 3*beta/alpha * qs**2), c = 'r')     # 2D quadratic
plt.xlim((-0.5,0.5))
plt.ylim((0,6))
#plt.clim((-10,7))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$\omega$[$\omega_{pe}$]")
plt.colorbar(sc)
plt.tight_layout()
#plt.savefig("./img/"+name+"_charge_omegas.png", dpi=200)
plt.close()

sc = plt.imshow(np.log10(Ophi+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
#plt.xlim((-10,10))
#plt.ylim((0,6))
#plt.clim((-10,7))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$\omega$[$\omega_{pe}$]")
plt.colorbar(sc)
plt.tight_layout()
#plt.savefig("./img/"+name+"_field_omegas.png", dpi=200)
plt.close()

qs = np.linspace(0,50,1000)
plt.scatter(range(totK)*dk+startK*dk, fit_Ms)
plt.plot(qs, np.sqrt(-scp.special.jv(1,qs)*qs), c = 'r')
plt.xlim((0,12.5))
plt.ylim((0,0.1))
plt.xlabel(r"$k a$")
plt.ylabel(r"growth-rate $\gamma/\omega_{coll}$")
plt.tight_layout()
plt.savefig("./img/"+name+"_growth_rates.png", dpi=200)
plt.close()
