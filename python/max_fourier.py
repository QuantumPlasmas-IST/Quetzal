# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
import scipy.optimize as opt
import mpmath as mpm

#############################
name = 'dispersion_2DMax'
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

# Creates movie
fig = plt.figure()

def line(x,a,b):
    return a*x+b

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

T = 0.01
u = 0.15
wpe = np.sqrt(2)
ge = np.pi
alpha = T * np.log(2+2*np.cosh(u/T))
beta = (T**2) * ((np.pi**2)/6 + (u/T)**2 / 2) 

qs = np.linspace(0,50,1000)

sc = plt.imshow(np.log10(Srho+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
plt.plot(qs, np.sqrt(wpe**2 + 3 * T * qs**2), c = 'r') # 3D quadratic
#plt.plot(qs, np.sqrt(ge * alpha * qs + 3*beta/alpha * qs**2), c = 'r')     # 2D quadratic
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
plt.plot(qs, np.sqrt(wpe**2 + 3 * T * qs**2), c = 'r') # 3D quadratic
#plt.plot(qs, np.sqrt(ge * alpha * qs + 3*beta/alpha * qs**2), c = 'r')     # 2D quadratic
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
    #return a*x+b*x*x
    return a + b*x*x

limit=20

if(name == 'dispersion_2DFermiQuadratic0'):
    limit = 9
if(name == 'dispersion_2DFermiQuadratic1'):
    limit = 9
if(name == 'dispersion_2DFermiQuadratic2'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic3'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic4'):
    limit = 10

if(name == 'dispersion_3DMaxwellQuadratic0'):
    limit = 9
if(name == 'dispersion_3DMaxwellQuadratic1'):
    limit = 9
if(name == 'dispersion_3DMaxwellQuadratic2'):
    limit = 10
if(name == 'dispersion_3DMaxwellQuadratic3'):
    limit = 10
if(name == 'dispersion_3DMaxwellQuadratic4'):
    limit = 10

params, cov = opt.curve_fit(disp_rel,kays[1:limit],omega[1:limit]**2)

plt.plot(kays, omega)
#plt.plot(kays[:limit], np.sqrt(params[0]*kays[:limit]+params[1]*kays[:limit]**2))
plt.plot(kays[:limit], np.sqrt(params[0]+params[1]*kays[:limit]**2))
plt.xlim((0,5))
plt.tight_layout()
plt.savefig("./img/"+name+"_omega.png", dpi=200)
plt.close()

print(ge*alpha,params[0])
print(3*beta/alpha,params[1])

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

times = np.linspace(0,Nt*dt,n_figs)

def maxwell_derivative(func_p, func_T):
    return -2/np.sqrt(2*np.pi*func_T) * func_p/func_T * np.exp(-func_p*func_p/(2*func_T))


def damping_rate(func_k, func_T):
    func_omegaP = 2
    func_beta = 3*T
    func_w = np.sqrt(func_omegaP + func_beta * func_k**2)
    func_p = func_w/func_k
    return np.pi * func_w / (2 * func_k**2) * maxwell_derivative(func_p,func_T)

exp_damps = np.zeros(50)

ptest=35

for it in range(50):

    slic = it+1
    k_slic = dk*slic

    signal_raw = np.log(Drho[:,(pos_num-4)//2+slic])
    signal_raw = signal_raw-np.average(signal_raw)

    b, a = scp.signal.butter(5, 0.003*slic+0.01, btype='low')
    signal_raw[0] = signal_raw[1]
    signal_filtered = scp.signal.medfilt(scp.signal.lfilter(b, a, signal_raw), 1)

    pro = 0.001
    n_peak=0

    divs = 2

    if(slic>=20):
        divs=2
    if(slic>=25):
        divs=3
    if(slic>=30):
        divs=4
    if(slic>=35):
        divs=6
    if(slic>=40):
        divs=10
    if(slic>=45):
        divs=30
    if(slic>=46):
        divs=40


    peaks,_ = scp.signal.find_peaks(signal_filtered, prominence=pro)
    valleys,_ = scp.signal.find_peaks(-signal_filtered, prominence=pro)

    print(slic)

    #if(times[peaks[0]]<times[valleys[0]]):
    #    exp_damps[it] = (signal_filtered[valleys[n_peak]]-signal_filtered[peaks[n_peak]])/(times[valleys[n_peak]]-times[peaks[n_peak]])*2*np.pi
    #else:
    #    exp_damps[it] = (signal_filtered[valleys[2]]-signal_filtered[peaks[1]])/(times[valleys[2]]-times[peaks[1]])*2*np.pi

    func_peak = signal_filtered[peaks[0]]
    if(slic<=40):
        params, covs = opt.curve_fit(line, times[1:n_figs//divs], signal_filtered[1:n_figs//divs])
    else:
        params, covs = opt.curve_fit(line, times[np.argmax(signal_filtered):np.argmin(signal_filtered)], signal_filtered[np.argmax(signal_filtered):np.argmin(signal_filtered)])
    exp_damps[it] = params[0]

    if(slic==ptest):
        plt.plot(times, np.exp(signal_raw), label="raw signal")
        plt.plot(times, np.exp(signal_filtered), label="filtered signal")
        #plt.plot(times, np.cos(np.sqrt(ge * alpha * k_slic + 3 * beta / alpha * k_slic**2)*times))
        #plt.plot(times, signal_filtered[peaks[0]] + damping_rate(k_slic,T)*(times-times[peaks[0]]), c='r')
        plt.plot(times, np.exp(params[1] + exp_damps[it]*(times-times[peaks[n_peak]])), c='g', label="decay fit")
        #plt.scatter(times[peaks], signal_filtered[peaks], c='r', zorder=99)
        #plt.scatter(times[valleys], signal_filtered[valleys], c='orange', zorder = 100)
        plt.xlabel(r"$t \cdot \omega_p$")
        plt.ylabel(r"$\tilde{\rho}(k_x,t)$")
        #plt.ylim((-10,5))
        plt.yscale('log')
        plt.legend()
        plt.title(r"$k_x =$ "+str(slic*dk)+r" [$n_0^{-1/3}$]")
        plt.tight_layout()
        plt.savefig("./img/"+name+"_decay.png", dpi=200)
        plt.close()


th_damps = np.zeros(200)
for i in range(200):
    th_damps[i] = damping_rate(np.linspace(1,50,200,endpoint=True)[i]*dk,T)

plt.plot(np.linspace(1,50,200,endpoint=True)*dk, -th_damps)
plt.scatter(np.linspace(1,50,50,endpoint=True)*dk, -exp_damps, c='r',s=100,marker='+',linewidths=1)
plt.scatter(np.array([ptest])*dk, -exp_damps[ptest-1])
plt.xlim((2,5))
#plt.yscale('log')
plt.xlabel(r"$k_x\cdot n_0^{1/3}$")
plt.ylabel(r"$-\gamma / \omega_p$")
#plt.yscale('log')
plt.tight_layout()
plt.savefig("./img/"+name+"_landau.png", dpi=200)
plt.close()
