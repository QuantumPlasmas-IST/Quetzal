# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
import scipy.optimize as opt
import mpmath as mpm

#############################
name = 'weibel_vac'
T = 0.00005
u = 0.05
wpe = np.sqrt(2)
ge = np.pi
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
Nt = 100
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

#print("dx = ",dx)
#print("dt = ",dt)
#print("Dt = ",Dt)
#print("wmax = ",wmax)
#print("kmax = ",kmax)
#print("dw = ", dw)
#print("dk = ", dk)

for i in range(n_figs):

    potential = file['/Fields'+str((i+1)*diag_f)][1]['real']
    charge1 = file['/Sources'+str((i+1)*diag_f)][0]['real']
    #charge2 = file['/Sources'+str((i+1)*diag_f)][1]['real']

    phi[i,:] = potential[padding:-padding]
    rho[i,:] = charge1[padding:-padding]#-charge2[padding:-padding]
    #rho[i,:] = charge1[padding:-padding]

    #print(np.average(phi[i,:]))

Frho = scp.fft.fftshift(scp.fft.fft2(rho))
Fphi = scp.fft.fftshift(scp.fft.fft2(phi))

Srho = np.abs(Frho)**2
Sphi = np.abs(Fphi)**2

alpha = T * np.log(2+2*np.cosh(u/T))
beta = (T**2) * ((np.pi**2)/6 + (u/T)**2 / 2) 

qs = np.linspace(0,50,1000)

sc = plt.imshow(np.log10(Srho+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
plt.plot(qs, np.sqrt(1 * qs + 3 * T * qs**2), c = 'r') # 3D quadratic
#plt.plot(qs, np.sqrt(ge * alpha * qs + 3*beta/alpha * qs**2), c = 'r')     # 2D quadratic
#plt.plot(qs, np.sqrt(ge * alpha / 2 * qs + 3/4 * qs**2), c = 'r')    # 2D linear
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
#plt.plot(qs, np.sqrt(wpe**2 + 3 * T * qs**2), c = 'r') # 3D quadratic
#plt.plot(qs, np.sqrt(ge * alpha * qs + 3*beta/alpha * qs**2), c = 'r')     # 2D quadratic
#plt.plot(qs, np.sqrt(ge * alpha / 2 * qs + 3/4 * qs**2), c = 'r')    # 2D linear
#plt.xlim((-10,10))
#plt.ylim((0,6))
plt.xlabel(r"$k_x$[$\omega_{pe} c^{-1}$]")
plt.ylabel(r"$\omega$[$\omega_{pe}$]")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_fieldFFT.png", dpi=200)
plt.close()

"""kays = np.linspace(0,kmax/2,(pos_num-4)//2)

omega = np.argmax(Srho[n_figs//2: , (pos_num-4)//2:], axis = 0) * dw

def disp_rel(x, a, b):
    #return a*x+b*x*x
    #return np.pi*a*x+3*b/a*x*x-9*b*b/(a*a*a*np.pi)*x*x*x + 54*b*b*b/(a**5 * np.pi*np.pi)*x*x*x*x
    return np.pi/2*a*x + 3/4*b*b*x*x - 9 * b**4 * x**3/(8 * a * np.pi) + 27 * b**6 * x**4/(8 * a*a * np.pi**2)

limit=10
if(name == 'dispersion_2DFermiQuadratic0'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic1'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic2'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic3'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic4'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic5'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic6'):
    limit = 10
if(name == 'dispersion_2DFermiQuadratic7'):
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
plt.plot(kays[:limit], np.sqrt(disp_rel(kays[:limit],params[0],params[1])))
plt.plot(kays[:limit], np.sqrt(disp_rel(kays[:limit],alpha,beta)))
plt.xlim((0,5))
plt.tight_layout()
plt.savefig("./img/"+name+"_omega.png", dpi=200)
plt.close()

print(params[0])
print(params[1])
print()
print(alpha)
print(beta)
"""
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

"""
times = np.linspace(Dt,Nt*dt,n_figs, endpoint=True)

def fermi_dirac_derivative(func_p, func_u, func_T):
    #return - func_p / np.cosh((func_p*func_p/2-u)/(2*T))**2 /(4*T)
    return np.sqrt(2*np.pi/func_T)* func_p * mpm.polylog(-1/2,-np.exp(-(func_p*func_p/2-func_u)/func_T))


def damping_rate(func_k, func_u, func_T):
    func_alpha = func_T * np.log(2+2*np.cosh(func_u/func_T))
    func_beta = (func_T**2) * ((np.pi**2)/6 + (func_u/func_T)**2 / 2) 
    func_w = np.sqrt(ge * func_alpha * func_k + 3 * func_beta / func_alpha * func_k**2)
    func_p = func_w/func_k
    return np.pi * func_w / (4 * func_k) * (fermi_dirac_derivative(func_p, func_u, func_T)+fermi_dirac_derivative(func_p, -func_u, func_T))

exp_damps = np.zeros(50)

ptest=15

for it in range(30):

    slic = it+1
    k_slic = dk*slic

    signal_raw = np.log(Drho[:,(pos_num-4)//2+slic])
    signal_raw = signal_raw-np.average(signal_raw)

    divs = 2

    if(slic>0):
        divs=10
        cutoff=0.005+0.003*slic
    if(slic>11):
        divs=20
    if(slic>14):
        divs=30
    if(slic>21):
        divs=40
        cutoff=0.05

    b, a = scp.signal.butter(5, cutoff, btype='low')
    signal_raw[0] = signal_raw[1]
    signal_filtered = scp.signal.medfilt(scp.signal.lfilter(b, a, signal_raw), 1)

    #exp_damps[it] = (signal_filtered[valleys[n_valley]]-signal_filtered[peaks[n_peak]])/(times[valleys[n_valley]]-times[peaks[n_peak]])
    
    cond_min = np.argmax(signal_filtered[:n_figs//divs])
    cond_max = np.argmin(signal_filtered[cond_min:n_figs//divs])
    time_min = times[cond_min]
    time_max = times[cond_max]

    if(slic<0):
        params, covs = opt.curve_fit(line, times[1:cond_max], signal_filtered[1:cond_max])
    else:
        params, covs = opt.curve_fit(line, times[cond_min:cond_max], signal_filtered[cond_min:cond_max])
    
    if(slic<5):
        params, covs = opt.curve_fit(line, times[1:cond_max], signal_filtered[1:cond_max])
    else:
        peaks,_ = scp.signal.find_peaks(signal_raw)
        peaks2,_= scp.signal.find_peaks(-signal_raw[peaks])
        #if(slic==ptest):
        #    print(peaks2)
        params, covs = opt.curve_fit(line, times[peaks[:peaks2[0]-1]], signal_raw[peaks[:peaks2[0]-1]])
    
    exp_damps[it] = params[0]


    if(slic==ptest):
        plt.plot(times, np.exp(signal_raw), label="raw signal")
        plt.plot(times, np.exp(signal_filtered), label="filtered signal")
        #plt.plot(times, np.cos(np.sqrt(ge * alpha * k_slic + 3 * beta / alpha * k_slic**2)*times))
        #plt.plot(times, np.exp((params[1] + damping_rate(k_slic,u,T)*times).astype(float)), c='r')
        plt.plot(times, np.exp(params[1] + exp_damps[it]*times), c='g', label="decay fit")
        plt.scatter(times[peaks], np.exp(signal_raw[peaks]), c='r', zorder=99,s=10)
        #plt.vlines([times[peaks[peaks2[0]]]],1e-10,1e10)
        #plt.scatter(times[valleys], np.exp(signal_filtered[valleys]), c='orange', zorder = 100)
        plt.ylim((1e-4,2e2))
        plt.xlim((0,300))
        plt.yscale('log')
        plt.xlabel(r"$t/t_0$")
        plt.ylabel(r"$\tilde{\rho}(k_x,t)$")
        plt.title(r"$k_x =$ "+str(slic*dk)+r" [$l_0^{-1}$]")
        plt.legend()
        plt.tight_layout()
        plt.savefig("./img/"+name+"_decay.png", dpi=200)
        plt.close()

    


th_damps = np.zeros(200)
for i in range(200):
    th_damps[i] = damping_rate(np.linspace(1,50,200,endpoint=True)[i]*dk,u,T)

plt.plot(np.linspace(1,50,200,endpoint=True)*dk, -th_damps)
plt.scatter(np.linspace(1,50,50,endpoint=True)*dk, -exp_damps,c='r',s=100,marker='+',linewidths=1)
plt.scatter(np.array([ptest])*dk, -exp_damps[ptest-1])
#plt.yscale('log')
plt.xlim((0,3))
plt.ylim((-0.02,0.2))
plt.xlabel(r"$k_x\cdot l_0$")
plt.ylabel(r"$-\gamma\cdot t_0$")
plt.tight_layout()
plt.savefig("./img/"+name+"_landau.png", dpi=200)
plt.close()
"""