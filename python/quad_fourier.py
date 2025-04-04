# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
import scipy.optimize as opt
import mpmath as mpm

#############################
name = 'dispersion_2DFerm1'
T = 0.01
u = 0.001
ge = np.pi

wp2 = 1
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 18})

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

#print("dx = ",dx)
#print("dt = ",dt)
#print("Dt = ",Dt)
#print("wmax = ",wmax)
#print("kmax = ",kmax)
#print("dw = ", dw)
#print("dk = ", dk)

for i in range(n_figs):

    potential = file['/Fields'+str((i+1)*diag_f)][0]['real']
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

gamma = 1/2*T*(mpm.polylog(2,-np.exp(u/T))/mpm.polylog(1,-np.exp(u/T))+mpm.polylog(2,-np.exp(-u/T))/mpm.polylog(1,-np.exp(-u/T)))
#gamma = 1/2*T*(mpm.polylog(2,-np.exp(u/T))/mpm.polylog(1,-np.exp(u/T))+1)


gamma1 = T*mpm.polylog(2,-np.exp(u/T))/mpm.polylog(1,-np.exp(u/T))

qs = np.linspace(0,50,1000)

sc = plt.imshow(np.log10(Srho+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
#plt.plot(qs, np.sqrt(wp2 * qs + 3 * T * qs**2), c = 'r')                  # 2D Maxwell
#plt.plot(qs, np.sqrt(qs + 3*gamma * qs**2), c = 'r')     # 2D Fermi
plt.plot(qs, qs, c = 'r')     # Damping threshold
plt.xlim((-5,5))
plt.ylim((0,np.pi))
plt.clim((-5,2))
plt.xlabel(r"$k_x \cdot l_0$")
plt.ylabel(r"$\omega$[$\omega_{pe}$]")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_chargeFFT.png", dpi=200)
plt.close()

sc = plt.imshow(np.log10(Sphi+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,-wmax/2-dw/2,wmax/2-dw/2), aspect='auto', origin = 'lower')
#plt.plot(-qs, np.sqrt(qs/2 + np.sqrt((qs/2)**2 + 3*gamma*qs**3)), c = 'g')                    # 2D Maxwell
plt.plot(qs, np.sqrt(qs + 3*gamma * qs**2), c = 'r')     # 2D Fermi
#plt.plot(qs, qs, c = 'r')     # Damping threshold
plt.xlim((-10,10))
plt.ylim((0,4))
plt.xlabel(r"$k_x \cdot l_0$")
plt.ylabel(r"$\omega\cdot t_0$")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_fieldFFT.png", dpi=200)
plt.close()

kays = np.linspace(0,kmax/2,(pos_num-4)//2)

omega = np.argmax(Srho[n_figs//2: , (pos_num-4)//2:], axis = 0) * dw

def disp_rel(x, a, b):
    #return a*x+b*x*x
    return a*x+3*b*x*x
    #return np.pi/2*a*x + 3/4*b*b*x*x - 9 * b**4 * x**3/(8 * a * np.pi) + 27 * b**6 * x**4/(8 * a*a * np.pi**2)

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
    limit = 8

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
plt.plot(kays[:limit], np.sqrt(disp_rel(kays[:limit],params[0],params[1])),zorder=1000)
plt.plot(kays[:2*limit], np.sqrt(disp_rel(kays[:2*limit],1,gamma)))
plt.xlim((0,5))
plt.tight_layout()
plt.savefig("./img/"+name+"_omega.png", dpi=200)
plt.close()

print("alpha =", params[0],"\pm", cov[0,0], cov[0,1])
print("gamma =", params[1],"\pm", cov[1,1], cov[1,0])
print()
print(gamma)
fit_alpha = params[0]
fit_beta = params[1]

Drho = np.abs(scp.fft.fftshift(scp.fft.fft(rho,axis=-1),axes=-1))
Dphi = np.abs(scp.fft.fftshift(scp.fft.fft(phi,axis=-1),axes=-1))

sc = plt.imshow(np.log(Drho+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,0,Nt*dt), aspect='auto', origin = 'lower')
plt.xlim((-15,15))
plt.ylim((0,500))
plt.xlabel(r"$k_x\cdot l_0$")
plt.ylabel(r"$t/t_0$")
plt.clim((-6,-3))
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_chargeDamp.png", dpi=200)
plt.close()

sc = plt.imshow(np.log(Dphi+1e-5), extent = (-kmax/2-dk/2,kmax/2-dk/2,0,Nt*dt), aspect='auto', origin = 'lower')
plt.xlim((-10,10))
plt.ylim((0,500))
plt.xlabel(r"$k_x\cdot l_0$")
plt.ylabel(r"$t/t_0$")
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_fieldDamp.png", dpi=200)
plt.close()

times = np.linspace(Dt,Nt*dt,n_figs, endpoint=True)

def fermi_dirac_derivative(func_p, func_u, func_T):
    #return - func_p / np.cosh((func_p*func_p/2-u)/(2*T))**2 /(4*T)
    return np.sqrt(2*np.pi/func_T)* func_p * mpm.polylog(-1/2,-np.exp(-(func_p*func_p/2-func_u)/func_T))/(2*np.pi*func_T*mpm.polylog(-1/2,-np.exp(func_u/func_T)))


def damping_rate(func_k, func_u, func_T):
    func_alpha = 1+0*func_T * np.log(2+2*np.cosh(func_u/func_T))
    #func_beta = (func_T**2) * ((np.pi**2)/6 + (func_u/func_T)**2 / 2) 
    func_w = np.sqrt(func_k+3*0.0118*func_k**2)
    func_p = func_w/func_k
    return np.pi * func_w / (4*func_k) * float(-mpm.fabs(fermi_dirac_derivative(func_p, func_u, func_T)+fermi_dirac_derivative(func_p, -func_u, func_T)))

def maxwell_derivative(func_p, func_T):
    return -2/np.sqrt(2*np.pi*func_T) * func_p/func_T * np.exp(-func_p*func_p/(2*func_T))

def damping_rate_max(func_k, func_T):
    func_omegaP = 1
    func_beta = 3*func_T
    func_w = np.sqrt(func_omegaP * func_k + func_beta * func_k**2)
    func_p = func_w/func_k
    return np.pi * func_w / (4 * func_k) * maxwell_derivative(func_p,func_T)

landau_size=90
exp_damps = np.zeros(landau_size)

ptest=70

for it in range(0,90):

    print(it)

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
    if(slic>16):
        divs=50
        cutoff=0.05

    b, a = scp.signal.butter(5, cutoff, btype='low')
    signal_raw[0] = signal_raw[1]
    signal_filtered = scp.signal.medfilt(scp.signal.lfilter(b, a, signal_raw), 1)

    #exp_damps[it] = (signal_filtered[valleys[n_valley]]-signal_filtered[peaks[n_peak]])/(times[valleys[n_valley]]-times[peaks[n_peak]])
    
    cond_min = np.argmax(signal_filtered[:n_figs//divs])
    cond_max = np.argmin(signal_filtered[cond_min:n_figs//divs])
    time_min = times[cond_min]
    time_max = times[cond_max]
    
    pro=0.5
    if(slic>50):
        pro=0.8
    if(slic>65):
        pro=1
    if(slic>75):
        pro=0.8
    peaks,_ = scp.signal.find_peaks(signal_raw, prominence = (pro))
    peaks2,_= scp.signal.find_peaks(-signal_raw[peaks])
    #peaks = np.append([0],peaks)
    threshold = 50
    if(slic>=50):
        threshold=30
    if(slic>=60):
        threshold=20
    if(slic>=65):
        nPeaks=6
    if(slic>=74):
        nPeaks=5
    if(slic>=80):
        nPeaks=4

    if(slic<65):
        params, covs = opt.curve_fit(line, times[peaks[times[peaks]<threshold]], signal_raw[peaks[times[peaks]<threshold]])
    else:
        params, covs = opt.curve_fit(line, times[peaks[:nPeaks]], signal_raw[peaks[:nPeaks]])
    
    exp_damps[it] = params[0]


    if(slic==ptest):
        plt.plot(times, np.exp(signal_raw), label="raw signal")
        #plt.plot(times, np.exp(signal_filtered), label="filtered signal")
        #plt.plot(times, np.cos(np.sqrt(ge * alpha * k_slic + 3 * beta / alpha * k_slic**2)*times))
        #plt.plot(times, np.exp((params[1] + damping_rate(k_slic,u,T)*times).astype(float)), c='r')
        plt.plot(times, np.exp(params[1] + exp_damps[it]*times), c='g', label="decay fit")
        plt.scatter(times[peaks], np.exp(signal_raw[peaks]), c='r', zorder=99,s=10)
        #plt.vlines([threshold],[0.000001],[10000000])
        #plt.vlines([times[peaks[peaks2[0]]]],1e-10,1e10)
        #plt.scatter(times[valleys], np.exp(signal_filtered[valleys]), c='orange', zorder = 100)
        plt.ylim((2e-3,2e2))
        plt.xlim((0,50))
        plt.yscale('log')
        plt.xlabel(r"$t/t_0$")
        plt.ylabel(r"$\tilde{\rho}(k_x,t)$")
        plt.title(r"$k_x\cdot l_0 =$ "+str(np.round(slic*dk,2)))
        plt.legend()
        plt.tight_layout()
        plt.savefig("./img/"+name+"_decay.png", dpi=200)
        plt.close()

    


th_damps = np.zeros(landau_size*20)
p_theo = np.linspace(1,landau_size*20,landau_size*20,endpoint=True)
p_exp = np.linspace(1,landau_size,landau_size,endpoint=True)
for i in range(landau_size*2):
    #th_damps[i] = damping_rate(np.linspace(1,landau_size,landau_size,endpoint=True)[i]*dk,u,T)
    th_damps[i] = damping_rate(p_theo[i]*dk,u,T)

plt.plot(p_theo*dk-2.2, -th_damps, linewidth = 3)
plt.scatter(p_exp*dk, -exp_damps,c='r',s=150,marker='+',linewidths=1, label = r"$\mu = 0.001$ $m_0 (l_0/t_0)^2$", zorder=10)
#plt.scatter(p_exp[ptest-1]*dk, -exp_damps[ptest-1])
plt.scatter(np.array([1000])*dk, -exp_damps[ptest-1], label = r"$T = 0.01$ $m_0 (l_0/t_0)^2$")
#plt.yscale('log')
plt.xlim((2,9))
plt.ylim((-0.02,0.5))
plt.xlabel(r"$k_x\cdot l_0$")
plt.ylabel(r"$-\gamma\cdot t_0$")
plt.legend(markerscale=0)
plt.tight_layout()
plt.savefig("./img/"+name+"_landau.png", dpi=400)
plt.close()
