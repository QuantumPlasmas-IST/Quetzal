# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

#############################
name = 'weibel_vac'
iter = 1
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

# Set the scaling for the potential and charge density plots
scale = 1

def func1(x):
    return x/scale
def func2(x):
    return x*scale

# Reads Test_Data.h5 file in ../data/ directory
file=h5py.File('./output/'+name+'.h5','r')

potential = file['/Fields'+str(iter)][1]['real']
pot_im = file['/Fields'+str(iter)][1]['imaginary']
potential2 = file['/Fields'+str(iter)][2]['real']
pot_im2 = file['/Fields'+str(iter)][2]['imaginary']

#charge = file['/Sources'+str(iter)][0]['real']
#charge_im = file['/Sources'+str(iter)][0]['imaginary']

pos_min = file.attrs['Position Min.'][0]
pos_max = file.attrs['Position Max.'][0]
mom_min = file.attrs['Momentum Min.'][0]
mom_max = file.attrs['Momentum Max.'][0]
pos_num = file.attrs['Number of Position Points'][0]
dx = file.attrs['Position Delta'][0]
dt = file.attrs['Time Delta']

kernel = np.zeros(512)
for i in range(253,259):
    kernel[i] = 1
#potential2 = np.convolve(charge,kernel,mode='same')

#plt.scatter(np.linspace(pos_min,pos_max,pos_num, endpoint=False), potential, c='r', s=1)
#plt.scatter(np.linspace(pos_min,pos_max,pos_num, endpoint=False), pot_im, c='b', s=1)
plt.plot(np.linspace(pos_min,pos_max,pos_num, endpoint=False), potential, c='r')
plt.plot(np.linspace(pos_min,pos_max,pos_num, endpoint=False), potential2, c='b')
plt.xlabel(r"$x$ [$c \omega_{pe}^{-1}$]")
plt.title(r"$t = $ "+str(iter*dt)+r" $\omega_{pe}^{-1}$") 
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_"+str(iter)+".png", dpi=500)
plt.close()

#print(np.amax(potential))
#print(np.argmax(potential))
print("////////////////")
point = int((2*pos_num)//5)
point=pos_num//2
#print(-charge[point])
#print((potential[point+1]-2*potential[point]+potential[point-1])/dx**2)

num=1024
ps = np.linspace(-55,55,num,endpoint=False)
Ey0 = scp.fft.ifft(scp.fft.fft(np.exp(-ps*ps/0.05)))
Bz0 = scp.fft.ifft(scp.fft.fft(0.00*np.sin(0.57119866 * ps)))
qs = 2*np.pi*scp.fft.fftfreq(num, 110./num)
I = 1j
dt = 0.1
fEy = np.copy(scp.fft.fft(Ey0))
fBz = np.copy(scp.fft.fft(Bz0))

for i in range(iter):
    # 1
    fBz1 = np.copy(fBz - dt*I*qs*fEy)
    fEy1 = np.copy(fEy - dt*I*qs*fBz)

    # 2
    fBz2 = np.copy(fBz1 - dt*I*qs*fEy1)
    fEy2 = np.copy(fEy1 - dt*I*qs*fBz1)

    # Av
    fEy3 = np.copy((fEy + fEy2)/2.)
    fBz3 = np.copy((fBz + fBz2)/2.)

    fEy = np.copy(fEy3)
    fBz = np.copy(fBz3)

Bz = np.copy(scp.fft.ifft(fBz))
Ey = np.copy(scp.fft.ifft(fEy))

plt.plot(ps, Ey, c='r')
plt.plot(ps, Bz, c='b')
plt.xlabel(r"$x$ [$c \omega_{pe}^{-1}$]")
plt.title(r"$t = $ "+str(iter*dt)+r" $\omega_{pe}^{-1}$") 
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_"+str(iter)+"_theo.png", dpi=500)
plt.close()

plt.scatter(scp.fft.fftshift(qs), np.real(scp.fft.fftshift(fEy)), c='r', s=1)
plt.scatter(scp.fft.fftshift(qs), np.imag(scp.fft.fftshift(fEy)), c='b', s=1)
plt.xlabel(r"$x$ [$c \omega_{pe}^{-1}$]")
plt.title(r"$t = $ "+str(iter*dt)+r" $\omega_{pe}^{-1}$") 
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_"+str(iter)+"_fftE.png", dpi=500)
plt.close()

plt.scatter(scp.fft.fftshift(qs), np.real(scp.fft.fftshift(fBz)), c='r', s=1)
plt.scatter(scp.fft.fftshift(qs), np.imag(scp.fft.fftshift(fBz)), c='b', s=1)
plt.xlabel(r"$x$ [$c \omega_{pe}^{-1}$]")
plt.title(r"$t = $ "+str(iter*dt)+r" $\omega_{pe}^{-1}$") 
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_"+str(iter)+"_fftB.png", dpi=500)
plt.close()