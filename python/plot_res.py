# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

#############################
name = 'anomalous'
iter = 0
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

electrons = file['/Species'+str(iter)][0]
positrons = file['/Species'+str(iter)][1]
species = positrons-electrons
potential = file['/Fields'+str(iter)][0]['real']
pot_im = file['/Fields'+str(iter)][0]['imaginary']
charge = file['/Sources'+str(iter)][0]['real']
charge_im = file['/Sources'+str(iter)][0]['imaginary']

pos_min = file.attrs['Position Min.'][0]
pos_max = file.attrs['Position Max.'][0]
mom_min = file.attrs['Momentum Min.'][0]
mom_max = file.attrs['Momentum Max.'][0]
pos_num = file.attrs['Number of Position Points'][0]
dx = file.attrs['Position Delta'][0]
dt = file.attrs['Time Delta']
padding=2

fig, axs = plt.subplots(layout='constrained')
sc = plt.imshow(species.T, extent = (pos_min-0.05,pos_max+0.05,mom_min-0.05,mom_max+0.05), aspect='auto', origin = 'lower')
#plt.tight_layout()
#plt.plot(np.linspace(pos_min,pos_max,pos_num), species[:,50], c='b')
print(np.shape(species))
plt.scatter(np.linspace(pos_min,pos_max,pos_num), scale*potential, c='r', s=0.2)
plt.plot(np.linspace(pos_min,pos_max,pos_num), scale*charge, c='y', linewidth=1)
secay = axs.secondary_yaxis('right', functions=(func1, func2))
plt.xlabel(r"$x$ [$c \omega_{pe}^{-1}$]")
plt.ylabel(r"$p$ [$m_e c$]")
plt.title(r"$t = $ "+str(iter*dt)+r" $\omega_{pe}^{-1}$") 
plt.colorbar(sc)
#plt.gca().set_aspect('equal')
plt.savefig("./img/"+name+"_"+str(iter)+".png", dpi=500)

#print(np.amax(potential))
#print(np.argmax(potential))
print("////////////////")
point = int((2*pos_num)//5)
point=pos_num//2
#print(-charge[point])
#print((potential[point+1]-2*potential[point]+potential[point-1])/dx**2)