# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

#############################
name = 'landau2'
iter =70000
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

# Reads Test_Data.h5 file in ../data/ directory
file=h5py.File('./output/'+name+'.h5','r')

species = file['/Species'+str(iter)][0]
potential = file['/Fields'+str(iter)][0]['real']
charge = file['/Sources'+str(iter)][0]['real']


pos_min = file.attrs['Position Min.']
pos_max = file.attrs['Position Max.']
mom_min = file.attrs['Momentum Min.']
mom_max = file.attrs['Momentum Max.']
pos_num = file.attrs['Number of Position Points']
mom_num = file.attrs['Number of Momentum Points']
dt = file.attrs['Time Delta']


potential = (potential-np.amin(potential)) / (np.amax(potential)-np.amin(potential)) * (mom_max[0]-mom_min[0]) + mom_min[0]
charge = (charge-np.amin(charge)) / (np.amax(charge)-np.amin(charge)) * (mom_max[0]-mom_min[0]) + mom_min[0]

#sc = plt.imshow(species.T, extent = (pos_min[0]-0.05,pos_max[0]+0.05,mom_min[0]-0.05,mom_max[0]+0.05), aspect='auto', interpolation='quadric', origin = 'lower')
plt.plot(np.linspace(mom_min[0],mom_max[0],mom_num[0]), np.sum(species, axis = 0)/pos_num, c='r')
#plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), potential, c='r')
#plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), charge, c='y')
#plt.gca().set_aspect('equal')
plt.xlabel(r"$p$ [$m_e c$]")
plt.ylabel(r"$f$")
plt.title(r"$t = $ "+str(iter*dt)+r" $\omega_{pe}^{-1}$") 
plt.tight_layout()
plt.savefig("./img/"+name+"_mom.png", dpi=100)

print(np.amax(species))
print(np.amin(species))
print(np.argmax(species))
print(np.argmin(species))