# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
from celluloid import Camera


#############################
name = 'landau'
iter = 29950
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

fig = plt.figure()
camera = Camera(fig)

# Reads Test_Data.h5 file in ../data/ directory
file=h5py.File('./output/'+name+'.h5','r')

pos_min = file.attrs['Position Min.']
pos_max = file.attrs['Position Max.']
mom_min = file.attrs['Momentum Min.']
mom_max = file.attrs['Momentum Max.']
pos_num = file.attrs['Number of Position Points']
mom_num = file.attrs['Number of Momentum Points']
Nt = file.attrs['Number of Timesteps']
diag_f = file.attrs['Diagnostic Frequency']
pos_num = file.attrs['Number of Position Points']
n_figs = (Nt)//diag_f

for i in range(n_figs):
    species = file['/Species'+str(i*diag_f)][0]
    potential = file['/Fields'+str(i*diag_f)][0]['real']
    charge = file['/Sources'+str(i*diag_f)][0]['real']

    #potential = (potential-np.amin(potential)) / (np.amax(potential)-np.amin(potential)) * (mom_max[0]-mom_min[0]) + mom_min[0]
    #charge = (charge-np.amin(charge)) / (np.amax(charge)-np.amin(charge)) * (mom_max[0]-mom_min[0]) + mom_min[0]

    #sc = plt.imshow(species.T, extent = (pos_min[0]-0.05,pos_max[0]+0.05,mom_min[0]-0.05,mom_max[0]+0.05), aspect='auto', interpolation='quadric', origin = 'lower')
    plt.plot(np.linspace(mom_min[0],mom_max[0],mom_num[0]), np.sum(species, axis = 0)/pos_num, c='r')
    #plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), potential, c='r')
    #plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), charge, c='y')
    #plt.gca().set_aspect('equal')
    plt.tight_layout()
    #plt.savefig("./img/"+name+"_mom.png", dpi=1000)
    camera.snap()
    print("timestep "+str(i*diag_f))

animation = camera.animate(interval=50)
animation_name = './img/'+name+'_mom.mp4'
animation.save(animation_name)