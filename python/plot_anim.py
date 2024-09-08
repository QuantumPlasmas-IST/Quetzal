# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
from celluloid import Camera

#############################
name = 'anomalousF'
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

# Creates movie
fig = plt.figure()
camera = Camera(fig)

# Reads Test_Data.h5 file in ../data/ directory
file=h5py.File('./output/'+name+'.h5','r')

pos_min = file.attrs['Position Min.']
pos_max = file.attrs['Position Max.']
mom_min = file.attrs['Momentum Min.']
mom_max = file.attrs['Momentum Max.']
Nt = file.attrs['Number of Timesteps']
#Nt=10050
diag_f = file.attrs['Full Diagnostic Frequency']
pos_num = file.attrs['Number of Position Points']
n_figs = (Nt)//diag_f

for i in range(n_figs):

    result = file['/Species'+str(i*diag_f)]
    positrons = result[0]
    electrons = result[1]
    potential = file['/Fields'+str(i*diag_f)][0]['real']
    charge = file['/Sources'+str(i*diag_f)][0]['real']

    #potential = (potential-np.amin(potential)) / (np.amax(potential)-np.amin(potential)) * (mom_max[0]-mom_min[0]) + mom_min[0]
    #charge = (charge-np.amin(charge)) / (np.amax(charge)-np.amin(charge)) * (mom_max[0]-mom_min[0]) + mom_min[0]

    print("timestep "+str(i*diag_f))

    sc = plt.imshow((positrons-electrons).T, extent = (pos_min[0]-0.05,pos_max[0]+0.05,mom_min[0]-0.05,mom_max[0]+0.05), aspect='auto', origin = 'lower')
    plt.xlabel(r"$x$ [$c \omega_{pe}^{-1}$]")
    plt.ylabel(r"$p$ [$m_e c$]")
    plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), potential, c='r')
    plt.ylim((mom_min[0]-0.05,mom_max[0]+0.05))
    #plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), 4*charge, c='y')
    #plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), 0*np.linspace(pos_min[0],pos_max[0],pos_num[0]))
    #plt.colorbar(sc)
    #plt.gca().set_aspect('equal')
    plt.tight_layout()
    camera.snap()

animation = camera.animate(interval=50)
animation_name = './img/'+name+'.mp4'
animation.save(animation_name)