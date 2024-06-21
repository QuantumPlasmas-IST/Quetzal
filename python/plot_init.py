# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

#############################
name = 'noise'
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

# Reads Test_Data.h5 file in ../data/ directory
file=h5py.File('./output/'+name+'.h5','r')

species = file['/Species-1'][0]

pos_min = file.attrs['Position Min.']
pos_max = file.attrs['Position Max.']
mom_min = file.attrs['Momentum Min.']
mom_max = file.attrs['Momentum Max.']
pos_num = file.attrs['Number of Position Points']

sc = plt.imshow(species.T, extent = (pos_min[0]-0.05,pos_max[0]+0.05,mom_min[0]-0.05,mom_max[0]+0.05), aspect='auto', interpolation='quadric', origin = 'lower')
plt.colorbar(sc)
plt.tight_layout()
plt.savefig("./img/"+name+"_init.png", dpi=1000)