# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

#############################
#name = 'weibel_linBandNew2'
#iter = 1000
name = 'weibel_linBandWall1'
iter = 10000
#name = 'wallTest'
#iter = 2200
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

potential = file['/Fields'+str(iter)][2]['real']
pot_im = file['/Fields'+str(iter)][0]['imaginary']
charge = file['/Sources'+str(iter)][0]['real']#- file['/Sources'+str(iter)][3]['real']
charge_im = file['/Sources'+str(iter)][0]['imaginary']

pos_min = file.attrs['Position Min.'][0]
pos_max = file.attrs['Position Max.'][0]
mom_min = file.attrs['Momentum Min.'][0]
mom_max = file.attrs['Momentum Max.'][0]
pos_num = file.attrs['Number of Position Points'][0]
dx = file.attrs['Position Delta'][0]
dt = file.attrs['Time Delta']

#sc = plt.imshow(charge.T, extent = (pos_min,pos_max,mom_min,mom_max), aspect='auto', origin = 'lower')
sc = plt.imshow(potential.T, extent = (pos_min,pos_max,mom_min,mom_max), aspect='auto', origin = 'lower')
plt.xlabel(r"$x / a$")
plt.ylabel(r"$y / a$") 
plt.title("t =" + str(dt*iter))
plt.colorbar()
plt.tight_layout()
plt.savefig("./img/"+name+"_charge_"+str(iter)+".png", dpi=500)
plt.close()