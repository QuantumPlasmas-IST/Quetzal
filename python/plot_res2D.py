# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap

#############################
name = 'weibel3'
iter = 8000
supress = 0
mode = 's' # s for sum; c for cut
plane = 10
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
#holes = file['/Species'+str(iter)][1]
potential = file['/Fields'+str(iter)][1]['real']
pot_im = file['/Fields'+str(iter)][0]['imaginary']
charge = file['/Sources'+str(iter)][0]['real']
charge_im = file['/Sources'+str(iter)][0]['imaginary']

species_l = electrons

pos_min = file.attrs['Position Min.']
pos_max = file.attrs['Position Max.']
mom_min = file.attrs['Momentum Min.']
mom_max = file.attrs['Momentum Max.']
pos_dim = file.attrs['Number of Position Points']
mom_dim = file.attrs['Number of Momentum Points']
pos_num = file.attrs['Number of Position Points']
mom_num = file.attrs['Number of Momentum Points']
dx = file.attrs['Position Delta']
dp = file.attrs['Momentum Delta']
dt = file.attrs['Time Delta']

shape = tuple(np.concatenate((pos_num,mom_num)))
dl = tuple(np.concatenate((dx,dp)))
maxs = tuple(np.concatenate((pos_max,mom_max)))
mins = tuple(np.concatenate((pos_min,mom_min)))

if mode == 's':
    species = np.sum(species_l, axis=supress) * dl[supress]
if mode == 'c':
    species = species_l[:,int(plane/dl[supress]),:]

fig, axs = plt.subplots(layout='constrained')

print(np.shape(charge))

if(supress == -1):
    sc = plt.imshow(species.T, extent = (mins[0]-2*dl[0], maxs[0]+dl[0], mins[1]-2*dl[1],maxs[1]+dl[1]), aspect='auto', origin = 'lower', cmap=my_cmap ,vmin=-np.amax(np.abs(species))*0.8,vmax=np.amax(np.abs(species))*0.8)
    #sc = plt.imshow(charge.T, extent = (mins[0], maxs[0], mins[1],maxs[1]), aspect='auto', origin = 'lower', cmap='RdBu' ,vmin=-np.amax(np.abs(charge))*0.8,vmax=np.amax(np.abs(charge))*0.8)
    plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), scale*potential, c='r')
    #secay = axs.secondary_yaxis('right', functions=(func1, func2))
if(supress == 0):
    sc = plt.imshow(species.T, extent = (mins[1]-2*dl[1], maxs[1]+dl[1], mins[2]-2*dl[2],maxs[2]+dl[2]), aspect='auto', origin = 'lower', cmap='RdBu',vmin=-np.amax(np.abs(species)),vmax=np.amax(np.abs(species)))
if(supress == 1):
    sc = plt.imshow(species.T, extent = (mins[0], maxs[0], mins[2], maxs[2]), aspect='auto', origin = 'lower', cmap='RdBu',vmin=-np.amax(np.abs(species)),vmax=np.amax(np.abs(species)))
    plt.plot(np.linspace(pos_min[0],pos_max[0],pos_num[0]), scale*potential, c='r')
    secay = axs.secondary_yaxis('right', functions=(func1, func2))

#plt.tight_layout()
plt.xlabel(r"$x$ [$l_0$]")
plt.ylabel(r"$y$ [$l_0$]")
plt.title(r"$t = $ "+str(iter*dt)+r" $t_0$") 
plt.colorbar(sc, extend = 'both')
#plt.gca().set_aspect('equal')
plt.savefig("./img/"+name+"_"+str(iter)+".png", dpi=500)
