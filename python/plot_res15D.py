# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap

#############################
name = 'Ze2'
iter = 599
axis = 2
mode = 'c' # s for sum; c for cut
plane = 0
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
holes = file['/Species'+str(iter)][1]
potential = file['/Fields'+str(iter)][0]['real']
pot_im = file['/Fields'+str(iter)][0]['imaginary']
charge = file['/Sources'+str(iter)][0]['real']
charge_im = file['/Sources'+str(iter)][0]['imaginary']

species_l = electrons#-holes

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
    if axis == 0:
        species = np.sum(species_l, axis=(2,3)) * dl[2] * dl[3]
    if axis == 1:
        species = np.sum(species_l, axis=(1,3)) * dl[1] * dl[3]
    if axis == 2:
        species = np.sum(species_l, axis=(1,2)) * dl[1] * dl[2]
if mode == 'c':
    if axis == 0:
        species = species_l[:, :, int(plane1/dl[2]), int(plane2/dl[3])]
    if axis == 1:
        species = species_l[:, int(plane1/dl[1]), :, int(plane2/dl[3])]
    if axis == 2:
        species = species_l[:, int(plane1/dl[1]), int(plane2/dl[2]), :]

fig, axs = plt.subplots(layout='constrained')

print(np.shape(charge))


sc = plt.imshow(species.T, extent = (mins[0], maxs[0], mins[axis],maxs[axis]), aspect='auto', origin = 'lower' )
plt.tight_layout()
plt.xlabel(r"$x$ [$l_0$]")
plt.ylabel(r"$y$ [$l_0$]")
plt.title(r"$t = $ "+str(iter*dt)+r" $t_0$") 
plt.colorbar(sc, extend = 'both')
plt.savefig("./img/"+name+"_"+str(iter)+".png", dpi=500)
