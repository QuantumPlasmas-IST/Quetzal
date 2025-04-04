# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap

#############################
name1 = 'twostreamtest1'
name2 = 'twostreamtest2'
iter = -1
axis = 2
mode = 'c' # s for sum; c for cut
plane1 = 0
plane2 = 3
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
file1=h5py.File('./output/'+name1+'.h5','r')
file2=h5py.File('./output/'+name2+'.h5','r')

electrons = file1['/Species'+str(iter)][0]
holes = file2['/Species'+str(iter)][0]
potential = file1['/Fields'+str(iter)][0]['real']-file2['/Fields'+str(iter)][0]['real']
#pot_im = file1['/Fields'+str(iter)][0]['imaginary']
#charge = file1['/Sources'+str(iter)][0]['real']-file2['/Sources'+str(iter)][0]['real']
#charge_im = file1['/Sources'+str(iter)][0]['imaginary']

species_l = electrons-holes

pos_min = file1.attrs['Position Min.']
pos_max = file1.attrs['Position Max.']
mom_min = file1.attrs['Momentum Min.']
mom_max = file1.attrs['Momentum Max.']
pos_dim = file1.attrs['Number of Position Points']
mom_dim = file1.attrs['Number of Momentum Points']
pos_num = file1.attrs['Number of Position Points']
mom_num = file1.attrs['Number of Momentum Points']
dx = file1.attrs['Position Delta']
dp = file1.attrs['Momentum Delta']
dt = file1.attrs['Time Delta']

shape = tuple(np.concatenate((pos_num,mom_num)))
dl = tuple(np.concatenate((dx,dp)))
maxs = tuple(np.concatenate((pos_max,mom_max)))
mins = tuple(np.concatenate((pos_min,mom_min)))

if mode == 's':
    if axis == 1:
        species = np.sum(species_l, axis=(2,3)) * dl[2] * dl[3]
    if axis == 2:
        species = np.sum(species_l, axis=(1,3)) * dl[1] * dl[3]
    if axis == 3:
        species = np.sum(species_l, axis=(1,2)) * dl[1] * dl[2]
if mode == 'c':
    if axis == 1:
        species = species_l[:, :, int(plane1/dl[2]), int(plane2/dl[3])]
    if axis == 2:
        species = species_l[:, int(plane1/dl[1]), :, int(plane2/dl[3])]
    if axis == 3:
        species = species_l[:, int(plane1/dl[1]), int(plane2/dl[2]), :]

fig, axs = plt.subplots(layout='constrained')

#print(np.shape(charge))


sc = plt.imshow(species.T, extent = (mins[0], maxs[0], mins[axis],maxs[axis]), aspect='auto', origin = 'lower' )
plt.tight_layout()
plt.xlabel(r"$x$ [$l_0$]")
plt.ylabel(r"$y$ [$l_0$]")
plt.title(r"$t = $ "+str(iter*dt)+r" $t_0$") 
plt.colorbar(sc, extend = 'both')
plt.savefig("./img/"+name1+"_diff_"+str(iter)+".png", dpi=500)
plt.close()


sc = plt.imshow(potential[:,:].T, extent = (pos_min[0],pos_max[0],pos_min[0],pos_max[0]), aspect='auto', origin = 'lower')
plt.xlabel(r"$x / a$")
plt.ylabel(r"$y / a$") 
plt.colorbar()
plt.tight_layout()
plt.savefig("./img/"+name1+"_diffField_"+str(iter)+".png", dpi=500)
plt.close()
