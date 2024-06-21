# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp
from celluloid import Camera
from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap

img = imread('../Pictures/Screenshots/colormap.png')
print(img.shape)
colors_from_img = img[::-1, 0, :]
my_cmap = LinearSegmentedColormap.from_list('my_cmap', colors_from_img, N=878)

#############################
name = 'anomalous_warm'
supress = -1
#############################

# Allows the use of LateX notation in labels
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})

# Creates movie
fig = plt.figure()
camera = Camera(fig)

# Reads Test_Data.h5 file in ../data/ directory
file=h5py.File('./output/'+name+'.h5','r')

Nt = file.attrs['Number of Timesteps']
diag_f = file.attrs['Diagnostic Frequency']
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

n_figs = (Nt)//diag_f

for i in range(n_figs):

    result1 = file['/Species'+str(i*diag_f)][0]
    result2 = file['/Species'+str(i*diag_f)][1]
    species_l = np.reshape(result1-result2, shape)
    species = np.sum(species_l, axis=supress) * dl[supress]
    potential = file['/Fields'+str(i*diag_f)][0]['real']
    charge = file['/Sources'+str(i*diag_f)][0]['real']

    #potential = (potential-np.amin(potential)) / (np.amax(potential)-np.amin(potential)) * (mom_max[0]-mom_min[0]) + mom_min[0]
    #charge = (charge-np.amin(charge)) / (np.amax(charge)-np.amin(charge)) * (mom_max[0]-mom_min[0]) + mom_min[0]

    if(supress == -1):
        sc = plt.imshow(species.T, extent = (mins[0]-2*dl[0], maxs[0]+dl[0], mins[1]-2*dl[1],maxs[1]+dl[1]), aspect='auto', origin = 'lower', cmap='seismic',vmin=-np.amax(np.abs(species)),vmax=np.amax(np.abs(species)))
    if(supress == 0):
        sc = plt.imshow(species.T, extent = (mins[1]-2*dl[1], maxs[1]+dl[1], mins[2]-2*dl[2],maxs[2]+dl[2]), aspect='auto', origin = 'lower', cmap=my_cmap,vmin=-np.amax(np.abs(species)),vmax=np.amax(np.abs(species)))
    if(supress == 1):
        sc = plt.imshow(species.T, extent = (mins[0]-2*dl[0], maxs[0]+dl[0], mins[2]-2*dl[2],maxs[2]+dl[2]), aspect='auto', origin = 'lower', cmap=my_cmap,vmin=-np.amax(np.abs(species)),vmax=np.amax(np.abs(species)))

    plt.xlabel(r"$x$ [$c \omega_{pe}^{-1}$]")
    plt.ylabel(r"$p$ [$m_e c$]")
    plt.tight_layout()
    camera.snap()

animation = camera.animate(interval=50)
animation_name = './img/'+name+'.mp4'
animation.save(animation_name)