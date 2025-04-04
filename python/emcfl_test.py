# import required libraries
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as scp

# Grid
Lx = 110
Nx = 256
dx = Lx/Nx
Nt = 10000
dt = 0.2
print("CFL = "+str(dt/dx))

# Init
B = np.zeros((Nx,Nt+1),dtype=np.complex128)
E = np.zeros((Nx,Nt+1),dtype=np.complex128)
Bf = np.zeros((Nx,6),dtype=np.complex128)
Ef = np.zeros((Nx,6),dtype=np.complex128)
qs = np.linspace(-np.pi/dx,np.pi/dx,Nx,endpoint=False)

method = 0
# 0: Yee
# 1: Spectral
# 2: 1st order space, 2nd order time FV
# 3: 1st order space, 1st order time FV
# 4: 2nd order space, 2nd order time FV
time = 400

# Initial Conditions
E[:,0] = np.sin(0.285599*(np.linspace(-55,55,Nx,endpoint=False)))
#E[Nx//2,0] = 0.001

print(E[10,0]**2)

if(method==1):
    for i in range(Nt):

        # Fourier Transform into Ef/Bf [0]
        Bf[:,0] = scp.fft.fftshift(scp.fft.fft(B[:,i]))
        Ef[:,0] = scp.fft.fftshift(scp.fft.fft(E[:,i]))
        # Create derivatives on #f/Bf [1]
        Bf[:,1] = Bf[:,0] - 1j * dt * qs[:] * Ef[:,0]
        Ef[:,1] = Ef[:,0] - 1j * dt * qs[:] * Bf[:,0]
        Bf[:,2] = Bf[:,1] - 1j * dt * qs[:] * Ef[:,1]
        Ef[:,2] = Ef[:,1] - 1j * dt * qs[:] * Bf[:,1]
        Bf[:,3] = (Bf[:,0] + Bf[:,2])/2
        Ef[:,3] = (Ef[:,0] + Ef[:,2])/2
        B[:,i+1] = scp.fft.ifft(scp.fft.ifftshift(Bf[:,3]))
        E[:,i+1] = scp.fft.ifft(scp.fft.ifftshift(Ef[:,3]))
if(method==0):
    for i in range(Nt):

        Bf[:,0] = B[:,i]
        Ef[:,0] = E[:,i]

        Bf[:,1] = Bf[:,0] - dt/(2*dx) * (np.roll(Ef[:,0],-1)-np.roll(Ef[:,0],1))
        Ef[:,1] = Ef[:,0] - dt/(2*dx) * (np.roll(Bf[:,0],-1)-np.roll(Bf[:,0],1))

        Bf[:,2] = Bf[:,1] - dt/(2*dx) * (np.roll(Ef[:,1],-1)-np.roll(Ef[:,1],1))
        Ef[:,2] = Ef[:,1] - dt/(2*dx) * (np.roll(Bf[:,1],-1)-np.roll(Bf[:,1],1))

        B[:,i+1] = (Bf[:,0]+Bf[:,2])/2
        E[:,i+1] = (Ef[:,0]+Ef[:,2])/2
if(method==2):
    for i in range(Nt):

        Bf[:,0] = B[:,i]
        Ef[:,0] = E[:,i]

        Bf[:,1] = Bf[:,0] - dt/(2*dx) * (np.roll(Ef[:,0],-1)-np.roll(Ef[:,0],1)-np.roll(Bf[:,0],-1)-np.roll(Bf[:,0],1)+2*Bf[:,0])
        Ef[:,1] = Ef[:,0] - dt/(2*dx) * (np.roll(Bf[:,0],-1)-np.roll(Bf[:,0],1)-np.roll(Ef[:,0],-1)-np.roll(Ef[:,0],1)+2*Ef[:,0])

        Bf[:,2] = Bf[:,1] - dt/(2*dx) * (np.roll(Ef[:,1],-1)-np.roll(Ef[:,1],1)-np.roll(Bf[:,1],-1)-np.roll(Bf[:,1],1)+2*Bf[:,1])
        Ef[:,2] = Ef[:,1] - dt/(2*dx) * (np.roll(Bf[:,1],-1)-np.roll(Bf[:,1],1)-np.roll(Ef[:,1],-1)-np.roll(Ef[:,1],1)+2*Ef[:,1])

        B[:,i+1] = (Bf[:,0]+Bf[:,2])/2
        E[:,i+1] = (Ef[:,0]+Ef[:,2])/2
if(method==3):
    for i in range(Nt):

        B[:,i+1] = B[:,i] - dt/(2*dx) * (np.roll(E[:,i],-1)-np.roll(E[:,i],1)-np.roll(B[:,i],-1)-np.roll(B[:,i],1)+2*B[:,i])
        E[:,i+1] = E[:,i] - dt/(2*dx) * (np.roll(B[:,i],-1)-np.roll(B[:,i],1)-np.roll(E[:,i],-1)-np.roll(E[:,i],1)+2*E[:,i])

def phi(x):
    return np.fmax(0,np.fmin((1+x)/2,np.fmin(2,2*x)))

if(method==4):
    for i in range(Nt):
        # slopes
        Bf[:,0] = phi((B[:,i]-np.roll(B[:,i],1))/(np.roll(B[:,i],-1)-B[:,i]))*(np.roll(B[:,i],-1)-B[:,i])
        Ef[:,0] = phi((E[:,i]-np.roll(E[:,i],1))/(np.roll(E[:,i],-1)-E[:,i]))*(np.roll(E[:,i],-1)-E[:,i])
        # minus values
        Bf[:,1] = B[:,i] - Bf[:,0]/2
        Ef[:,1] = E[:,i] - Ef[:,0]/2
        # plus values
        Bf[:,2] = B[:,i] + Bf[:,0]/2
        Ef[:,2] = E[:,i] + Ef[:,0]/2
        # fluxes at j+1/2
        Bf[:,3] = 1/2*(Ef[:,2]+np.roll(Ef[:,1],-1)+Bf[:,2]-np.roll(Bf[:,1],-1))
        Ef[:,3] = 1/2*(Bf[:,2]+np.roll(Bf[:,1],-1)+Ef[:,2]-np.roll(Ef[:,1],-1))
        # next step
        Bf[:,4] = B[:,i] - dt/dx*(Bf[:,3]-np.roll(Bf[:,3],1))
        Ef[:,4] = E[:,i] - dt/dx*(Ef[:,3]-np.roll(Ef[:,3],1))

        # slopes
        Bf[:,0] = phi((np.roll(Bf[:,4]-Bf[:,4],1))/(np.roll(Bf[:,4],-1)-Bf[:,4]))*(np.roll(Bf[:,4],-1)-Bf[:,4])
        Ef[:,0] = phi((np.roll(Ef[:,4]-Ef[:,4],1))/(np.roll(Ef[:,4],-1)-Ef[:,4]))*(np.roll(Ef[:,4],-1)-Ef[:,4])
        # minus values
        Bf[:,1] = Bf[:,4] - Bf[:,0]/2
        Ef[:,1] = Ef[:,4] - Ef[:,0]/2
        # plus values
        Bf[:,2] = Bf[:,4] + Bf[:,0]/2
        Ef[:,2] = Ef[:,4] + Ef[:,0]/2
        # fluxes at j+1/2
        Bf[:,3] = 1/2*(Ef[:,2]+np.roll(Ef[:,1],-1)+Bf[:,2]-np.roll(Bf[:,1],-1))
        Ef[:,3] = 1/2*(Bf[:,2]+np.roll(Bf[:,1],-1)+Ef[:,2]-np.roll(Ef[:,1],-1))
        # next step
        Bf[:,4] = Bf[:,4] - dt/dx*(Bf[:,3]-np.roll(Bf[:,3],1))
        Ef[:,4] = Ef[:,4] - dt/dx*(Ef[:,3]-np.roll(Ef[:,3],1))

        # average
        B[:,i+1] = (B[:,i]+Bf[:,4])/2
        E[:,i+1] = (E[:,i]+Ef[:,4])/2


dgNt= int(time/dt)

test=np.array([1,2,3,4,5,6,7,8,9])
print(np.roll(test,-1)-np.roll(test,1))

print(E[10,dgNt]**2+B[10,dgNt]**2)

plt.plot(np.linspace(-55,55,Nx),E[:,dgNt],c='r')
plt.plot(np.linspace(-55,55,Nx),B[:,dgNt],c='b')
plt.title("t="+str(dgNt*dt))
plt.show()

    

