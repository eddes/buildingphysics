import numpy as np
import matplotlib.pyplot as plt

n=50+2 # number of nodes
K=np.eye(n,n,k=-1)*1 + np.eye(n,n)*-2 + np.eye(n,n,k=1)*1
# modify the terms related to boundary conditions
K[0,0],K[0,1]=0,0
K[-1,-1],K[-1,-2]=0,0

# time
sim_time=600 #s
dt=1 #s
# geometry and thermal properties
L=0.1 # m
dx=L/n # m
alpha=1e-7 #m2/s
Fo=alpha*dt/dx**2 #compute Fourier's number
# initial conditions
T=np.zeros(n)
# boundary conditions
T[0]=0
T[n-1]=10

# stability check for Euler explicit in 1D
if Fo>0.5:
	print("stability issue")
# time
t=0
while t < sim_time:
	# matrix multiplication K*T
	T_plus=Fo*np.dot(K,T) +T
	T=T_plus # replace 
	t+=dt

# plotting
x_pos=np.arange(0,L,dx)
plt.xlabel("x position [m]")
plt.ylabel("Temperature [°C]")
plt.plot(x_pos, T_plus, 'o-', alpha=0.65)
plt.show()
