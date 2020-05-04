import numpy as np
import matplotlib.pyplot as plt

# number of nodes
n=100
# prepare the T and T+ arrays
T_plus=np.zeros(n)
T=np.zeros(n)
# time
sim_time=600
dt=1

#geometry and thermal props
L=0.1 # m
dx=L/n #
alpha=1e-7 #m2/s
Fo=alpha*dt/dx**2 # Fourier

# stability check
if Fo>0.5:
	print("stability issue")
# time
t=0
while t < sim_time:
	# boundary conditions
	T[0]=0
	T[n-1]=10
	# inside the domain 
	for i in range(1,n-1):
		T_plus[i]=T[i]*(1-2*Fo)+Fo*(T[i+1]+T[i-1])
	T=T_plus # replace T by T+
	t+=dt
# plotting
x_pos=np.arange(0,L,dx)
plt.xlabel("x position [m]")
plt.ylabel("Temperature [°C]")
plt.plot(x_pos, T_plus, 'o-', alpha=0.65)
plt.show()
# plt.savefig("./FDM_t"+str(sim_time)+".pdf") # optionally save
