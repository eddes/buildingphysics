import numpy as np
import matplotlib.pyplot as plt

n=100
T_plus=np.zeros(n)
T=np.zeros(n)

t=0
sim_time=600
dt=1

L=0.1 # m
alpha=1e-7 #m2/s
dx=L/n # 

Fo=alpha*dt/dx**2
#Fo=np.ones(n)*0.25 # >> local Fo for PCMs
# 
if Fo>0.5:
	print("stability issue")
# time 
while t < sim_time:
	# boundary conditions
	T[0]=0
	T[n-1]=10
	# inside the domain 
	for i in range(1,n-1):
		T_plus[i]=T[i]*(1-2*Fo)+Fo*(T[i+1]+T[i-1])
	T=T_plus # replace 
	t+=dt
	
x_pos=np.arange(0,L,dx)
plt.xlabel("x position [m]")
plt.ylabel("Temperature [°C]")
plt.plot(x_pos, T_plus, 'o-', alpha=0.65)
plt.savefig("./FDM_t"+str(sim_time)+".png",dpi=200)
