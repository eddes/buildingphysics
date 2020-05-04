import numpy as np
import matplotlib.pyplot as plt

n=50+2 # number of nodes
K=np.eye(n,n,k=-1)*1 + np.eye(n,n)*-2 + np.eye(n,n,k=1)*1

K[0,0]=0
K[0,1]=0
K[-1,-1]=0
K[-1,-2]=0


T_plus=np.zeros(n)
T=np.zeros(n)


t=0
sim_time=600
dt=1

L=0.1 # m
alpha=1e-7 #m2/s
dx=L/n # 

Fo=alpha*dt/dx**2
T[0]=0
T[n-1]=10
# stability check
if Fo>0.5:
	print("stability issue")
# time 
while t < sim_time:
	# matrix multiplicatio K*T
	T_plus=Fo*np.dot(K,T) +T
	T=T_plus # replace 
	t+=dt

x_pos=np.arange(0,L,dx)
plt.xlabel("x position [m]")
plt.ylabel("Temperature [°C]")
plt.plot(x_pos, T_plus, 'o-', alpha=0.65)
plt.savefig("./EEM_t"+str(sim_time)+".png",dpi=200)
