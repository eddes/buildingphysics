import numpy as np
import matplotlib.pyplot as plt

def fc_distrib(T,dTf,Tf):
	return np.exp(- (T*(T-Tf)**2/dTf)/(np.sqrt(np.pi*dTf**2) ))
def fc_fraction(T,dTf,Tf):
	if T<(Tf-dTf): f=0 # all solid
	elif T>(Tf+dTf): f=1 # all liquid
	else:f=(T-Tf+dTf)/(2*dTf)
	return f
def fc_Cp_apparent(T,dTf,Tf,Lf,Cp_s,Cp_l):
	fraction=fc_fraction(T,dTf,Tf)
	distrib=fc_distrib(T,dTf,Tf)
	return Lf*distrib + Cp_s + fraction*(Cp_l-Cp_s)

# Cps in J/kg/K
Cp_s=1800
Cp_l=2400
Lf=188000 # Lf
dTf=0.01
Tmin=20
Tmax=30
Tf=27
T_init=Tf

# material props
L=0.1# m
k=0.9
rho=800

# domain properties and initial conditions
n=10+2
K=np.eye(n,n,k=-1)*1 + np.eye(n,n)*-2 + np.eye(n,n,k=1)*1
K[0,0]=0
K[0,1]=0
K[-1,-1]=0
K[-1,-2]=0
dx=L/n #
T_plus,Cp_t=np.zeros(n),np.zeros(n)
T=np.ones(n)*T_init # initialize at Tmin

# simulation time and time step
t=0 
heures=0.1
sim_time=heures*3600 #s
dt=5 # s

 #intialize the local Fourier number for PCMs
Fo=np.zeros(n)
for i in range(len(Fo)):
	Fo[i]=k*dt/(rho*fc_Cp_apparent(T_init,dTf,Tf,Lf,Cp_s,Cp_l)*dx**2)
	# local stability check
	if Fo[i]>0.5:
		print("stability issue... i=",i)
		dt_min=0.5*dx**2*rho*Cp_s/k
		print("minimum time step =", round(dt_min,2))
		dt=0.9*dt_min
		print("changing to =", round(dt,2))
# time 
while t < sim_time:
	# boundary conditions
	T[0]=Tmin
#	T[n-1]=Tmin+(Tmax-Tmin)*t/sim_time
	T[n-1]=Tmax
	# inside the domain 
	T_plus=Fo*np.dot(K,T) +T
	# update the local apparent Cp
	for i in range(n):
		Cp_t[i] = fc_Cp_apparent(T_plus[i],dTf,Tf,Lf,Cp_s,Cp_l)
	Fo=k*dt/(rho*Cp_t*dx**2)
	if len(np.argwhere(Fo > 0.5))>0 : 
		indice=np.argwhere(Fo > 0.5)
		print("stability issue during simulation... Fo =", Fo[indice[0]])
		print("occurred for... t=", round(t,1), " / i=")
	t+=dt #time increment
	T=T_plus # T turns into T_plus to allow the calculation of the next T_plus

x_pos=np.arange(0,L,dx)
plt.plot(x_pos, np.ones(n)*Tf,color=coule[-1],linestyle="-",alpha=0.9,marker='')
plt.xlabel("x position [m]")
plt.ylabel("Temperature [°C]")
plt.plot(x_pos, T_plus, alpha=0.65, linestyle="--",marker='')
