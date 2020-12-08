# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import matplotlib.cm as cm
from tqdm import tqdm

staggered=True# pipe alignment
color = cm.RdYlBu # RdYlBu_r
#############################################
#
#	Domain setup : geometry, properties...
#
n=100
m=n
L=40 # m
dx=L/n#
depth=15 # m of pipes

# time
t=0
days=200
sim_time=days*24*3600
 
# soil properties
rho=2150
Cp=2000
alpha=1.8/(rho*Cp) #m2/s
# water properties
rho_w=1000
Cp_w=4182
v_w=0.4 # underground water velocity in m/day
v_w= v_w/(24*3600) # m/day > m/s
# ratio of heat content per degree [J/K]
r=(rho_w*Cp_w)/(rho*Cp)
# for the COP computation
eta_carnot=0.35 # efficiency versus the theoretical Carnot COP (own fit against manufacturer data)
# temperatures
T_init=10 # initialisation
T_inlet=10 # underground water temperature entering the domain
T=np.ones((n,m))*T_init
T_temp=np.zeros((n,m)) # 

# computation of the stable dt for the explicit scheme
dt_max=1./( (2*r*v_w)/dx  + (4*alpha)/dx**2 )
dt=int(0.95*dt_max) # stay on the safe side
print("stable dt = ", dt, ' s')
# take it easy with stability/precision compromise > do not exceed 1h
if dt>3600:
	dt=3600
	print("chosen dt = ", dt, ' s')
# Courant and Fourier adim numbers
Co=v_w*dt/dx
Fo=alpha*dt/dx**2

# interpolate outdoor conditions at the chosen time step
Tout=np.load('./Ta_out.npy')
t_epw=np.arange(0,len(Tout))*3600
t_new= np.arange(0,t_epw[-1], dt)
x,y = t_epw,Tout
y_interp = scipy.interpolate.interp1d(x, y)
Tout=y_interp(t_new)

# prepare the computation of P_l the power taken to the ground
Pmax=10000 # max capacity for heating [W]
dTmax=(19--6) # sizing temperature delta

#############################################
#
#	Setup of the geothermal pipes
#
n_pipes=9
n_rows=3 # use multiples of n_pipes,  otherwise a pipe may miss out...
n_ppr=int(n_pipes/n_rows) # pipes per row
ij_pipes=[] # list of coordinates

# determination of the pipe position
for j in range(1,n_rows+1):
	for k in range(1,n_ppr+1):
		i_pipes = int( j*(n/4)/(n_rows+1))
		if staggered==True and j%2==0:
			j_pipes= int(m/4)+int(k*(m/2)/((n_ppr+1)))+ int(m/2/((n_ppr+1))/2)
		else:
			j_pipes = int(m/4)+int(k*(m/2)/((n_ppr+1)))
		ij_pipes.append([i_pipes, j_pipes])

T_avg=T_init # average pipe temp for the 1st COP computation

ii=0 # increment for the outdoor temperature
for t in tqdm(np.arange(0,sim_time,dt)):
	# compute the heating departure temperature
	T_hot=-1.05*Tout[ii]+33.7
	# heating power required
	P_h   = (19-Tout[ii])/dTmax*Pmax 
	# compute the heating mode COP
	COP_h = eta_carnot*(T_hot+273.15+5)/(T_hot+5-T_avg-5)
	# compute the cooling mode COP
	COP_c = COP_h - 1
	# compute the electrica
	Pcomp = P_h/COP_h
	# compute the cooling power in the ground
	P_c = Pcomp*COP_c
	# update the heat pump power (in W/m)
	P_l=-P_c/n_pipes/depth
	source_HP=P_l*dt/(rho*Cp*dx**2)
	# in/out boundary conditions
	for j in range(0,m):
		# inlet
		T[0,j]=T_inlet
		T_temp[0,j]=T[0,j]
		# outlet
		T[n-1,j]=T[n-2,j]
		T_temp[n-1,j]=T[n-1,j]
	# adiabatic boundary conditions
	for i in range(1,n-1):
		# j=0 
		T[i,0]= T[i,0]*(1-r*Co-3*Fo)\
				+ Fo * (T[i+1,0] + T[i-1,0]) \
				+ Fo * (T[i,1])\
				+ r*Co*(T[i-1,0])
		T_temp[i,0]=T[i,0]
		# j=m-1
		T[i,m-1]= T[i,m-1]*(1-r*Co-3*Fo)\
				+ Fo * (T[i-1,m-1] + T[i+1,m-1] ) \
				+ Fo * (T[i,m-2]) \
				+ r*Co*(T[i-1,m-2])
		T_temp[i,m-1]=T[i,m-1]
	# inside the domain
	for i in range(1,n-1):
		for j in range(1,m-1):
			if [i,j] in ij_pipes:
				source=source_HP
			else:
				source=0
			T_temp[i,j]=T[i,j]*(1-r*Co-4*Fo)\
					+ Fo * (T[i+1,j]+T[i-1,j]) \
					+ Fo * (T[i,j+1]+T[i,j-1]) \
					+ r*Co*(T[i-1,j]) \
					+ source
	T=T_temp
	T_plus=T_temp
	t+=dt
	# compute the average pipe temperature
	T_avg=0
	for i,elt in enumerate(ij_pipes):
		T_avg+=T[elt[0],elt[1]]
		
	T_avg=T_avg/len(ij_pipes)
	ii+=1
###################################
#
#	plotting/post-proc
#
plt.clf()
fig = plt.figure()
ax0 = fig.add_subplot()
surf = ax0.imshow(T_plus,cmap=color,interpolation='bicubic',extent=[0,L,0,L])
plt.colorbar(surf,label='Temperature [°C]')
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.tight_layout()
plt.savefig('./geothermal_2D.pdf')
