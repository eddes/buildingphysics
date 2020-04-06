# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

#quelques couleurs
rouge_A='#C60C2E'
vert1_A='#005157'
vert2_A='#627D77'
vert3_A='#9EB28F'
vert4_A='#C5E5A4'
gris1_A='#595A5C'
coule=[rouge_A,vert1_A,vert2_A,vert3_A,vert4_A,gris1_A]

import numpy as np
import matplotlib.pyplot as plt
# distribution 'Dirac' function
def fc_distrib(T,dTf,Tf):
	return np.exp(- (T*(T-Tf)**2/dTf)/(np.sqrt(np.pi*dTf**2) ))
# function for the liquid fraction
def fc_fraction(T,dTf,Tf):
	if T<(Tf-dTf): f=0 # all solid
	elif T>(Tf+dTf): f=1 # all liquid
	else:f=(T-Tf+dTf)/(2*dTf)
	return f
# function for the apparent Cp 
def fc_Cp_apparent(T,dTf,Tf,Lf,Cp_s,Cp_l):
	fraction=fc_fraction(T,dTf,Tf)
	distrib=fc_distrib(T,dTf,Tf)
	return Lf*distrib + Cp_s + fraction*(Cp_l-Cp_s)

# Cps in J/kg/K
Cp_s=800
Cp_l=2000
Lf=188000 # Lf
dTf=0.01
Tmin=20
Tmax=30
T_init=(Tmin+Tmax)/2
Tf=27 # fusion temperature
T_init=Tf # initialize 
# material props
L=0.1# m
k=0.9 # conductivity (supposed equal for solid & liquid)
rho=800 # density (same)
# domain properties and initial conditions
n=10+2 # add 2 lines to the solid domain for the boundary conditions
K=np.eye(n,n,k=-1)*1 + np.eye(n,n)*-2 + np.eye(n,n,k=1)*1
K[0,0]=0
K[0,1]=0
K[-1,-1]=0
K[-1,-2]=0
dx=L/(n-2) # compute the actual dx (n minus two boundary condition nodes)
T_plus,Cp_t=np.zeros(n),np.zeros(n)
T=np.ones(n)*T_init # initialize
# simulation time and time step
t=0 
hours=0.1
sim_time=hours*3600 #s
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
plt.plot(x_pos, T_plus, color=coule[0], alpha=0.65, linestyle="--",marker='')
