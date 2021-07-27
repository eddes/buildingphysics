# -*- coding: utf-8 -*-
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
from  scipy.optimize import fsolve
from tqdm import tqdm

import time
start_time = time.time()
############################################################
#
#	Thermophysical functions
#
# vapour pressure at saturation
def fc_pvsat(T):
	return np.power(10, 7.625*T/(241+T)+2.7877)
# vapour diffusivity depending on water content 
def fc_deltav(w):
	a=1.1*1e-7
	b=-1.57*1e-7
	return a+b*w
# thermal conductivity depending on water content 
def fc_lambda(w):
	a=0.23
	b=6
	return a+b*w
# sorption curve
def fc_w_phi(phi):
	a=700 #1000
	b=250*1e-6 # 146*1e-6
	c=2.9#1.59
	res=a/np.power( 1- np.log(phi/100)/b, 1/c)
	return res*1e-3
# sortpion curve the other way around
def fc_phi_w(w):
	a=700 #1000
	b=250*1e-6 # 146*1e-6
	c=2.9#1.59
	phi =np.zeros(len(w))
	phi=np.exp(b*(1-np.power((a/(1000*w)),c)))*100
	return phi

# update local Fourier 
def update_Fo(w_vec,k,rho,Cp,dt,dx):
	k=fc_lambda(w) # variable properties
	Cp_vec=Cp+w_vec/rho
	Fo_vec=k/rho/Cp_vec*dt/dx**2
	return Fo_vec

# update local mass Fourier 
def update_Fow(w_vec,dt,dx):
	deltav=fc_deltav(w) # variable properties
	Fow=deltav*dt/dx**2
	return Fow

# update local Fourier 
def update_Cm(w,phi,T):
	epsilon=0.001*np.min(w)
	wp=w+epsilon
	phip=fc_phi_w(wp)
	dw=abs(wp-w)
	dphi=abs(phip-phi)/100
	pvs=fc_pvsat(T)
	Cm=dw/dphi/pvs
	return Cm

############################################################
#
#	Equation solving
# 
def fc_coupled_HAM(vec_Tpvp,vec_Tpv,K,Fo,Fow,dt,rho,Cp,Cm,Lv):
	# split array 
	n=int(len(vec_Tpv)/2) # half index
	T,pv   = vec_Tpv[0:n] ,vec_Tpv[n:]
	Tp,pvp = vec_Tpvp[0:n],vec_Tpvp[n:]
	# we need phi in order to actualise w
	phi=pv/fc_pvsat(T)*100
	# compute w 
	w=fc_w_phi(phi)
	# update the properties 
	Cm = update_Cm(w,phi,T)
	Fo  = update_Fo(w,k,rho,Cp,dt,dx)
	Fow = update_Fow(w,dt,dx)
	# Crank-Nicolson scheme
	# explicit and implicit parts for T
	exp_T=0.5 * (Fo * np.dot(K,T) +  Lv * Fow/(rho * Cp) * np.dot(K,pv))
	imp_T=0.5 * (Fo * np.dot(K,Tp) + Lv * Fow/(rho * Cp) * np.dot(K,pvp))
	T_term = -Tp + T + exp_T + imp_T
	# explicit and implicit parts for pv
	exp_pv=0.5 * (Fow/Cm * np.dot(K,pv))
	imp_pv=0.5 * (Fow/Cm * np.dot(K,pvp))
	pv_term = -pvp + pv + exp_pv + imp_pv
	# send back as one array
	return np.hstack([T_term,pv_term])

# domain size
n_solid=15
n=n_solid+2
L=0.1 # m
# heat transfer matrix
K=np.eye(n,n,k=-1)*1 + np.eye(n,n)*-2 + np.eye(n,n,k=1)*1
K[0,0],K[0,1],K[-1,-1],K[-1,-2]=0,0,0,0

#####################
# Time & storage
t=0
dt=60
period=24 #
period_sec=period*3600
nb_period=0.1#
sim_time=nb_period*24*3600 # seconds
modulo_storage=int(0.1*3600) #sim_time/dt/100
# preparing post-process
store_Text,store_pvext=[],[]
store_w,store_phi,store_T,store_pv=[],[],[],[]
# arrays
pv,dw,w,T=np.ones(n),np.ones(n),np.ones(n),np.ones(n)

# physical props
Lv=2400*1e3 # J/kg
k=1.6
rho=2800
Cp=1000
alpha=k/rho/Cp #1e-7 #m2/s
dx=L/(n_solid) #
Fo=alpha*dt/dx**2
# boundary conditions
Tleft=20
pvleft=1200
Tright=Tleft
pvright=pvleft
# initial conditions
pv_init=np.ones(n)*1300
T_init=np.ones(n)*15
pv=pv_init
T=T_init
phi=pv/fc_pvsat(T)*100
w=fc_w_phi(phi)
w_init=w
phi_init=phi

Cm = update_Cm(w,phi,T)
Fo  = update_Fo(w,k,rho,Cp,dt,dx)
Fow = update_Fow(w,dt,dx)
i=0
pbar=tqdm(total=sim_time) #set up a progress par
# time loop
while t <= sim_time:
	# update boundary conditions
	T[0]  = Tleft
	T[-1] = Tleft +20
	pv[0] = pvleft
	pv[-1]= pvleft +300

	# solve the coupled, non-linear system
	result_array=fsolve(fc_coupled_HAM,
			 np.hstack([T,pv]),
			 args=(np.hstack([T,pv]), K, Fo, Fow, dt, rho, Cp, Cm, Lv))
	
	# split the result into T and pv
	T_plus,pv_plus=result_array[0:n], result_array[n:]
	
	# compute phi
	phi=pv/fc_pvsat(T)*100
	# compute water content
	w=fc_w_phi(phi)

	# do some storage for plotting
	if (int(t) % modulo_storage)==0 and t!=0:
		store_w.append(w[1:-1]*1000)
		store_phi.append(phi[1:-1])
		store_T.append(T[1:-1])
		store_pv.append(pv[1:-1])
	# update variables
	pv = pv_plus
	T=T_plus
	t+=dt
	i+=1
	pbar.update(dt) # update progress bar
pbar.close() # close it

elapsed_time= (time.time() - start_time) # in seconds
print('computational effort ',round(elapsed_time/(sim_time/3600),2), ' [seconds/hour simulation]')
print('elapsed time :', round(elapsed_time/3600), 'h', round((elapsed_time % 3600)/60), 'min')
x_pos=np.arange(int(L/dx)+2)*dx
x_pos=x_pos-dx/2
x_pos=x_pos[1:-1]

print("\n#################\nPlotting (can be long)")
stop=int(len(store_phi))
start=0
############################################################
##############################
plt.subplot(121)
plt.xlabel("x position [m]")
plt.ylabel(r"$\varphi$ [%]")
plt.plot(x_pos,phi_init[1:-1], '--',color=coule[1], alpha=0.5)
for i in range(start,stop):
	plt.plot(x_pos, store_phi[i], '-',color=coule[1], alpha=0.35)

plt.subplot(122)
plt.xlabel("x position [m]")
plt.ylabel("water content [g/m$^3$]")
plt.plot(x_pos,w_init[1:-1]*1000, '--',color=coule[3], alpha=0.5)
for i in range(start,stop):
	plt.plot(x_pos, store_w[i], '-',color=coule[3], alpha=0.15)

plt.tight_layout()
plt.savefig("./graph_duo_phi_w.pdf")
