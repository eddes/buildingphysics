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
import scipy.interpolate
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

# Cp depending on water content
def fc_Cp_w(w_vec,rho,Cp):
	Cp_vec=np.ones(len(w_vec))
	Cpw=4182
	return Cp*Cp_vec+w/rho*Cpw

def fc_deltav(w):
	a=1.1*1e-7
	b=-1.57*1e-10
	return a+b*w

def fc_lambda(w):
	a=0.23
	b=0.006
	return a+b*w
def fc_w_phi(phi):
	a=700 #1000
	b=250*1e-6 # 146*1e-6
	c=2.9#1.59
	return a/np.power( 1- np.log(phi/100)/b, 1/c)

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
def update_Cpm(w,phi,T):
	epsilon=0.005*min(w)
	wp=w+epsilon
	phip=fc_w_phi(wp)
	dw=wp-w
	dphi=phip-phi
	pvs=fc_pvsat(T)
	Cpm=abs(dw/dphi)/pvs
	return Cpm,dw,dphi

############################################################
#
#	Equation solving
# 
def fc_coupled_HAM(vec_Tpvp,vec_Tpv,K,Fo,Fow,dt,rho,Cp,Cpm,Lv):
	rcp=rho*Cp
	# split array 
	n=int(len(vec_Tpv)/2) # number of size classes
	T,pv   = vec_Tpv[0:n] ,vec_Tpv[n:]
	Tp,pvp = vec_Tpvp[0:n],vec_Tpvp[n:]
	# Crank-Nicolson scheme
	T_term  = -Tp  +T   +0.5*(Fo*np.dot(K,T) + Lv/rcp*Fow*np.dot(K,pv) ) +0.5*(Fo*np.dot(K,Tp)+Lv/rcp*Fow*np.dot(K,pvp))
	pv_term = -pvp + pv +0.5*(Fow/Cpm*np.dot(K,pv) ) +0.5*(Fow/Cpm*np.dot(K,pvp))
	# unite array
	return np.hstack([T_term,pv_term])


def fc_coupled_HAM2(vec_Tpvp,vec_Tpv,K,Fo,Fow,dt,rho,Cp,Cpm,Lv):

	# split array 
	n=int(len(vec_Tpv)/2) # number of size classes
	T,pv   = vec_Tpv[0:n] ,vec_Tpv[n:]
	Tp,pvp = vec_Tpvp[0:n],vec_Tpvp[n:]
	
	# prepare updating the props
	phi=pv/fc_pvsat(T)*100
	phi[phi>100] = 99.5
	phi[phi<1] = 0.5
	# compute water content
	w=fc_w_phi(phi)
	# update the properties
	Cpm,dw,dphi=update_Cpm(w,phi,T)
	Fo=update_Fo(w,k,rho,Cp,dt,dx)
	Fow=update_Fow(w,dt,dx)
	
	rcp=rho*Cp
	# Crank-Nicolson scheme
	T_term  = -Tp  +T   +0.5*(Fo*np.dot(K,T) + Lv/rcp*Fow*np.dot(K,pv) ) +0.5*(Fo*np.dot(K,Tp)+Lv/rcp*Fow*np.dot(K,pvp))
	pv_term = -pvp + pv +0.5*(Fow/Cpm*np.dot(K,pv) ) +0.5*(Fow/Cpm*np.dot(K,pvp))
	# unite array
	return np.hstack([T_term,pv_term])

# domain size
n_solid=15
n=n_solid+2
L=0.1 # m
# heat transfer matrix
K=np.eye(n,n,k=-1)*1 + np.eye(n,n)*-2 + np.eye(n,n,k=1)*1
K[0,0],K[0,1],K[-1,-1],K[-1,-2]=0,0,0,0

#####################
# TIME  TIME TIME TIME TIME
dt=100
##################

# physical props
Lv=2400*1e3 # J/kg
k=1.6
rho=2800
Cp=1000
alpha=k/rho/Cp #1e-7 #m2/s
dx=L/(n_solid) #
Fo=alpha*dt/dx**2

# boundary conditions
Tleft=10
Tright=20
pvleft=0.8*fc_pvsat(Tleft)
pvright=0.5*fc_pvsat(Tright)
phileft=pvleft/(fc_pvsat(Tleft))*100
phiright=pvright/(fc_pvsat(Tright))*100

## vectors
pv=np.ones(n)
dw=np.zeros(n)
w=np.ones(n)
T=np.ones(n)

# initial conditions
pv=np.ones(n)*(pvleft+pvright)/2
T=np.ones(n)*(Tleft+Tright)/2

Cp_test=fc_Cp_w(w,rho,Cp)

phi=pv/fc_pvsat(T)*100
w=fc_w_phi(phi)
Cpm,dw,dphi=update_Cpm(w,phi,T)
Fo=update_Fo(w,k,rho,Cp,dt,dx)
Fow=update_Fow(w,dt,dx)

#-------------------------------------------------
# 	Preparing the boundary conditions
#
# load pressure and temperature
pvleft_epw=np.load("pv_out.npy")
Tleft_epw=np.load("Ta_out.npy")
# number of seconds
t_epw=np.arange(0,len(pvleft_epw))*3600
# resample pv
x,y = t_epw,pvleft_epw
y_interp = scipy.interpolate.interp1d(x, y)
t_new= np.arange(0,t_epw[-1], dt)
pvleft=y_interp(t_new)
# resample temperature
x,y = t_epw,Tleft_epw
y_interp = scipy.interpolate.interp1d(x, y)
Tleft=y_interp(t_new)

store_Text,store_pvext=[],[]
store_w,store_phi,store_T,store_pv=[],[],[],[]
t=0
period=24 # pour la sinusoide
period_sec=period*3600
nb_period=50
sim_time=nb_period*24*3600 # seconds
modulo_storage=3600 #sim_time/dt/100

i=0
# time loop
pbar = tqdm(total=sim_time)
while t < sim_time:
	# update boundary conditions
	T[0]  = Tleft[i]
	T[n-1]= Tright + 1.5*np.sin(t*2*np.pi/period_sec)
	pv[0] = pvleft[i]
	pv[n-1]=pvright + 100*np.sin(t*2*np.pi/period_sec)

	# solve the coupled, non-linear system
	result_array=fsolve(fc_coupled_HAM2,
			 np.hstack([T,pv]),
			 args=(np.hstack([T,pv]), K, Fo, Fow, dt, rho, Cp, Cpm, Lv))
	
	# split the result into T and pv
	T_plus,pv_plus=result_array[0:n], result_array[n:]
	
	# compute phi
	phi=pv/fc_pvsat(T)*100
	phi[phi>100] = 99.5
	phi[phi<1] = 0.5
	# compute water content
	w=fc_w_phi(phi)
	# update the properties
	Cpm,dw,dphi=update_Cpm(w,phi,T)
	Fo=update_Fo(w,k,rho,Cp,dt,dx)
	Fow=update_Fow(w,dt,dx)
  
	# update variables
	pv = pv_plus
	T=T_plus
	t+=dt
	i+=1
	pbar.update(dt)
pbar.close()

elapsed_time= (time.time() - start_time) # in seconds
print('computational effort ',round(elapsed_time/(sim_time/3600),2), ' [seconds/hour simulation]')
