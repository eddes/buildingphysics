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
# vapour diffusivity depending on water content 
def fc_deltav(w):
	a=1.1*1e-7
	b=-1.57*1e-10
	return a+b*w
# thermal conductivity depending on water content 
def fc_lambda(w):
	a=0.23
	b=0.006
	return a+b*w
# sorption curve
def fc_w_phi(phi):
	a=700 #1000
	b=250*1e-6 # 146*1e-6
	c=2.9#1.59
	# let's try to avoid numerical problems 
	# with invalid log values...
	try: 
		res=a/np.power( 1- np.log(phi/100)/b, 1/c)
	# ... hence the upper limit for w
	except: 
		print(phi)
		res=600
	return res

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
def fc_coupled_HAM(vec_Tpvp,vec_Tpv,K0,dt,rho,Cp,Lv):
	# split array 
	n=int(len(vec_Tpv)/2) # half index
	T,pv   = vec_Tpv[0:n] ,vec_Tpv[n:]
	Tp,pvp = vec_Tpvp[0:n],vec_Tpvp[n:]
	# we need phi in order to actualise w
	phi=pv/fc_pvsat(T)*100
	# compute w 
	w=fc_w_phi(phi)
	# update the properties 
	Cpm,dw,dphi=update_Cpm(w,phi,T)
	Fo=update_Fo(w,k,rho,Cp,dt,dx)
	Fow=update_Fow(w,dt,dx)
	K,Kw,KwT = fc_update_Kmatrices(K0,Fo,Fow,Cpm)
	# Crank-Nicolson scheme
	# explicit and implicit parts for T
	exp_T=0.5 * ( np.dot(K,T)  + Lv /(rho * Cp) * np.dot(KwT,pv))
	imp_T=0.5 * ( np.dot(K,Tp) + Lv /(rho * Cp) * np.dot(KwT,pvp))
	T_term = exp_T + imp_T
	# explicit and implicit parts for pv
	exp_pv=0.5 * np.dot(Kw,pv)
	imp_pv=0.5 * np.dot(Kw,pvp)
	pv_term = exp_pv + imp_pv
	# send back as one array
	return np.hstack([T_term,pv_term])

# domain size
n_solid=15
n=n_solid+2
L=0.1 # m
# heat transfer matrix
K0=np.eye(n,n,k=-1)*1 + np.eye(n,n)*-2 + np.eye(n,n,k=1)*1
K0[0,0],K0[0,1],K0[-1,-1],K0[-1,-2]=0,0,0,0

#####################
# 
# TIME  TIME TIME TIME TIME
#
t=0
dt=60
period=24 # pour la sinusoide
period_sec=period*3600
nb_period=0.5 # 365
sim_time=nb_period*24*3600 # seconds
modulo_storage=int(0.1*3600) #sim_time/dt/100

##################

# physical props
Lv=2400*1e3 # J/kg
k=1.6
rho=2800
Cp=1000
alpha=k/rho/Cp #1e-7 #m2/s
dx=L/(n_solid) #
Fo=alpha*dt/dx**2

# sur,face transfer coefficients
ha,hb=8,25 # W/m2/K
hv = 3*1e-8 # m2/s

# boundary conditions
Tleft=10
pvleft=0.8*fc_pvsat(Tleft)
phileft=pvleft/(fc_pvsat(Tleft))*100

Tright=Tleft
pvright=pvleft
phiright=pvright/(fc_pvsat(Tright))*100

## vectors
pv=np.ones(n)
dw=np.zeros(n)
w=np.ones(n)
T=np.ones(n)

# initial conditions
pv=np.ones(n)*pvleft*0.95
T=np.ones(n)*Tleft*0.95

def fc_update_Kmatrices(K0,Fo,Fow,Cpm):
	# HEAT - coefficients for K
	K =Fo*K0
	K[0,0],K[0,1],K[-1,-1],K[-1,-2]=0,0,0,0
	# exchange with air
	Fo_eqa,Fo_eqb=dt/(dx*(1/ha+dx/(2*fc_lambda(w[1])))), dt/(dx*(1/hb+dx/(2*fc_lambda(w[1]))))
	K[1,0],	K[1,1]= Fo_eqa, -Fo_eqa-Fo[1]
	# exchange with air
	K[-2,-1], K[-2,-2]= Fo_eqb,-Fo_eqb-Fo[-2]
  
	# MASS - coefficients for Kw
	Kw =Fow/Cpm*K0
	Kw[0,0],Kw[0,1],Kw[-1,-1],Kw[-1,-2]=0,0,0,0
	# exchange with air
	Fow_eqa, Fow_eqb = dt/(dx*(1/hv+dx/(2*fc_deltav(w[1])))), dt/(dx*(1/hv+dx/(2*fc_deltav(w[-2]))))
	Kw[1,0], Kw[-2,-1] = Fow_eqa/Cpm[1], Fow_eqb/Cpm[-2]
	Kw[1,1], Kw[-2,-2] =(-Fow_eqa-Fow[1])/Cpm[1],(-Fow_eqb-Fow[-2])/Cpm[-2]
	
	# COUPLING matrix - coefficients for KwT
	KwT =Fow/Cpm*K0
	KwT[0,0],KwT[0,1],KwT[-1,-1],KwT[-1,-2]=0,0,0,0
	# exchange with air
	KwT[1,0], KwT[-2,-1] =  Fow_eqa, Fow_eqb
	KwT[1,1], KwT[-2,-2]= -Fow_eqa-Fow[1], -Fow_eqb-Fow[-2]
	
  return K,Kw,KwT # three K's for the dot product

# coefficients for heat and mass transfer
phi=pv/fc_pvsat(T)*100
w=fc_w_phi(phi)

#-------------------------------------------------
# 	Preparing the boundary conditions
#
# load pressure and temperature
pvleft_epw=np.load("pv_out.npy")
Tleft_epw=np.load("Ta_out.npy")

pvleft_epw=pvleft_epw[0:int(24*(nb_period+1))]
Tleft_epw=Tleft_epw[0:int(24*(nb_period+1))]

#number of seconds
t_epw=np.arange(0,len(pvleft_epw))*3600
# resample pv
x,y = t_epw,pvleft_epw
y_interp = scipy.interpolate.interp1d(x, y)
t_new= np.arange(0,t_epw[-1], dt)
pvleft=y_interp(t_new)
#resample temperature
x,y = t_epw,Tleft_epw
y_interp = scipy.interpolate.interp1d(x, y)
Tleft=y_interp(t_new)


i=0
# time loop
pbar = tqdm(total=sim_time)
while t <= sim_time:
	# update boundary conditions
	T[0]  = Tleft[i]
	T[n-1]= Tright # + 1.5*np.sin(t*2*np.pi/period_sec)
	pv[0] = pvleft[i]
	pv[n-1]=pvright #+ 100*np.cos(t*2*np.pi/period_sec)

	# solve the coupled, non-linear system
	result_array=fsolve(fc_coupled_HAM,
			 np.hstack([T,pv]),
			 args=(np.hstack([T,pv]), K0, dt, rho, Cp, Lv))
	
	# split the result into T and pv
	T_plus,pv_plus=result_array[0:n], result_array[n:]
	
	# update variables
	pv = pv_plus
	T=T_plus
	t+=dt
	i+=1
	pbar.update(dt)
pbar.close()

elapsed_time= (time.time() - start_time) # in seconds
print('computational effort ',round(elapsed_time/(sim_time/3600),2), ' [seconds/hour simulation]')
print('elapsed time :', round(elapsed_time/3600), 'h', round((elapsed_time % 3600)/60), 'min')
