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

def fc_IAQ_coupled_classes(vec_CLp, vec_CL, tau,delta,dt,rho,S,V,Ce):
	n=int(len(vec_CL)/2) # number of size classes
#	print(vec_CL)
	C,L = vec_CL[0:n], vec_CL[n:] # split into C and L
	Cp,Lp = vec_CLp[0:n], vec_CLp[n:]
	# solve with delta being a vector this time
	C_term = -Cp + C + 0.5*dt*( tau*(Ce-C) -delta*C+ rho*S*L/V) + 0.5*dt*( tau*(Ce-Cp) -delta*Cp+ rho*S*Lp/V)
	L_term = -Lp + L + 0.5*dt*( delta*V*C/S - rho*S*L/S ) + 0.5*dt*( delta*V*Cp/S - rho*S*Lp/S )
	return np.hstack([C_term,L_term])

# concentrations
Cinit=20
Linit=100
Ce_base=20
# enclosure properties
# enclosure
L,l,h=5,5,3# dimensions
S=2*( L*l+l*h+ h*L) # m2
V=L*l*h # m3 
tau=0.0# vol/h
qv=V*tau#m3/h
tau=qv/V
# particle size classes
d=[0.001,0.01,0.1,1,10]
# coeffs for distribution behaviour
rho_base=0 # resuspension rate
# deposition coeffs
delta=[2,0.01,0.005,0.01,4]
delta=[2,0.01,0.005,0.01,4]

dist_C_frac=[0.05,0.15,0.05,0.15,0.6]

# pour deux classes de taileles
#d=[0.01,10]
#dist_C_frac=[0.5,0.5]
#delta=[0.005,2]

delta=np.asarray(delta) # for array term by term mult
dist_C_frac=np.asarray(dist_C_frac)

n_classes=len(d)
rho = rho_base*np.ones(n_classes)
C = Cinit*np.ones(n_classes)*dist_C_frac
C0 = C # for plotting
L = Linit*np.ones(n_classes)*dist_C_frac
Ce = Ce_base*np.ones(n_classes)*dist_C_frac
concentration,deposition,time=[],[],[]
rhoC,deltaC=[],[]
rhop=[]

period=12
nb_period=1
dt=0.01 #h
sim_time=nb_period*24 # hour

sim_time=1 #
t=0 # hour
matrice_C,matrice_L,temps=[],[],[]
while t < sim_time:
	# store for plotting
	matrice_C.append(C)
	matrice_L.append(L)
	temps.append(round(t,2))
	# store for plotting
	concentration.append(C)
	deposition.append(L)
	time.append(round(t,2))
	rhop.append(rho)
	rhoC.append(rho*S*L/V)
	deltaC.append(delta*C)
	
	# variable Ce and rho
	Ce=dist_C_frac*(Ce_base+ 5*abs(np.cos(t*2*np.pi/period)))
	rho=rho_base#*abs(np.sin(t*2*np.pi/period))
	# constant Ce
	# np.hstack is used to flatten vector
	vec_CL =fsolve(fc_IAQ_coupled_classes, np.hstack([C,L]), args=(np.hstack([C,L]),tau,delta,dt,rho,S,V,Ce))
	C_plus,L_plus = vec_CL[0:n_classes], vec_CL[n_classes:] # split into C and L
	C,L=C_plus,L_plus
	t+=dt

plotit=True
if plotit==True:
	# mass fraction
	plt.clf()	
	fig = plt.figure()
	ax = fig.add_subplot(121)
	plt.bar(np.arange(0,n_classes,1),dist_C_frac, color=coule[-1])
	plt.xticks(np.arange(0,n_classes,1), d)
	plt.xlabel("Particle size class [µm]")
	plt.ylabel("Mass fraction repartition [-]")
	plt.tight_layout()
	# mass fraction
	ax = fig.add_subplot(122)
	
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.plot(d, delta, color=coule[-1],marker='^',ls='--')
	plt.xlabel("Particle size class [µm]")
	plt.ylabel(r"Deposition coefficient $\delta$ [h$^{-1}$]")
	plt.tight_layout()
	
	
	# 
	plt.clf()
	plt.subplot(121)
	plt.xlabel("Time [h]")
	plt.ylabel(r"Concentration [µg/m$^3$]")
	plt.plot(temps, C_tot, '-',color=coule[0], alpha=0.65,label='C')
	plt.legend()
	
	plt.subplot(122)
	plt.xlabel("Time [h]")
	plt.ylabel(r"Mass on surfaces [µg/m$^2$]")
	plt.plot(temps, L_tot, '-',color=coule[-1], alpha=0.65,label='L')
	plt.legend()
	plt.tight_layout()
