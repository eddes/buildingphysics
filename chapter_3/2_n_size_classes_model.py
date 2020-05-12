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

# Crank-Nicolson scheme for polydisperse aerosol
def fc_IAQ_coupled_classes(vec_CLp, vec_CL, tau,delta,dt,rho,S,V,Ce):
	n=int(len(vec_CL)/2) # number of size classes
	C,L = vec_CL[0:n], vec_CL[n:] # array splitting into C and L
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
tau=0.1# vol/h
qv=V*tau#m3/h

# particle size classes
d=[0.001,0.01,0.1,1,10]
#initial mass distribution among the size classes (the sum must equal one)
dist_C_frac=[0.05,0.15,0.05,0.15,0.6]

# coeffs for distribution behaviour
rho_base=0.1 # resuspension rate
delta=[2,0.01,0.005,0.01,4] # deposition coeffs

# convert to numpy array for term by term multiplication
delta=np.asarray(delta)
dist_C_frac=np.asarray(dist_C_frac)

n_classes=len(d)
rho = rho_base*np.ones(n_classes) # same resuspension rate for all size classes
C = Cinit*np.ones(n_classes)*dist_C_frac #initialise with the mass distribution
C0 = C # for plotting
L = Linit*np.ones(n_classes)*dist_C_frac  #initialise with the mass distribution
Ce = Ce_base*np.ones(n_classes)*dist_C_frac
# prepare storage for plotting
concentration,deposition,time=[],[],[]
rhoC,deltaC=[],[]
rhop=[]

nb_period=2
period=24
dt=0.1 #h
sim_time=nb_period*48 # hours

t=0 # hour
matrice_C,matrice_L,temps=[],[],[]
while t < sim_time:
	# store for plotting
	matrice_C.append(C)
	matrice_L.append(L)
	temps.append(round(t,2))
	concentration.append(C)
	deposition.append(L)
	time.append(round(t,2))
	rhop.append(rho)
	rhoC.append(rho*S*L/V)
	deltaC.append(delta*C)
	# variable Ce and rho
	Ce=dist_C_frac*(Ce_base+ 5*abs(np.cos(t*2*np.pi/period)))
	rho=rho_base
	# solve for C+ and L+ using the array splitting method
	vec_CL =fsolve(fc_IAQ_coupled_classes, np.hstack([C,L]), args=(np.hstack([C,L]),tau,delta,dt,rho,S,V,Ce))
	C_plus,L_plus = vec_CL[0:n_classes], vec_CL[n_classes:] # split into C and L
	C,L=C_plus,L_plus
	t+=dt


# we want to sum over the size-classes
matrice_C=np.asarray(matrice_C)
matrice_L=np.asarray(matrice_L)
time=len(matrice_C[:,0])
C_tot,L_tot=np.zeros(time),np.zeros(time)
for k in range(time):
	C_tot[k]=sum(matrice_C[k,:])
	L_tot[k]=sum(matrice_L[k,:])

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
plt.show()