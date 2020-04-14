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

# governing equations as a reminder
#Cp = C + dt*( tau*(Ce-C) -delta*C+ rho*S*L/V)
#Lp = L + dt*( rho*V*C/S - rho*S*L/S )
def fc_IAQ_coupled(vec_CLp, vec_CL, tau,delta,dt,rho,S,V,Ce):
	C,L = vec_CL[0],vec_CL[1]
	Cp,Lp = vec_CLp[0],vec_CLp[1]
	C_term = -Cp + C + 0.5*dt*( tau*(Ce-C) -delta*C+ rho*S*L/V) + 0.5*dt*( tau*(Ce-Cp) -delta*Cp+ rho*S*Lp/V)
	L_term = -Lp + L + 0.5*dt*( delta*V*C/S - rho*S*L/S ) + 0.5*dt*( delta*V*Cp/S - rho*S*Lp/S )
	return [C_term,L_term]

# enclosure
L,l,h=5,5,3# dimensions
S=2*( L*l+l*h+ h*L) # m2
V=L*l*h # m3 
tau=0.1# vol/h
qv=V*tau#m3/h
tau=qv/V

Cinit=20
Linit=100
Ce_base=20
delta=0.15#	1/h
rho_base=0.5 # initial efficiency

C,L=Cinit,Linit
concentration,deposition,time=[],[],[]
rhoC,deltaC=[],[]
rhop=[]

nb_period=2
period=24
dt=0.1 #h
sim_time=nb_period*24 # hour
t=0 # hour

while t < sim_time:
	Ce=Ce_base + 5*(np.cos(t*2*np.pi/period))
	rho=rho_base*abs(np.sin(t*2*np.pi/period))
	C_plus,L_plus=fsolve(fc_IAQ_coupled, [C,L], args=([C,L],tau,delta,dt,rho,S,V,Ce))
	C,L=C_plus,L_plus
	t+=dt
	concentration.append(C)
	deposition.append(L)
	time.append(t)
	rhop.append(rho)
	rhoC.append(rho*S*L/V)
	deltaC.append(delta*C)
	
	
plt.subplot(121)
plt.xlabel("Time [h]")
plt.ylabel(r"Concentration [µg/m$^3$]")
plt.plot(time, concentration, '-',color=coule[0], alpha=0.65,label='C')
#plt.plot(time, rhop, '--',color=coule[1], alpha=0.65,label='rho')
plt.legend()

plt.subplot(122)
plt.xlabel("Time [h]")
plt.ylabel(r"Mass on surfaces [µg/m$^2$]")
plt.plot(time, deposition, '-',color=coule[-1], alpha=0.65,label='L')
plt.legend()

plt.tight_layout()


plt.clf()
plt.subplot(121)
plt.xlabel("Time [h]")
plt.ylabel(r"Resuspension rate [µg/m$^3$]")
plt.plot(time, rhop, '-',color=coule[0], alpha=0.65,label=r'$\rho$')
plt.legend()

plt.subplot(122)
plt.xlabel("Time [h]")
plt.ylabel(r"Transfer to air [µg/m$^{3}$] and surfaces [µg/m$^{2}$]")
plt.plot(time, rhoC, '--',color=coule[3], alpha=0.65,label=r'$ \frac{\rho L S}{V}$')
plt.plot(time, deltaC, '-',color=coule[2], alpha=0.65,label=r'$\delta C$')
plt.legend()

plt.tight_layout()
