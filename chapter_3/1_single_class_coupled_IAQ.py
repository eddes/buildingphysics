# -*- coding: utf-8 -*-
# a few colors
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

# enclosure features
# dimensions
L,l,h=5,5,3 #m3
# surface
S=2*( L*l+l*h+ h*L) # m2
#volume
V=L*l*h # m3
# air change rate
tau=0.1# vol/h
# flow rate
qv=V*tau#m3/h

# definitions
Cinit=20 # initial C [\mu g/m3]
Linit=100 # initial L [\mu g/m2]
Ce_base=20 # base outdoor concentratio n
delta=0.15#	1/h deposition rate
rho_base=0.5 # base resuspension rate

# initialize
C,L=Cinit,Linit
concentration,deposition,time=[],[],[] # for storage and plotting
rhoC,deltaC=[],[] # same
rhop=[] # same same

nb_period=2 # how many periods are to be simulated ?
period=24 # for the periodic resuspension rate
dt=0.1 # h (this time dt is in [hours])
sim_time=nb_period*24 # hour
t=0 # initialize
# time loop
while t < sim_time:
	Ce=Ce_base + 5*(np.cos(t*2*np.pi/period)) # update outdoor concentration
	rho=rho_base*abs(np.sin(t*2*np.pi/period)) # update resuspension rate
	# compute C+ and L+
	C_plus,L_plus=fsolve(fc_IAQ_coupled, [C,L], args=([C,L],tau,delta,dt,rho,S,V,Ce))
	C,L=C_plus,L_plus
	t+=dt
	# store for plotting
	concentration.append(C)
	deposition.append(L)
	time.append(t)
	rhop.append(rho)
	rhoC.append(rho*S*L/V)
	deltaC.append(delta*C)
	
# plotting
plt.subplot(221)
plt.xlabel("Time [h]")
plt.ylabel(r"Concentration [µg/m$^3$]")
plt.plot(time, concentration, '-',color=coule[0], alpha=0.65,label='C')
plt.legend()

plt.subplot(222)
plt.xlabel("Time [h]")
plt.ylabel(r"Mass on surfaces [µg/m$^2$]")
plt.plot(time, deposition, '-',color=coule[-1], alpha=0.65,label='L')
plt.legend()

plt.subplot(223)
plt.xlabel("Time [h]")
plt.ylabel(r"Resuspension rate [µg/m$^3$]")
plt.plot(time, rhop, '-',color=coule[0], alpha=0.65,label=r'$\rho$')
plt.legend()

plt.subplot(224)
plt.xlabel("Time [h]")
plt.ylabel(r"Transfer to air and surfaces [µg/m$^{3}/s$]")
plt.plot(time, rhoC, '--',color=coule[3], alpha=0.65,label=r'$ \frac{\rho L S}{V}$')
plt.plot(time, deltaC, '-',color=coule[2], alpha=0.65,label=r'$\delta C$')
plt.legend()

plt.tight_layout()
plt.show()
