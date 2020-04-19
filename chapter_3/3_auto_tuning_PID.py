# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:15:51 2020

@author: walthere
"""
#quelques couleurs
rouge_A='#C60C2E'
vert1_A='#005157'
vert2_A='#627D77'
vert3_A='#9EB28F'
vert4_A='#C5E5A4'
medblue='royalblue'
gris1_A='#595A5C'
coule=[rouge_A,vert1_A,vert2_A,vert3_A,vert4_A,medblue,gris1_A]
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import minimize

def fc_CN(hp,h,dt,A,B):
	return  -hp + h +0.5*(dt*(A-B*np.sqrt(h)) ) + 0.5*( dt*(A-B*np.sqrt(hp)) )


def fc_valve_curve(characteristic, valve_pos,Qmax):
	if characteristic == "linear":
		Qsupply = Qmax * valve_pos
	elif characteristic == "equal_pct":
		Qsupply = Qmax * np.exp(3.5 * (valve_pos - 1))
	elif characteristic == "quadratic":
		Qsupply = Qmax * valve_pos**2
	return Qsupply

characteristic="linear"
#characteristic="equal_pct"
#characteristics="quadratic"


# Geometrical data
Pi = 3.14159269 #nombre pi
Dr = 1 #m tank diameter
Ds = 0.04 #m outlet diameter
Sr = Pi * Dr**2 / 4 #m² tank surface 
Ss = Pi * Ds**2 / 4 #m² outlet surface
H0 = 1 #m tank height

# PID parameters
prop_band = 0.5 # m
Kp=1/prop_band
Tn = 500
Td = 0
Qmax = 15 # L/s
KP = 1 / prop_band # proportional gain 
Qmax = Qmax / 1000 #conversion to m3/s
f_init=0.
# time and timestep
t = 0
dt=1 # [s] time stem /sampling rate
B = Ss / Sr * np.sqrt(2 * 9.81) # precompute B, constant in the equation for h

def tu_bluffes_martoni(x):
	BP,Tn,Td=x[0],x[1],x[2]
	Kp=1/BP
	# initial values
	h = f_init*H0 #m height =set point height
	Qsupply = 0  # no flow rate
	sim_time=200 # simulation duration s
	time,height,v_pos=[],[],[]
	t=0
	#Remise à zero des paramètres de calcul de la regulation
	sum_error = 0
	delta_error = 0
	valve_position = 0
	deltaT_previous = 0
	
	time.append(t)
	height.append(h)
	v_pos.append(valve_position )
	
	while t <= sim_time:
		#Euler explicit for water height
		A = Qsupply / Sr #compute A 
		h=fsolve(fc_CN, h, args=(h,dt,A,B))
		#integral action
		sum_error = dt / Tn * (H0 - h) + sum_error
		#derivative action
		delta_error = (H0 - h) - deltaT_previous
		#computation of the valve position with PID
		valve_position = Kp * ((H0 - h) + sum_error + Td * delta_error / dt)
	
		#Control for valve opening
		if valve_position < 0:
			valve_position = 0
			Qsupply = 0
		elif valve_position > 1:
			valve_position = 1
			Qsupply = Qmax
		else:
			Qsupply=fc_valve_curve(characteristic, valve_position,Qmax)
		deltaT_previous = H0 - h
		t+=dt
		time.append(t)
		height.append(h)
		v_pos.append(valve_position)
	# computation of the objective function
	v_pos=np.asarray(v_pos)
	height_set= H0*np.ones(len(height))
	diff=np.asarray(height)-height_set
	diff=abs(diff)
	mean_diff=np.mean(diff)
	return mean_diff, time, height, v_pos

def fc_to_minimize(x):
	BP,Tn,Td=x[0],x[1],x[2]
	Kp=1/BP
	# initial values
	h = f_init*H0 #m height =set point height
	Qsupply = 0  # no flow rate initially
	sim_time=200 # simulation duration s
	time,height,v_pos=[],[],[]
	t=0
	#Remise à zero des paramètres de calcul de la regulation
	sum_error = 0
	delta_error = 0
	valve_position = 0
	deltaT_previous = 0
	while t <= sim_time:
		#Crank-Nicolson semi-implicit for water height
		A = Qsupply / Sr #compute A 
		h=fsolve(fc_CN, h, args=(h,dt,A,B))
		#integral action
		sum_error = dt / Tn * (H0 - h) + sum_error
		#derivative action
		delta_error = (H0 - h) - deltaT_previous
		#computation of the valve position with PID
		valve_position = Kp * ((H0 - h) + sum_error + Td * delta_error / dt)
#		print(valve_position)
		#Control for valve opening
		if valve_position < 0.:
			valve_position = 0
			Qsupply = 0
		elif valve_position > 1:
			valve_position = 1
			Qsupply = Qmax
		else:
			Qsupply=fc_valve_curve(characteristic, valve_position,Qmax)
		deltaT_previous = H0 - h
		t+=dt
		time.append(t)
		height.append(h)
		v_pos.append(valve_position)
	v_pos=np.asarray(v_pos)
	
	# penalize the result if 
	penalty=0
	idx=np.where(v_pos==0)
	if len(idx[0])>10:penalty=100
		
	height_set= H0*np.ones(len(height))
	diff=np.asarray(height)-height_set
	diff=abs(diff)
	mean_diff=np.mean(diff)
	return mean_diff+penalty

#           Bp      integration  derivation
bnds = ((0.01, 0.6), (1, 400) , (10, 40))
k=0.5
x0=[k*bnds[0][0]+bnds[0][1],
    k*bnds[1][0]+bnds[1][1],
    k*bnds[2][0]+bnds[2][1]]
#sol = minimize(fc_to_minimize, x0, bounds=bnds, method='Nelder-Mead', tol=1e-3)
#sol = minimize(fc_to_minimize, x0, bounds=bnds, method='SLSQP', tol=1e-3)
#sol = minimize(fc_to_minimize, x0, bounds=bnds, method='L-BFGS-B', tol=1e-3)
sol = minimize(fc_to_minimize, x0, bounds=bnds, method='TNC', tol=1e-3)


BP,Tn,Td=sol.x[0],sol.x[1],sol.x[2]
objective, time, height, v_pos=tu_bluffes_martoni([BP,Tn,Td])
print("BP =",round(BP,2))
print("Tn =",round(Tn,0))
print("Td =", round(Td,0))
print("fc obj =", round(objective[0],3))

#plt.subplot(121)
#plt.ylim(0,1)
plt.xlabel("Time [s]")
plt.ylabel("Water height [m]")
plt.plot(time, height,color=coule[-2],linestyle="-",alpha=1,marker='',label='water height')
plt.plot(time, v_pos,color=coule[-1],linestyle="--",alpha=1,marker='',label='valve position')
plt.plot(time, np.ones(len(time))*H0,color='black',linestyle="--",alpha=0.9985,marker='',label='set water height')

plt.legend()
plt.tight_layout()
plt.savefig("./PID_optimum_values.pdf")