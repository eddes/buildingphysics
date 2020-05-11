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

def fc_valve_curve(characteristic, valve_pos,Qmax):
	if characteristic == "linear":
		Qalim = Qmax * valve_pos
	elif characteristic == "equal_pct":
		Qalim = Qmax * np.exp(3.5 * (valve_pos - 1))
	elif characteristic == "quadratic":
		Qalim = Qmax * valve_pos**2
	return Qalim

characteristic="linear"
#characteristic="equal_pct"
#characteristics="quadratic"

# Geometry of the problem
Pi = 3.14159269 #nombre pi
Dr = 1 #m tank diameter
Ds = 0.04 #m outlet diameter
Sr = Pi * Dr**2 / 4 #m2 tank surface 
Ss = Pi * Ds**2 / 4 #m2 outlet surface
H0 = 1 #m tank height

# PID parameters
prop_band = 0.04 # m proportional band
Tn = 10 # integration time
Td = 10 # derivation time
Qmax = 100 # L/s max flow rate 
KP = 1 / prop_band # proportional gain 
Qmax = Qmax / 1000 #conversion to m3/s
# initial values
h = 0.1*H0 #m height =set point heit
Qalim = 0  # no flow rate
sim_time= 300 # simulation duration s
# time and timestep
t = 0
dt = 1 # [s] time stem /sampling rate
B = Ss / Sr * np.sqrt(2 * 9.81) # precompute B, constant in the equation for h

#initialise to zero
sum_error = 0
delta_error = 0
valve_position = 0
deltaT_previous = 0

pcent_P,pcent_I,pcent_D=[],[],[]
time,height=[],[]
v_pos=[]

while t <= sim_time:
	#Euler explicit for water height
	A = Qalim / Sr #compute A 
	h = dt * (A - B * np.sqrt(h)) + h
	#integral action
	sum_error = dt / Tn * (H0 - h) + sum_error
	#derivative action
	delta_error = (H0 - h) - deltaT_previous
	#computation of the valve position with PID
	valve_position = KP * ((H0 - h) + sum_error + Td * delta_error / dt)

	#Control for valve opening
	if valve_position < 0:
		valve_position = 0
		Qalim = 0
	elif valve_position > 1:
		valve_position = 1
		Qalim = Qmax
	else:
		Qalim=fc_valve_curve(characteristic, valve_position,Qmax)

	# let's see and plot who does what at each time step
	if valve_position > 0:
		pcent_P.append(abs (H0 - h) * KP / valve_position)
		pcent_I.append( sum_error+(H0 - h) * (dt / Tn) / valve_position)
		pcent_D.append((Td * delta_error / dt) / valve_position)
	else:
		pcent_P.append(0)
		pcent_I.append(0)
		pcent_D.append(0)
	v_pos.append(valve_position)
	
	# compute the difference to the set value for the derivative term
	deltaT_previous = H0 - h 
	t+=dt
	time.append(t)
	height.append(h)

plt.subplot(121)
plt.xlabel("Time [s]")
plt.ylabel("Water height [m]")
plt.plot(time, height,color=coule[-2],linestyle="-",alpha=1,marker='')
plt.plot(time, np.ones(len(time))*H0,color=coule[-1],linestyle="--",alpha=0.85,marker='')

plt.subplot(122)
plt.xlabel("Time [s]")
plt.ylabel("Share of P, I and D in the valve opening [-]")
plt.plot(time, pcent_P,color=coule[0],linestyle="-",alpha=1,marker='',label="P")
plt.plot(time, pcent_I,color=coule[1],linestyle="-",alpha=1,marker='',label="I")
if Td != 0:
	plt.plot(time, pcent_D,color=coule[2],linestyle="-",alpha=1,marker='',label="D")
plt.plot(time, v_pos,color=coule[-1],linestyle="--",alpha=0.95,marker='',label="valve position")

plt.legend()
plt.tight_layout()
plt.show()