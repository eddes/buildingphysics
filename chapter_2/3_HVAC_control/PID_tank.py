import numpy as np
import matplotlib.pyplot as plt
from tools.plot_tools import colors

def fc_valve_curve(characteristic, valve_pos,Qmax):
	if characteristic == "linear":
		Qalim = Qmax * valve_pos
	elif characteristic == "equal_pct":
		Qalim = Qmax * np.exp(3.5 * (valve_pos - 1))
	elif characteristic == "quadratic":
		Qalim = Qmax * valve_pos**2
	return Qalim

#valve characteristic ="equal_pct" or "quadratic" or "linear"
characteristic="linear" 

# Geometrical data
Pi = 3.14159269 #nombre pi
Dr = 1 #m tank diameter
Ds = 0.04 #m outlet diameter
Sr = Pi * Dr**2 / 4 #m² tank surface 
Ss = Pi * Ds**2 / 4 #m² outlet surface
H0 = 1 #m tank height
B = Ss / Sr * np.sqrt(2 * 9.81) # precompute B, constant in the equation for h

# PID parameters
prop_band = 0.5 # m
Tn = 15
Td = 0
Qmax = 10 # L/s supply flowrate
KP = 1 / prop_band # proportional gain 
Qmax = Qmax / 1000 #conversion to m3/s

# initial values
h = H0 #m height =set point heit
Qalim = 0  # no flow rate
sim_time= 600 # simulation duration s

# time and timestep
t = 0
dt = 0.05 # [s] time stem /sampling rate

#Remise à zero des paramètres de calcul de la regulation
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
	v_pos.append(valve_position)
	# let's see and plot who does what at each time step
	if valve_position > 0:
		pcent_P.append(abs (H0 - h) * KP / valve_position)
		pcent_I.append( sum_error + KP * (H0 - h) * (dt / Tn) / valve_position)
		pcent_D.append(KP * (Td * delta_error / dt) / valve_position)
	else:
		pcent_P.append(0)
		pcent_I.append(0)
		pcent_D.append(0)
	deltaT_previous = H0 - h
	t+=dt
	time.append(t)
	height.append(h)

plt.subplot(121)
plt.xlabel("Time [s]")
plt.ylabel("Water height [m]")
plt.plot(time, height,color=colors[-2],linestyle="-",alpha=1,marker='')
plt.plot(time, np.ones(len(time))*H0,color=colors[-1],linestyle="--",alpha=0.85,marker='')

plt.subplot(122)
plt.xlabel("Time [s]")
plt.ylabel("Share of P, I and D in the valve opening [-]")
plt.plot(time, pcent_P,color=colors[0],linestyle="-",alpha=1,marker='',label="P")
plt.plot(time, pcent_I,color=colors[1],linestyle="-",alpha=1,marker='',label="I")
if Td != 0:
	plt.plot(time, pcent_D,color=colors[2],linestyle="-",alpha=1,marker='',label="D")
plt.plot(time, v_pos,color=colors[-1],linestyle="--",alpha=0.95,marker='',label="valve position")

plt.legend()
plt.tight_layout()
plt.savefig("./PID_tank.png",dpi=200,bbox_inches='tight')
