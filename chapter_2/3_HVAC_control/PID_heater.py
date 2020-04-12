# -*- coding: utf-8 -*-
from tools.plot_tools import colors
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import math


def fc_valve_curve(caracteristique, valve_pos):
	if caracteristique == "linear":
		Kv_ratio =  valve_pos
	elif caracteristique == "equal_pct":
		Kv_ratio = np.exp(3.5 * (valve_pos - 1))
	elif caracteristique == "quadratic":
		Kv_ratio =  valve_pos**2
	elif caracteristique == "square_root":
		Kv_ratio = np.power(valve_pos,0.5)
	return Kv_ratio 

caracteristique="linear"
#caracteristique="equal_pct"
#caracteristique="quadratic"


# Geometrical data
m_Cp = 1.6 * 10 * 2000 #kJ/K = Cp*V*rho
S_ext = 50 # m2 envelope's surface to the outdoor
U_mur = 1 # W/m2/K
US_ext = U_mur * S_ext
Cp = 4182 # kJ/kg/K pour l#eau

# 
Pmax = 1500 # heat power 
US_radiateur = Pmax / np.power(50,1.3) #DTLM reference = 50K
qm = Pmax / (Cp * 20) # flowrate for 20K temperature difference

Tmax = 60 #Tmax inlet water
T_in = Tmax #
T_out = T_in-10 # initialize T_out with 10K temperature drop
T_set=20
T_init=10
T_amb = T_init
T_ext=5 #outdoor

# PID parameters
prop_band = 10# K
Kp = 1 / prop_band  # passage au gain
Tn = 30
Td = 0
auth=0.5 #valve authority

# initial values
sim_time= 3600 # simulation duration s
# time and timestep
t=0
dt=10 # [s] time step /sampling rate

# initialize PID control values
sum_error = 0
delta_error = 0
valve_position = 0
deltaT_precedent = 0

pcent_P,pcent_I,pcent_D=[],[],[]
time,T_ambl,T_inl,T_outl=[],[],[],[]
heater,losses=[],[]
v_pos=[]

def fc_solve_T_out(T_out,T_in,T_amb):
	DTLM=(T_in-T_out)/math.log((T_in-T_amb)/(T_out-T_amb))
	return  qm*Cp*(T_in-T_out)-US_radiateur*np.power(DTLM,1.3)

P=0.5*Pmax
valve_position=1
deltaT_previous=0
sum_error=0
QA=0
while t <= sim_time:
	# if not heating required
	if T_out==T_amb:
#		print("froid", T_in,T_out,T_amb,QA)
		P_radiateur=0
		# Euler explicit
		T_amb = T_amb + dt / m_Cp * (P_radiateur - US_ext * (T_amb - T_ext))
		# I action
		sum_error = dt / Tn * (T_set - T_amb) + sum_error
		# D action
		delta_error = (T_set - T_amb) - deltaT_previous
		# valve pos
		valve_position = Kp * ((T_set - T_amb) + sum_error + Td * delta_error / dt)
	# heating required
	else:
#		print("chaud", T_in,T_out,T_amb,QA)
		# avoid div by zero
		if (T_in-T_amb)==(T_out-T_amb):
			P_radiateur=qm*Cp*(T_in-T_out)
		else:
			# solve for T_out with initial guess T_out ~ T_in-P/(qm*Cp)
			T_out= fsolve(fc_solve_T_out, T_in-P/(qm*Cp), args=(T_in,T_amb))
			# DTLM
			DTLM=(T_in-T_out)/math.log((T_in-T_amb)/(T_out-T_amb))
			P_radiateur=US_radiateur * np.power(DTLM,1.3)
		# Euler explicit
		T_amb = T_amb + dt / m_Cp * (P_radiateur - US_ext * (T_amb - T_ext))
		# I action
		sum_error = dt / Tn * (T_set - T_amb) + sum_error
		# D action
		delta_error = (T_set - T_amb) - deltaT_previous
		# valve pos
		valve_position = Kp * ((T_set - T_amb) + sum_error + Td * delta_error / dt)
		#Control for valve opening
		if valve_position < 0:
			valve_position = 0
			QA=0
			T_out=T_amb
		elif valve_position > 1:
			valve_position = 1
			QA=1
		else:
			Kv_ratio=fc_valve_curve(caracteristique, valve_position)
			QA = np.power( (auth * np.power(Kv_ratio, -2) + 1 - auth),-0.5)
		# compute the 
		T_in = (QA * Tmax + (1-QA) * T_out) # mixing temperature
		
	v_pos.append(valve_position)
	# let#s see and plot who does what at each time step
	if valve_position > 0:
		pcent_P.append(abs (T_set - T_amb) * Kp / valve_position)
		pcent_I.append( sum_error + Kp*abs (T_set - T_amb)* (dt / Tn) / valve_position)
		pcent_D.append(Kp * (Td * delta_error / dt) / valve_position)
	else:
		pcent_P.append(0)
		pcent_I.append(0)
		pcent_D.append(0)
	
	deltaT_previous = T_set - T_amb
	t+=dt
	time.append(t)
	T_ambl.append(T_amb)
	T_inl.append(T_in)
	T_outl.append(T_out)
	heater.append(US_radiateur * np.power(DTLM,1.3))
	losses.append(US_ext * (T_amb - T_ext))

plt.subplot(121)
plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.plot(time, T_ambl,color=colors[3],linestyle="-",alpha=1,marker='',label='room')
plt.plot(time, T_inl,color=colors[0],linestyle="-",alpha=1,marker='',label='heater inlet')
plt.plot(time, T_outl,color=colors[1],linestyle="-",alpha=1,marker='',label='heater outlet')
plt.plot(time, np.ones(len(time))*T_set,color=colors[-1],linestyle="--",alpha=0.85,marker='')
plt.legend()

plt.subplot(122)
plt.xlabel("Time [s]")
plt.ylabel("Share of P, I and D in the valve opening [-]")
plt.plot(time, pcent_P,color=colors[0],linestyle="-",alpha=1,marker='',label="P")
plt.plot(time, pcent_I,color=colors[1],linestyle="-",alpha=1,marker='',label="I")
if Td != 0:
	plt.plot(time, pcent_D,color=colors[2],linestyle="-",alpha=1,marker='',label="D")


plt.legend()
plt.tight_layout()
plt.savefig("./PID_heater.png",dpi=200,bbox_inches='tight')
