#quelques couleurs
rouge_A='#C60C2E'
vert1_A='#005157'
vert2_A='#627D77'
vert3_A='#9EB28F'
vert4_A='#C5E5A4'
medblue='royalblue'
gris1_A='#595A5C'
coule=[rouge_A,vert1_A,vert2_A,vert3_A,vert4_A,medblue,gris1_A]

from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import math
import sys

# function for Crank-Nicolson's scheme
def fc_CN(Tp,T, dt, m_Cp, P_heater,US_ext, T_amb,T_ext):
	return  -Tp +T +0.5*(dt/m_Cp*(P_heater-US_ext*(T-T_ext))) + 0.5*(dt/m_Cp*(P_heater-US_ext*(Tp-T_ext)))


def fc_valve_curve(characteristic, valve_pos):
	if characteristic == "linear":
		Kv_ratio =  valve_pos
	elif characteristic == "equal_pct":
		Kv_ratio = np.exp(3.5 * (valve_pos - 1))
	elif characteristic == "quadratic":
		Kv_ratio =  valve_pos**2
	elif characteristic == "square_root":
		Kv_ratio = np.power(valve_pos,0.5)
	return Kv_ratio 

#characteristic="linear"
characteristic="equal_pct"
#characteristic="quadratic"


# Geometrical data
m_Cp = 1006 * 5*5*2.5* 1 #kJ/K = Cp*V*rho
S_ext = 50 #m2 envelope surface to outdoors
U_wall = 1 #W/m2/K
US_ext = U_wall*S_ext
Cp = 4182 #kJ/kg/K water heat capacity

# 
Pmax = 3000 # maximum heater power 
US_heater = Pmax / np.power(50,1.3) #DTLM nominal = 50K
qm = Pmax / (Cp * 20) # flowrate

Tmax = 60 #Tmax eau
T_in = Tmax #
T_out = T_in-5 #initialisation of the temperature heater outlet temp.
T_set=20
T_init=16
T_amb = T_init
T_ext=5 #outdoor

# PID parameters
prop_band =5# K
Kp = 1/prop_band  # passage au gain
Tn = 3000
Td = 50
auth=0.5 #valve authority

# initial values
sim_time= 3600 # simulation duration s
# time and timestep
t=0
dt=10 # [s] time step /sampling rate

#Remise à zero des paramètres de calcul de la regulation
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
	return  qm*Cp*(T_in-T_out)-US_heater*np.power(DTLM,1.3)

# initializations
P_heater=0.5*Pmax 
valve_position=1
deltaT_previous=0
delta_error=0
sum_error=0
QA=0
# time loop
while t <= sim_time:
	# if not heating required
	if T_amb>T_set:
		QA=0
		T_out=T_amb  #aesthetical fill
		T_in = T_out #aesthetical fill
		P_heater=0
		# Crank Nicolson
		T_amb=fsolve(fc_CN, T_amb, args=(T_amb, dt, m_Cp, P_heater,US_ext, T_amb,T_ext))
		# I action
		sum_error = dt / Tn * (T_set - T_amb) + sum_error
		# D action
		delta_error = (T_set - T_amb) - deltaT_previous
		# valve pos
		valve_position = Kp * ((T_set - T_amb) + sum_error + Td * delta_error / dt)
		# store the absolute contribution of all actions for plotting:
		v_share = Kp * (abs(T_set - T_amb) + abs(sum_error) + Td *abs( delta_error )/ dt)
		#Control for valve opening
		if valve_position < 0:
			valve_position = 0
			QA=0
		elif valve_position > 1:
			valve_position = 1
			QA=1
		else:
			#valve curve
			Kv_ratio=fc_valve_curve(characteristic, valve_position)
			#flow rate
			QA = np.power( (auth * np.power(Kv_ratio, -2) + 1 - auth),-0.5)
	# heating required
	else:
		T_out= fsolve(fc_solve_T_out, T_in-P_heater/(qm*Cp), args=(T_in,T_amb))
		# DTLM
		DTLM=(T_in-T_out)/math.log((T_in-T_amb)/(T_out-T_amb))
		# avoid div by zero and NaN when heating required after non-heating period
		if np.isnan(DTLM):
			T_out=max((T_in-0.1),(T_amb-0.1)) # well, you gotta start somewhere...
			DTLM=(T_in-T_out)/math.log((T_in-T_amb)/(T_out-T_amb)) # avoid math error
		P_heater=US_heater * np.power(DTLM,1.3)

		# Crank Nicolson
		T_amb=fsolve(fc_CN, T_amb, args=(T_amb, dt, m_Cp, P_heater,US_ext, T_amb,T_ext))
		# I action
		sum_error = dt / Tn * (T_set - T_amb) + sum_error
		# D action
		delta_error = (T_set - T_amb) - deltaT_previous
		# valve pos
		valve_position = Kp * ((T_set - T_amb) + sum_error + Td * delta_error / dt)
		v_share = Kp * (abs(T_set - T_amb) + abs(sum_error) + Td *abs( delta_error )/ dt)
		
		#Control for valve opening
		if valve_position < 0:
			valve_position = 0
			QA=0
			T_out=T_amb
		elif valve_position > 1:
			valve_position = 1
			QA=1
		else:
			# valve curve
			Kv_ratio=fc_valve_curve(characteristic, valve_position)
			# flow rate
			QA = np.power( (auth * np.power(Kv_ratio, -2) + 1 - auth),-0.5)
	
	# compute the inlet temperature
	T_in = (QA * Tmax + (1-QA) * T_out) # mixing temperature

#		sys.exit(0)
	v_pos.append(valve_position)
	# let#s see and plot who does what at each time step
	if valve_position > 0:
		pcent_P.append(abs (T_set - T_amb) * Kp / v_share)
		pcent_I.append(Kp*(abs(sum_error))/ v_share)
		pcent_D.append(Kp * (Td *abs(delta_error) / dt) / v_share)
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
	heater.append(P_heater) # US_heater * np.power(DTLM,1.3)
	losses.append(US_ext * (T_amb - T_ext))

plt.subplot(121)

plt.xlabel("Time [s]")
plt.ylabel("Temperature [°C]")
plt.plot(time, T_ambl,color=coule[3],linestyle="-",alpha=1,marker='',label='room')
plt.plot(time, T_inl,color=coule[0],linestyle="-",alpha=1,marker='',label='heater inlet')
plt.plot(time, T_outl,color=coule[1],linestyle="-",alpha=1,marker='',label='heater outlet')
plt.plot(time, np.ones(len(time))*T_set,color=coule[-1],linestyle="--",alpha=0.85,marker='')
plt.legend()

plt.subplot(122)
plt.xlabel("Time [s]")
plt.ylabel("Share of P, I and D in the valve opening [-]")
plt.plot(time[1:], v_pos[1:],color=coule[-1],linestyle="--",alpha=1,marker='',label="valve position")
plt.plot(time[1:], pcent_P[1:],color=coule[0],linestyle="-",alpha=1,marker='',label="P")
plt.plot(time[1:], pcent_I[1:],color=coule[1],linestyle="-",alpha=1,marker='',label="I")
if Td != 0:
	plt.plot(time[1:], pcent_D[1:],color=coule[2],linestyle="-",alpha=1,marker='',label="D")

plt.legend()
plt.tight_layout()
plt.show()
