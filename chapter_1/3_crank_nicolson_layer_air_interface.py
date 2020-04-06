import numpy as np
import matplotlib.pyplot as plt
from  scipy.optimize import fsolve

# function for Crank-Nicolson's scheme
def fc_CN(Tp,T,K,beta):
	return  Tp -T - beta*np.dot(K,T) - beta*np.dot(K,Tp)

beta=0.5
n_solid=19
n=n_solid+2
t=0
hours=4
sim_time=hours*3600
dt=50

L_tot=0.9 # m
dx=L_tot/(n_solid) # m
x_pos=np.arange(0,L_tot+2*dx,dx)

# material 1 and 2 (1 = on the side of i=0)
lambda1,rho1,Cp1 =1.6,2200,1100 # 8.1*1e-7
lambda2,rho2,Cp2 =  0.04,300,900 # 1.2*1e-7

alpha1 = lambda1/(rho1*Cp1)
alpha2 = lambda2/(rho2*Cp2)

# fourier for both layers
Fo1=alpha1*dt/dx**2
Fo2=alpha2*dt/dx**2

# fourier-biot
hi=7.7 # W/m2/K heat transfer side i
FoBi_i=dt/(dx*rho1*Cp1)*np.power((1/hi+dx/lambda1/2),-1)
he=25 # W/m2/K heat transfer side e
FoBi_e=dt/(dx*rho2*Cp1)*np.power((1/he+dx/lambda2/2),-1)

Tmin=-5
Tmax=20
T_init=Tmin #0.5*(Tmin+Tmax)
T=np.ones(n)*T_init
T[0]=Tmin
T[n-1]=Tmax

Fo=np.zeros(n)
L_layer1= 3*L_tot/4 + 0.003 # add epsilon to avoid falling straight on a node

# for the interface property
x_interface=L_tot-L_layer1
i_interface=int(x_interface/dx)
dx1=max(i_interface*dx-x_interface,x_interface-i_interface*dx)
dx2=dx-dx1
# lambda_eq at the interface
k_eq=1/(dx1/lambda1+dx2/lambda2)
# Fo_eq
Fo_eq=dt*k_eq/(rho1*Cp1*dx**2)

# let's build the conductivity matrix
K=np.eye(n,n,k=-1)*1 + np.eye(n,n)*-2+ np.eye(n,n,k=1)*1
K[0,0],K[0,1],K[-1,-1],K[-1,-2]=0,0,0,0 # ghost lines for BC

coeffs_Fo=np.ones(len(T))
coeffs_Fo[0:i_interface+1]=Fo1 # prepare the upper part of the matrix
coeffs_Fo[i_interface+1:]=Fo2 # ... the lower part of the matrix
K=coeffs_Fo*K
# diagonal terms around interface
K[i_interface,i_interface]=-(Fo1+Fo_eq)
K[i_interface+1,i_interface+1]=-(Fo2+Fo_eq)
# sub and supra diag terms around interface
K[i_interface,i_interface+1]=Fo_eq
K[i_interface+1,i_interface]=Fo_eq
# exchange with air
K[1,0]=FoBi_i
K[1,1]=-FoBi_i-Fo1
# exchange with air
K[-2,-1]=FoBi_e
K[-2,-2]=-FoBi_e-Fo2


t=0
while t < sim_time:
	T_plus = fsolve(fc_CN, T, args=(T,K,beta))
	T = T_plus
	t+=dt

plt.xlabel("x position [m]")
plt.ylabel("Temperature [°C]")
plt.plot(x_pos, T, 'o-', alpha=0.65,label='crank-nicolson')
plt.legend()
plt.arrow(L_layer1,Tmin,0,abs(Tmax-Tmin))










