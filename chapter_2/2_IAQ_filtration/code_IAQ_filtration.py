# -*- coding: utf-8 -*-
from tools.plot_tools import colors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#this function computes the average week from a "dataframe" (a better routine may well exist!)
def fc_semaine_type_df(df, type_analyse):
	if type_analyse=="semaine":
		Cmoy_semaine=np.zeros(168)
		for i in range(7):
			df_temp=pd.DataFrame()
			df_temp=df[df.index.dayofweek==i] # dont forget square brackets
			j_moy_chaque=df_temp.groupby(df_temp.index.hour).mean()
			j_moy=j_moy_chaque.mean(axis=1) #average calculation
			for j in range(24):
				Cmoy_semaine[i*24+j]=round(j_moy[j],2) # rounding 
		return Cmoy_semaine
	else:
		j_moy=df.groupby(df.index.hour).mean()
		return j_moy.mean(axis=1)

# compute average day/week from a "Serie"
def fc_semaine_type_series(df,type_analyse):
	if type_analyse=="semaine":
		Cmoy_semaine=np.zeros(168)
		for i in range(7):
			df_temp=pd.DataFrame()
			df_temp=df[df.index.dayofweek==i]
			j_moy=df_temp.groupby(df_temp.index.hour).mean()
			for j in range(24):
				Cmoy_semaine[i*24+j]=round(j_moy[j],2) # rounding
		return Cmoy_semaine
	else:
		df_temp=pd.DataFrame() # to perform stats
		df_temp=df #copy the original df
		return df_temp.groupby(df_temp.index.hour).mean() # hourly average

#  returns the filter efficiency
def fc_eta(m_filt): # m_filt given in grams
	if m_filt !=0:
		return (-0.0117*m_filt**2- 0.3033*m_filt+ 80.4)/100
	else: return 0.8

# returns the pressure drop
def fc_pressure_drop(m_filt):
	return 0.145*m_filt**2 - 0.81*m_filt+ 47.2

tau=1# vol/h
delta=0.15#	1/h
eta=0.8 # initial efficiency
V=2000 #m3
qv=2000 #m3/h
dt=1	#h

m_limit=30 # max mass before maintenance (g)
nb_filt=0 # number of filters used
# load the PM2.5 data
Cext=np.loadtxt("PM25.txt")
nb=len(Cext)
time=np.arange(0,nb,1)
# prepare the vectors for computing
Csupply=np.zeros(nb)
C=np.zeros(nb)
m_filter=np.zeros(nb) # initialy nothing in filter
eta_filter=np.ones(nb)*eta # initial eta for eta_filter[0] actually (will be updated in the time loop)
pdc_filter=np.zeros(nb)

Csupply[0]=Cext[0] # aesthetical fill for plotting
m_filter[0]=1e-4 # against division by !0
pdc_ref=fc_pressure_drop(0) # reference pressure drop for clean filter
# time loop
for i in range(1,nb-1):
	# compute efficiency
	eta_filter[i]=fc_eta(m_filter[i-1])
	# compute the additionnal pressure drop compared to clean filter
	pdc_filter[i]=fc_pressure_drop(m_filter[i-1])-pdc_ref 
	# supply PM2.5 concentration
	Csupply[i]=(1-eta_filter[i])*Cext[i]
	# mass in filter, converted to grams
	m_filter[i]=(Cext[i]-Csupply[i])*qv/1e6 + m_filter[i-1]
	# check if we change the filter
	if m_filter[i]>m_limit:
		m_filter[i]=0
		nb_filt+=1 # count the number of filter change 
	#explicit scheme for C+
	C[i]=C[i-1] + dt*tau*(Csupply[i]-C[i-1])
# compute the energy use 
watts_hourly=pdc_filter*qv/3600
total_kWh=sum(watts_hourly)/1000

# post process and plotting
a=pd.date_range('1/1/2018', periods=8760, freq='h')
dfC = pd.DataFrame(index=a)

dfC["PM25"]=Csupply
dfC["PM25_ext"]=Cext
dfC["eta"]=eta_filter
dfC["mass"]=m_filter
dfC["pdc"]=pdc_filter

# get rid off the last empty index
dfC=dfC.drop(dfC.index[-1])

typical_week_filt=fc_semaine_type_series(dfC["PM25"],"semaine")
typical_week_ext=fc_semaine_type_series(dfC["PM25_ext"],"semaine")

plt.clf()
plt.xlabel("Average week [h]]")
plt.ylabel(r"$PM_{2.5}$[µg/m$^3$]")
plt.plot(typical_week_ext, color=colors[-1], linestyle="-", alpha=0.9, marker='', label='outdoor')
plt.plot(typical_week_filt, color=colors[0], linestyle="-", alpha=0.9, marker='', label='indoor')
plt.legend()
plt.savefig("./filter_typical_ouik.pdf",dpi=200,bbox_inches='tight')

plt.clf()
plt.ylabel("Efficiency [-]")
dfC["eta"].plot(color=colors[2], linestyle="--", alpha=1, marker='')
plt.savefig("./filter_eta.pdf",dpi=200,bbox_inches='tight')

plt.clf()
plt.ylabel("Mass in filter [g]")
dfC["mass"].plot(color=colors[3], linestyle="--", alpha=1, marker='')
plt.savefig("./filter_mass.pdf",dpi=200,bbox_inches='tight')

plt.clf()
plt.ylabel(r"Additional pressure drop of filter [$\Delta$ Pa]")
dfC["pdc"].plot(color=colors[4], linestyle="--", alpha=1, marker='')
plt.savefig("./filter_pdc.pdf",dpi=200,bbox_inches='tight')

plt.clf()
plt.ylabel(r"$PM_{2.5}$[µg/m$^3$]")
dfC["PM25"].plot(color=colors[-1], linestyle="--", alpha=0.5, marker='')
plt.savefig("./PM25_year.pdf",dpi=200,bbox_inches='tight')
