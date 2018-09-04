#!/usr/bin/env python -i

import csv, sys, os, array, warnings
import input_sub as input_sub
import numpy as np
from astropy.io import ascii
from astropy import constants as consts
from matplotlib import pyplot as plt
import v_hydro_sub as vhs
from scipy.optimize import brentq
from matplotlib import rc
from matplotlib import gridspec
import pickle
from scipy.interpolate import RectBivariateSpline as RBS
from astropy import units as u
from astropy import constants as c
import plotter_sub2 as ps
from astropy.convolution import convolve, Box1DKernel
import subprocess




		
#t_e,r,theta=ps.py_to_rtheta(ascii.read(root+".te.dat"))
#fe25,r,theta=ps.py_to_rtheta(ascii.read(root+".ioncFe25.dat"))
#fe26,r,theta=ps.py_to_rtheta(ascii.read(root+".ioncFe26.dat"))
#ne,r,theta=ps.py_to_rtheta(ascii.read(root+".ne.dat"))
#density,r,theta=ps.py_to_rtheta(ascii.read(root+".rho.dat"))
#h1,r,theta=ps.py_to_rtheta(ascii.read(root+".ioncH1.dat"))
#h2,r,theta=ps.py_to_rtheta(ascii.read(root+".ioncH2.dat"))
#vx,r,theta=ps.py_to_rtheta(ascii.read(root+".vx.dat"))
#vz,r,theta=ps.py_to_rtheta(ascii.read(root+".vz.dat"))

index=41   #80 
#index=29   #70 

e1=6952.751953
lambda1=1.783442
Z=26
oscillator1=0.139000



e2=6973.968262
lambda2=1.778016
Z=26
oscillator2=0.277000



lmin=1.773
lmax=1.785

fnames=['0_04_edd','0_1_edd','0_3_edd','0_6_edd']


flux=[]
ew_array=[]


fig1=plt.figure()
ax1=fig1.add_subplot(111)
ax2=ax1.twinx()

fig2=plt.figure()
ax3=fig2.add_subplot(111)
ax4=ax3.twinx()



cols=['r','g','b','k','c','m','y','--r','--g','--b','--k','--c','--m','--y']
cols2=['--r','--g','--b','--k','c','m','y','--r','--g','--b','--k','--c','--m','--y']

i=0

for root in fnames:
	print root
	cmdline="py_wind82e "+root+" < pywind_cmds 1> output 2>output2"  
	subprocess.check_call(cmdline,shell=True)

	density,r,theta=ps.py_to_rtheta(ascii.read(root+".rho.dat"))
	fe25,r,theta=ps.py_to_rtheta(ascii.read(root+".ioncFe25.dat"))
	fe26,r,theta=ps.py_to_rtheta(ascii.read(root+".ioncFe26.dat"))
	t_e,r,theta=ps.py_to_rtheta(ascii.read(root+".te.dat"))
	vx,r,theta=ps.py_to_rtheta(ascii.read(root+".vx.dat"))
	vz,r,theta=ps.py_to_rtheta(ascii.read(root+".vz.dat"))
	v_r=np.sqrt(vx**2+vz**2)
	
	ax1.loglog(r,v_r[index]/100./1000.,cols[i])
	ax2.loglog(r,fe25[index],cols2[i])

	ax3.loglog(r,v_r[index]/100./1000.,cols[i])
	ax4.loglog(r,t_e[index],cols2[i])

	i=i+1
	lamb,flux_temp1,EW1=ps.abs_calc(fe26[index],v_r[index],t_e[index],density[index],r,lambda1,Z,oscillator1,lambdamin=lmin,lambdamax=lmax,dlambda=1e-5)
	lamb,flux_temp2,EW2=ps.abs_calc(fe26[index],v_r[index],t_e[index],density[index],r,lambda2,Z,oscillator2,lambdamin=lmin,lambdamax=lmax,dlambda=1e-5)
	print root,EW1+EW2
	
	flux_temp=np.minimum(flux_temp1,flux_temp2)

	flux.append(flux_temp)
	ew_array.append(EW1+EW2)

#labels=['4\% Ledd','10\% Ledd','20\% Ledd','30\% Ledd','40\% Ledd','50\% Ledd','60\% Ledd','70\% Ledd','80\% Ledd','90\% Ledd','100\% Ledd']
labels=['4\% Ledd','10\% Ledd','30\% Ledd','60\% Ledd']

#ps.plot_mult(lamb,[flux[0],flux[1],flux[2],flux[3],flux[4],flux[5],flux[6],flux[7],flux[8],flux[9],flux[10]],lambda0,"80 degrees",lambdamin=lmin,lambdamax=lmax,labels=labels)
ps.plot_mult(lamb,[flux[0],flux[1],flux[2],flux[3]],1.78,"80_degrees_fe26",lambdamin=lmin,lambdamax=lmax,labels=labels)

'''
ax1.set_ylabel("v_r (km/s)")
ax2.set_ylabel("fe25 --")

ax3.set_ylabel("v_r (km/s)")
ax4.set_ylabel("t_e --")

fig1.savefig("v_r.png")
fig2.savefig("t_e.png")
			
			
	
e=6701.269531
lambda0=1.850400
Z=26
oscillator=0.798000

lmin=1.845
lmax=1.852	
	
'''
print 'fin'




		