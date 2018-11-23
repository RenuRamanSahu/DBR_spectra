#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 17:00:36 2018

@author: dilu
"""

import numpy as np
import matplotlib.pyplot as plt
import core_lib as cl
import E_lib as el
import E_field as ef


 


wl=np.arange(450,1000, 1)
theta_i=0*np.pi/180.0
wl_p=168.26 # plasma wavelength
wl_c=8934.20 # collision wavelength
d_metal=30
metal_pos='l'
#pol='s'  #s polarisation

#Guess parameters
N=8
n1=1.5
n2=2.874
d1=50
d2=90
n_med=1.0

metal_coat=[True, wl_p, wl_c, d_metal, metal_pos]
metal_coat0=[False, wl_p, wl_c, d_metal, metal_pos]
n_input, d_list=cl.n_and_d(n1,n2, d1, d2, N)



#for photonic crystal
wvl_0s, ph_r_0s, R_0s, T_0s=cl.plot_RT_vs_wavelength(wl, 's', theta_i, n_med, n_input, d_list, metal_coat0)
wvl_0p, ph_r_0p, R_0p, T_0p=cl.plot_RT_vs_wavelength(wl, 'p', theta_i, n_med, n_input, d_list, metal_coat0)


#for s polarisation
wvl_s, ph_r_s, R_s, T_s=cl.plot_RT_vs_wavelength(wl, 's', theta_i, n_med, n_input, d_list, metal_coat)

#for p polarisation
wvl_p, ph_r_p, R_p, T_p=cl.plot_RT_vs_wavelength(wl, 'p', theta_i, n_med, n_input, d_list, metal_coat)


fig_spectrum, [ax1, ax2]=plt.subplots(2,1, sharex=True)
#ax1.title('s polarisation spectrum')
#ax2.title('p polarisation spectrum')

ax1.plot(wvl_0s, R_0s, label=' PhC')
ax2.plot(wvl_0p, R_0p, label=' PhC')

ax1.plot(wvl_s, R_s, label=' metal coat')
ax2.plot(wvl_p, R_p, label=' metal coat')
ax2.set_xlabel('wavelength(nm)')
plt.legend()
plt.show()



