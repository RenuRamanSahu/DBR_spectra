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


 


wl=np.arange(400,825, 1)
theta_i=60*np.pi/180.0
wl_p=168.26 # plasma wavelength
wl_c=8934.20 # collision wavelength
d_metal=30
metal_pos='l'
pol='s'  #s polarisation




#Guess parameters
N=8
n1=1.5
n2=2.874
d1=50
d2=90
n_med=1.0

metal_coat=[True, wl_p, wl_c, d_metal, metal_pos]
n_input, d_list=cl.n_and_d(n1,n2, d1, d2, N)




def get_wvl(event):
    wvl=event.xdata # gets wavelength
    print "Plotting Intensity at", int(wvl), " nm ..."
    ef.plot_field(wvl, pol, theta_i, n_med, n_input, d_list, metal_coat)    
    

#for s polarisation

wvl, ph_r, R, T=cl.plot_RT_vs_wavelength(wl, 's', theta_i, n_med, n_input, d_list, metal_coat)

#for p polarisation
#wvl_p, ph_r_p, R_p, T_p=cl.plot_RT_vs_wavelength(wl, 'p', theta_i, n_med, n_input, d_list, metal_coat)


fig_spectrum=plt.figure()
plt.plot(wvl, R, label='s polarisation')
cid=fig_spectrum.canvas.mpl_connect('button_press_event', get_wvl)

#plt.plot(wvl_p, R_p, label='p polarisation')
plt.legend()
plt.show()



"""
#Plots Trnsmittivity
plt.figure()
plt.plot(wvl, T)
plt.show()


"""
