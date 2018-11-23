#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 15:07:01 2018

@author: dilu
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 19:07:38 2018

@author: dilu
"""

import numpy as np
import matplotlib.pyplot as plt
import core_lib as cl
#For s polarisation
#input Ef and Eb

# Required point inputs
# 1) wavelength
# 2) angle of incidence
# 3) Amplitude of forward going wave
# 4) Amplitude of backward going wave



n_med=1.0



def x_and_n(n_input, d_input):
    if len(n_input)!=len(d_input):
        print('x and n list should be of same length')
        return 0
    ni=[]
    for i in range(len(n_input)):
        ni=np.concatenate([ni,n_input[i]*np.ones(d_input[i])])
    
    nf=n_med*np.ones(10)
    nb=n_med*np.ones(1)
    n=np.concatenate([nb, ni, nf])
    x=np.arange(0, len(n), 1)
    return x, n

def get_theta(x_list, n_list, theta_i):
    n=n_list
    theta_lst=np.zeros(len(x_list), dtype=complex)
    theta_lst[0]=theta_i
    for i in np.arange(1, len(x_list)-1, 1):
        theta_lst[i]=np.arcsin((n[i]/n[i+1])*np.sin(theta_lst[i-1]))
    return theta_lst
    

def detect_interface(n_list):
    detect_i=np.zeros(len(n_list))
    for i in range(len(n_list)-1):
        if n_list[i]==n_list[i+1]:
            detect_i[i]=0
        else:
            detect_i[i]=1
    return detect_i
    
def BC(n1, n2, theta1, theta2):
    BC_array1=[[1.0, 1.0],[n1*np.cos(theta1), -n1*np.cos(theta1) ]]
    BC_array2=[[1.0, 1.0],[n2*np.cos(theta2), -n2*np.cos(theta2) ]]
    D=np.matmul(np.linalg.inv(BC_array2), BC_array1)
    return D





def R_and_T(wl, theta_i, n_input, d_input):
    x_list, n_list= x_and_n(n_input, d_input)
    is_interface=detect_interface(n_list)
    dx=x_list[1]-x_list[0]

    theta=np.zeros(len(x_list), dtype=complex)
    theta=get_theta(x_list, n_list, theta_i)

    k_list=(2*np.pi/wl)*n_list
    kx_list=k_list*np.cos(np.array(theta))

    phase_shift=np.zeros(len(x_list), dtype=complex)
    phase_shift=kx_list*dx

#PP stands for phase propagation
    PP_list=np.zeros((len(x_list), 2, 2), dtype=complex)
    for i in range(len(x_list)):
        PP_list[i]=[[np.exp((1j)*phase_shift[i]), 0],[0, np.exp((-1j)*phase_shift[i])]]
    
#transfer matrix at ith point
    Tmi=np.zeros((len(x_list), 2, 2), dtype=complex)
    Tmi=PP_list
    for i in range(len(Tmi)-1):
        if is_interface[i]==1:
            B_condition=BC(n_list[i], n_list[i+1], theta[i], theta[i+1])
            Tmi[i]=np.matmul(B_condition, PP_list[i])


    TM=[[1.0, 0], [0.0, 1.0]]
    for i in range(len(Tmi)):
        TM=np.matmul(TM, Tmi[i])
    
    r=TM[1][0]/TM[0][0]
    R=np.abs(r)**2
    t=1/TM[0][0]        
    T=np.abs(t)**2
    return R, T
    







def get_E_field(wl,pol,  theta_i, n_input, d_input, Ef, Eb):
    x_list, n_list= x_and_n(n_input, d_input)
    is_interface=detect_interface(n_list)
    dx=x_list[1]-x_list[0]

    Ep=np.zeros(len(x_list), dtype=complex)
    Ep[0]=Ef
    Em=np.zeros(len(x_list), dtype=complex)
    Em[0]=Eb
    
    theta=np.zeros(len(x_list), dtype=complex)
    theta=cl.get_theta(theta_i, n_med, n_list)

    k_list=(2*np.pi/wl)*n_list
    kx_list=k_list*np.cos(np.array(theta))

    phase_shift=np.zeros(len(x_list), dtype=complex)
    phase_shift=kx_list*dx

#PP stands for phase propagation
    PP_list=np.zeros((len(x_list), 2, 2), dtype=complex)
    for i in range(len(x_list)):
        PP_list[i]=[[np.exp((-1j)*phase_shift[i]), 0],[0, np.exp((1j)*phase_shift[i])]]
    
#transfer matrix at ith point
    Tmi=np.zeros((len(x_list), 2, 2), dtype=complex)
    Tmi=PP_list
    for i in range(len(Tmi)-1):
        if is_interface[i]==1:
            B_condition=BC(n_list[i], n_list[i+1], theta[i], theta[i+1])
            Tmi[i]=np.matmul(B_condition, PP_list[i])


    for i in range(len(x_list)-1):
        Ep[i+1], Em[i+1]=np.matmul(Tmi[i], [Ep[i], Em[i]])


       
    Ez=Ep+Em
  
    
    #In case of p polarisation, the H-field is calculated by this procedure
    #The variable in which the H-field is stored is Ez
    
    Ez    
    #returns only the relative field amplitude
    return x_list, n_list, Ez


"""
wl=645.0
theta_i=0.0
d_input=[50, 90, 50, 90, 50, 90, 50, 90, 50, 90]
n_input=[1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5, ]
Ef=0.75
Eb=0.0
x_list, n_list, E, Ep, Em = get_E_field(wl, theta_i, n_input, d_input, Ef, Eb)


I=np.abs(E)**2
plt.plot(x_list, I)


"""
