#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 19:24:42 2018

@author: dilu
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 16:21:00 2018

@author: dilu
"""

import numpy as np
#import matplotlib.pyplot as plt
import sys

#e stands for epsilon: The difference between 1 and the 
#next posiible float point number in python
e=sys.float_info.epsilon


"""
metal_coat('Bool', wl_p, wl_c, d_metal, metal_pos)
Bool = This variable tells if metal is coated or not
wl_p = This is the plasma wavelength of metal
wl_c= This is the collision wavelength of metal
d_metal = This is the thickness of metal coat
metal_pos='l' or 'r'.....if 'l' the metal is on the incident surface
                         if 'r' the metal is on the opposite side
"""


def n_and_d(n1, n2, d1, d2, N):
    n=[]
    d=[]
    for i in range(N):
        n=np.concatenate([[n1,n2], n])
        d=np.concatenate([[d1, d2], d])
   # n=np.concatenate([n,[1]] )
   # d=np.concatenate([d,[10]])
    return n, d
    
    

   
 


def get_theta(theta_i, n_med, n_input):
    n0_sin_theta0 =n_med*np.sin(theta_i)
    l=len(n_input)
    theta_list=np.zeros(l, dtype=complex)
    
    for i in range(l):
        n=n_input[i]
        re_n=np.real(n)
        im_n=np.imag(n)
        theta=np.arcsin(n0_sin_theta0/n)
        ncos_theta=np.sqrt(n**2 - (n0_sin_theta0)**2)
        chk_absrp=re_n*im_n
        assert (chk_absrp>=-100*e), "Ambiguous case: Can't handle active material"
        if  np.real(ncos_theta) <-e*100 :
            theta=np.pi-theta
        
        if chk_absrp<=np.abs(100*e):
            assert (re_n>-100*2), "Ambiguous Case"
            if re_n>n0_sin_theta0:
                if np.imag(ncos_theta)<-100*e:
                    theta=np.pi-theta
            else:
                if np.real(ncos_theta)<-100*e:
                    theta=np.pi-theta  

        theta_list[i]=theta
     
    return theta_list






def r12_s(n1, n2, th1, th2):
    num=n1*np.cos(th1)-n2*np.cos(th2)
    denom=n1*np.cos(th1)+n2*np.cos(th2)
    return num/denom

def t12_s(n1, n2, th1, th2):
    num=2*n1*np.cos(th1)
    denom=n1*np.cos(th1)+n2*np.cos(th2)
    return num/denom


def r12_p(n1, n2, th1, th2):
    num=-(n2*np.cos(th1)-n1*np.cos(th2))
    denom=n2*np.cos(th1)+n1*np.cos(th2)
    return num/denom

def t12_p(n1, n2, th1, th2):
    num=2*n1*np.cos(th1)
    denom=n2*np.cos(th1)+n1*np.cos(th2)
    return num/denom

def n_metal(wl, wl_p, wl_c):
    #o stands for omega
    op=1.0/wl_p
    o=1.0/wl
    oc=1.0/wl_c
    num=op**2
    denom=o**2 + (1j)*o*oc
    nsq=1-num/denom
    n=np.sqrt(nsq)
    return n

 
def M(wl, pol,  theta_i, n_med, n_input, d_list, metal_coat):
    
    if metal_coat[0]==True:
        wl_p=metal_coat[1]
        wl_c=metal_coat[2]
        d_met=metal_coat[3]
        n_met=n_metal(wl, wl_p, wl_c)
        if metal_coat[4]=='r':
            n_input=np.concatenate([n_input, [n_met], [1]])
            d_list=np.concatenate([d_list, [d_met], [10]])
        else:
            n_input=np.concatenate([[n_met],[n_input[1]], n_input,  [1]])
            d_list=np.concatenate([[d_met], [d_list[1]-d_met], d_list,  [10]])
    
    else:
        n_input=np.concatenate([n_input, [1]])
        d_list=np.concatenate([d_list, [10]])
    
    ln=len(n_input)
    
   
    k_list=np.zeros(ln, dtype=complex)
    kz_list=np.zeros(ln, dtype=complex)
    del_list=np.zeros(ln, dtype=complex)
    r_n_np1=np.zeros(ln-1, dtype=complex)
    t_n_np1=np.zeros(ln-1, dtype=complex)

    k_list=(2*np.pi/wl)*np.copy(n_input)
    theta_list=get_theta(theta_i, n_med, n_input)
    kz_list=k_list*np.cos(np.array(theta_list))
    del_list=kz_list*d_list
    
    if pol=='s':
         for i in range(ln-1):
             r_n_np1[i]=r12_s(n_input[i], n_input[i+1], theta_list[i], theta_list[i+1])
             t_n_np1[i]=t12_s(n_input[i], n_input[i+1], theta_list[i], theta_list[i+1])
    else:
         for i in range(ln-1):
             r_n_np1[i]=r12_p(n_input[i], n_input[i+1], theta_list[i], theta_list[i+1])
             t_n_np1[i]=t12_p(n_input[i], n_input[i+1], theta_list[i], theta_list[i+1])
    
         



    
    d=np.copy(del_list)
    r=np.copy(r_n_np1)
    t=np.copy(t_n_np1)
    
    l=ln-1
    M=np.zeros((l, 2, 2), dtype=complex)

    TM=[[1, 0],[0, 1]]
    for i in range(l):
        dn=d[i] 
        rn=r[i]
        tn=t[i]
        pp=[[np.exp((-1j)*dn), 0],[0, np.exp((1j)*dn)]]
        rt=[[1/tn, rn/tn],[rn/tn, 1/tn]]
        M[i]=np.matmul(pp, rt)
        TM=np.matmul(TM, M[i] )
    
    r01=r12_s(n_med, n_input[0], theta_i, theta_list[0])
    t01=t12_s(n_med, n_input[0], theta_i, theta_list[0])
    M0=[[1/t01, r01/t01],[r01/t01, 1/t01]]
    TM=np.matmul(M0, TM)
    return  TM


def phase_r(r):
    re= np.real(r)
    im=np.imag(r)
    phase=np.arctan(im/re)
    return phase

def plot_RT_vs_wavelength(wl,pol,  theta_i, n_med, n_input, d_list, metal_coat):
    wvl=np.copy(wl)
    l=len(wl)
    r=np.zeros(l, dtype=complex)
    t=np.zeros(l, dtype=complex)
    ph_r=np.zeros(l, dtype=complex)

    R=np.zeros(l)
    T=np.zeros(l)
    
    for i in range(l):
        TM=M(wvl[i],pol, theta_i, n_med, n_input, d_list, metal_coat)
        ri=TM[1][0]/TM[0][0]
        ti=1/TM[0][0]
        r[i]=ri
        t[i]=ti
        ph_r[i]=phase_r(ri)
        R[i]=np.abs(ri)**2
        T[i]=np.abs(ti)**2
    return wvl, ph_r, R, T





    
