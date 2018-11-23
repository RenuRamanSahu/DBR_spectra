"""
This program calculates the E_field for a particular wavelength
"""

import numpy as np
import matplotlib.pyplot as plt
import core_lib as cl
import E_lib as el


"""
theta_i=0.
n_med=1.0


n_input, d_list= cl.n_and_d(1.5, 2.5, 50, 90, 3)

metal_coat=[True, 168.26, 8920.34, 30, 'l']

wl=645  
"""







    

def plot_field(wl,pol, theta_i, n_med, n_input, d_list, metal_coat):
    
    if metal_coat[0]==True:
        wl_p=metal_coat[1]
        wl_c=metal_coat[2]
        d_met=metal_coat[3]
        n_met=cl.n_metal(wl, wl_p, wl_c)
        if metal_coat[4]=='r':
            n_input=np.concatenate([n_input, [n_met], [1]])
            d_list=np.concatenate([d_list, [d_met], [10]])
        else:
            n_input=np.concatenate([[n_met], n_input,  [1]])
            d_list=np.concatenate([[d_met],  d_list,  [10]])
    
    else:
        n_input=np.concatenate([n_input, [1]])
        d_list=np.concatenate([d_list, [10]])
        
    int_dlist=np.zeros(len(d_list), dtype=int)    
      
    for i in range(len(d_list)):
        int_dlist[i]=int(d_list[i])

    d_list=int_dlist    
    #print type(d_list[0])
    
    Ef=0.75
    Eb=0
    x_list, n_list, Ex, Ey, Ez = el.get_E_field(wl, pol, theta_i, np.flipud(n_input), np.flipud(d_list), Ef, Eb)
    
    
   
    n_max=max(n_list)
    
    n_0=np.flipud(n_list)
    Ex_0=np.flipud(Ex)
    Ey_0=np.flipud(Ey)
    Ez_0=np.flipud(Ez)


    
    I=np.abs(Ex_0)**2   + np.abs(Ey_0)**2+  np.abs(Ez_0)**2
    I_max=max(I)
    I=n_max*I/I_max
    
    
    wvl_char=str(int(wl))
    field_title="Intensity at "+wvl_char+" nm"
    plt.figure()
    plt.title(field_title)
    
     #Draws lines of interfaces
    xcoords = []
    for i in range(len(n_0)-1):
        if n_0[i]!=n_0[i+1]:
            xcoords=np.concatenate([xcoords, [x_list[i]]])
    
   
    
    for xc in xcoords:
        plt.axvline(x=xc, linestyle=':',  color='r')
    
    plt.plot(x_list, I, '-', color='b')
    plt.plot(x_list, n_0, color='g')
    plt.xlabel(' (nm)')
    plt.legend()
    plt.show()
       
#plot_field(wl, theta_i, n_med, n_input, d_list, metal_coat)    
    
    
    
    
    
    
    
    
    
    
    
    
  
