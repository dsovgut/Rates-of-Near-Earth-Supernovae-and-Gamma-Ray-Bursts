#!/usr/bin/env python
# coding: utf-8

# In[8]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# In[9]:


distance, LGRB_l, LGRB,LGRB_u,CCSN_l, CCSN,CCSN_u,SGRB_l, SGRB,SGRB_u,TypeIa_l,TypeIa,TypeIa_u = np.loadtxt(fname = 'final_rates2.csv', dtype = 'float', unpack= 'True', usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13), skiprows = 1, delimiter=',')

#splitting TypeIa supernovae into thick and thin disk
thin, thick = np.loadtxt(fname = 'scale_rates2.csv', dtype = 'float', unpack= 'True', usecols = (2,3), skiprows = 1, delimiter=',')
R_sn_l = 0.006 #*(1/u.yr)
R_sn = 0.0139 #*(1/u.yr) #type Ia sn rate 
R_sn_u = 0.028 #*(1/u.yr)
sn_l = []
sn = []
sn_u= []
for thin1,thick1 in zip(thin, thick):
    sn_l.append((thin1*R_sn_l*0.5)+(thick1*R_sn_l*0.5))
    sn.append((thin1*R_sn*0.5)+(thick1*R_sn*0.5))
    sn_u.append((thin1*R_sn_u*0.5)+(thick1*R_sn_u*0.5))


# In[10]:


plt.plot(distance, CCSN, color = 'blue', label = 'CCSN') #red and blue, m, green 
plt.plot(distance, sn, color = 'red', label = 'Type Ia') #red and blue, m, green 
plt.plot(distance, SGRB, color = 'green', label = 'Short GRB')
plt.plot(distance, LGRB, color = 'm', label = 'Long GRB')

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, sn_u,sn_l,color='red',alpha=.09, linestyle = ':')
plt.fill_between(distance, SGRB_l,SGRB_u,color='green',alpha=.09, linestyle = '--')
plt.fill_between(distance, LGRB_l,LGRB_u,color='m',alpha=.09, linestyle = ':')

plt.ylim(1e-10,10e-1)
plt.xlim(10,29915.75721530859)

plt.title("Rate vs Distance (Trilegal Model)", fontsize = 15)

plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance (Trilegal Model)", fontsize = 15)

plt.legend()
plt.savefig('Rate_vs_Distance_witherrors.png', dpi=300)

plt.show()


# # Comparison of Approximations

# In[4]:


import numpy as np
from astropy.table import Table, Column
from astropy import units as u
import pandas as pd 

import scipy
import scipy.special as sps
from scipy import integrate

H0_thin = 95 #*u.pc #the sacale height of the galaxy      LGRB&CCSN
R0_thin = 2.9e3 #*u.pc # the scale radius     LGRB&CCSN

H0_thick = 800#*u.pc #the sacale height of the galaxy        SGRB and Ia
R0_thick = 2.4e3#*u.pc # the scale radius     SGRB and Ia


z = 24#*u.pc #galaxy mid-plane height of the sun
R= 8700#*u.pc #the placement of the sun 
R_sun = 8700#*u.pc 

#Green Rate
a = 1.09
b = 3.87
h = 95#*u.pc #the scale height
R_sun = 8700#*u.pc 
z = 24 #*u.pc #galaxy mid-plane height of the sun


#kill radius
Rad_sn=10 #*u.pc
Rad_sgrb=91.42 #*u.pc
Rad_lgrb=2044 #*u.pc

R_sgrb_l = 1.51e-07 #* (1/u.yr)
R_sgrb = 3.14e-06 #* (1/u.yr)
R_sgrb_u=0.000101 #* (1/u.yr)

R_lgrb_l = 1.59e-08 #* (1/u.yr)
R_lgrb = 1.22e-07 #* (1/u.yr)
R_lgrb_u=6.95e-07 #* (1/u.yr)

R_ccsn_l = (0.006) #*(1/u.yr)#Milky Way galactic CCSN rate lower limit
R_ccsn = (0.032) #*(1/u.yr)#Milky Way galactic CCSN rate 
R_ccsn_u = (0.105) #*(1/u.yr)#Milky Way galactic CCSN rate upper limit

R_sn_l = 0.006 #*(1/u.yr)
R_sn = 0.0139 #*(1/u.yr) #type Ia sn rate 
R_sn_u = 0.028 #*(1/u.yr)

z_max = 71 #*u.pc the average z_max from Bahcal


# In[5]:


def Adams_uniform(r, R_sn, H0,R0): 
    q0 = (R_sn)/(4*np.pi*H0*R0**2) #normalization 
    q_adams = q0*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)
    return q_adams*((4/3)*np.pi*(r**3))


# In[6]:


def Adams_expheight(r, R_sn, H0, R0): 
    q0 = (R_sn)/(4*np.pi*H0*R0**2) #normalization 
    q_adams = q0*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)
    factor_exp=r**2*H0+2*r*H0**2*np.exp(-r/H0)-2*H0**3*(1-np.exp(-r/H0))
    return 2*np.pi*q_adams*factor_exp


# In[7]:


def integral_Adams(rhigh, R_sn, H0, R0):
    
    def q_adams(R,z):
        return ((R_sn)/(4*np.pi*H0*R0**2))*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)    
    
    def integrand(b,l,r): 
        R = (R_sun**2+(r*np.cos(b))**2-2*R_sun*r*np.cos(b)*np.cos(l))**(1/2)
        z = r*np.sin(b)+(2*z_max/(np.pi))
        q=q_adams(R,z)
        return 2*q*r**2*np.cos(b) 
        
    value1, error1 = scipy.integrate.tplquad(integrand, 0,rhigh, lambda r: 0, lambda r: np.pi/2, #l
                                lambda r, l: -np.pi/36, lambda r, l: np.pi/36) #b
    
    value2, error2 = scipy.integrate.tplquad(integrand, 0,rhigh, lambda r: np.pi/2, lambda r: np.pi,
                                lambda r, l: -np.pi/36, lambda r, l: np.pi/36)
    
    value3, error3 = scipy.integrate.tplquad(integrand, 0,rhigh, lambda r: 0, lambda r: np.pi,
                                lambda r, l: np.pi/36, lambda r, l: np.pi/2)
    
    value4, error4 = scipy.integrate.tplquad(integrand, 0,rhigh, lambda r: 0, lambda r: np.pi,
                                lambda r, l: -np.pi/2, lambda r, l: -np.pi/36) 
    
    answer = value1+value2+value3+value4
    
    return (answer)


# In[13]:


integral_Adams(10, R_ccsn, H0_thick, R0_thick)


# In[14]:


Adams_expheight(10, R_ccsn, H0_thick, R0_thick)


# In[8]:


uniform = np.empty(0, dtype = 'float64')
expheight = np.empty(0, dtype = 'float64')

uniform_u = np.empty(0, dtype = 'float64')
expheight_u = np.empty(0, dtype = 'float64')

uniform_l = np.empty(0, dtype = 'float64')
expheight_l = np.empty(0, dtype = 'float64')


for i in distance:
    #CCSN actual value
    
    result = Adams_uniform(i, R_ccsn, H0_thin, R0_thin) 
    uniform = np.append(uniform, [result], axis=0)
    
    result2 = Adams_expheight(i, R_ccsn, H0_thin, R0_thin) 
    expheight = np.append(expheight, [result2], axis=0)
    
    #CCSN upper limit
    
    result3 = Adams_uniform(i, R_ccsn_u, H0_thin, R0_thin) 
    uniform_u = np.append(uniform_u, [result3], axis=0)
    
    result4 = Adams_expheight(i, R_ccsn_u, H0_thin, R0_thin) 
    expheight_u = np.append(expheight_u, [result4], axis=0)

    #CCSN lower limit

    result5 = Adams_uniform(i, R_ccsn_l, H0_thin, R0_thin) 
    uniform_l = np.append(uniform_l, [result5], axis=0)
    
    result6 = Adams_expheight(i, R_ccsn_l, H0_thin, R0_thin) 
    expheight_l = np.append(expheight_l, [result6], axis=0)


    


# In[11]:


plt.plot(distance, uniform, color = 'green', label = 'Local Approximation')
plt.plot(distance, expheight, color = 'red', label = 'Scale Height Approximation')
plt.plot(distance, CCSN, color = 'blue', label = 'Full Integration') #red and blue, m, green 

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, uniform_u,uniform_l,color='green',alpha=.09, linestyle = ':')
plt.fill_between(distance, expheight_u,expheight_l,color='red',alpha=.09, linestyle = '--')


plt.ylim(10**(-8),10**(2))

plt.xlim(50,30000)



plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [events/yr]')
plt.title("Comparing Appriximations for the CCSN", fontsize = 15)

plt.xscale('log')
plt.yscale('log')


plt.legend()
plt.savefig('CCSN_approximation_comparison_full_range_with_errors2.png', dpi=300)
plt.show()


# In[ ]:




