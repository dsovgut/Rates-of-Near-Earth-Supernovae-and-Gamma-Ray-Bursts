#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
R_sun = 8700#*u.pc #we used 8700 bc what was in adams paper not what 8500 what says in the blog

#Green Rate
a = 1.09
b = 3.87
h = 95#*u.pc #the scale height
R_sun = 8700#*u.pc #we used 8700 bc what was in adams paper not what 8500 what says in the blog
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

def integral2_Adams(rhigh, H0, R0):
    
    def q_adams(R,z):
        return ((1)/(4*np.pi*H0*R0**2))*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)    
    
    def integrand(b,l,r): 
        R = (R_sun**2+(r*np.cos(b))**2-2*R_sun*r*np.cos(b)*np.cos(l))**(1/2)
        z = r*np.sin(b)+(2*z_max/(np.pi))
        q=q_adams(R,z)
        return 2*q*r**2*np.cos(b) 
    
    rlow = 0
    value1, error1 = scipy.integrate.tplquad(integrand, rlow,rhigh, lambda r: 0, lambda r: np.pi/2, #l
                                lambda r, l: -np.pi/36, lambda r, l: np.pi/36) #b
    
    value2, error2 = scipy.integrate.tplquad(integrand, rlow,rhigh, lambda r: np.pi/2, lambda r: np.pi,
                                lambda r, l: -np.pi/36, lambda r, l: np.pi/36)
    
    value3, error3 = scipy.integrate.tplquad(integrand, rlow,rhigh, lambda r: 0, lambda r: np.pi,
                                lambda r, l: np.pi/36, lambda r, l: np.pi/2)
    
    value4, error4 = scipy.integrate.tplquad(integrand, rlow,rhigh, lambda r: 0, lambda r: np.pi,
                                lambda r, l: -np.pi/2, lambda r, l: -np.pi/36) 
    
    answer = value1+value2+value3+value4
    
    return (answer)


# In[2]:


r1 = np.logspace(0, 4.4759, num=100, base = 10)
r = r1.tolist()


# In[3]:


R_al=7000
H_al=800
R_fe=3500
H_fe=300


# In[4]:


Al = []
Fe = []
distance = []

for i in r:
    distance.append(i)
    Al.append(integral2_Adams(i, H_al,R_al))
    Fe.append(integral2_Adams(i, H_fe,R_fe))
    
data = {'distance':distance,'26Al':Al, '60Fe':Fe} 
df = pd.DataFrame(data) 
#print(df[['distance','26Al','60Fe']]) 
df.to_csv('scale_rates_aluminum_ferum.csv')


lgrb_al = []
sn_al = []
ccsn_al = []

lgrb_fe = []
sn_fe = []
ccsn_fe = []

for i in Al:
    lgrb_al.append(i*R_lgrb)
    ccsn_al.append(i*R_ccsn)
    sn_al.append(i*0.5*R_sn)
    

for i in Fe:
    lgrb_fe.append(i*R_lgrb)
    ccsn_fe.append(i*R_ccsn)
    sn_fe.append(i*0.5*R_sn)
    
data = {'distance':distance, 'lgrb_al':lgrb_al, 'ccsn_al':ccsn_al, 'sn_al':sn_al, 'lgrb_fe':lgrb_fe, 'ccsn_fe':ccsn_fe,'sn_fe':sn_fe} 
df = pd.DataFrame(data) 
#print(df[['distance','Thin','Thick']]) 
df.to_csv('rates-Al-Fe.csv')


# In[ ]:




