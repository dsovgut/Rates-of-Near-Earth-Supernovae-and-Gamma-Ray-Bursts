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


#z = 24#*u.pc #galaxy mid-plane height of the sun
R= 8700#*u.pc #the placement of the sun 
R_sun = 8700#*u.pc 

#Green Rate
a = 1.09
b = 3.87
h = 95#*u.pc #the scale height
R_sun = 8700#*u.pc 
#z = 24 #*u.pc #galaxy mid-plane height of the sun


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


def integral2_Adams(rhigh,R_sn, H0, R0, z_sun):
    
    def q_adams(R,z):
        return ((R_sn)/(4*np.pi*H0*R0**2))*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)    
    
    def integrand(b,l,r): 
        R = np.sqrt(R_sun**2+(r*np.cos(b))**2-2*R_sun*r*np.cos(b)*np.cos(l))
        z = r*np.sin(b)+z_sun
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


# # CCSN and LGRB rates with solar heights [z] = 0; 93

# In[2]:


l = integral2_Adams(10,R_ccsn_l, H0_thin, R0_thin, 0)*10**9
a = integral2_Adams(10,R_ccsn, H0_thin, R0_thin, 0)*10**9
u = integral2_Adams(10,R_ccsn_u, H0_thin, R0_thin, 0)*10**9

print('lower value - %s' %(a-l))
print('average value %s' %(a))
print('upper value + %s' %(u-a))


# In[3]:


l = integral2_Adams(10,R_ccsn_l, H0_thin, R0_thin, 93)*10**9
a = integral2_Adams(10,R_ccsn, H0_thin, R0_thin, 93)*10**9
u = integral2_Adams(10,R_ccsn_u, H0_thin, R0_thin, 93)*10**9

print('lower value - %s' %(a-l))
print('average value %s' %(a))
print('upper value + %s' %(u-a))


# In[4]:


l = integral2_Adams(Rad_lgrb,R_lgrb_l, H0_thin, R0_thin, 0)*10**9
a = integral2_Adams(Rad_lgrb,R_lgrb, H0_thin, R0_thin, 0)*10**9
u = integral2_Adams(Rad_lgrb,R_lgrb_u, H0_thin, R0_thin, 0)*10**9

print('lower value - %s' %(a-l))
print('average value %s' %(a))
print('upper value + %s' %(u-a))


# In[5]:


l = integral2_Adams(Rad_lgrb,R_lgrb_l, H0_thin, R0_thin, 93)*10**9
a = integral2_Adams(Rad_lgrb,R_lgrb, H0_thin, R0_thin, 93)*10**9
u = integral2_Adams(Rad_lgrb,R_lgrb_u, H0_thin, R0_thin, 93)*10**9

print('lower value - %s' %(a-l))
print('average value %s' %(a))
print('upper value + %s' %(u-a))


# In[11]:


z_sun = [0, 25, 75, 100, 150]  
r1 = np.logspace(0, 3.0, num=30, base = 10)
r = r1.tolist()
print(r)


# In[12]:


CCSN_z0 = []
CCSN_z25 = []
CCSN_z75 = []
CCSN_z100 = []
CCSN_z150 = []

distance = []


# # Generating Dataset for CCSN Solar Heights

# In[13]:


for i in r:
    distance.append(i)
    CCSN_z0.append(integral2_Adams(i, R_ccsn, H0_thin, R0_thin, z_sun[0])) #R_ccsn - global rate
    CCSN_z25.append(integral2_Adams(i, R_ccsn, H0_thin, R0_thin, z_sun[1]))
    CCSN_z75.append(integral2_Adams(i, R_ccsn, H0_thin, R0_thin, z_sun[2]))
    CCSN_z100.append(integral2_Adams(i, R_ccsn, H0_thin, R0_thin, z_sun[3]))
    CCSN_z150.append(integral2_Adams(i, R_ccsn, H0_thin, R0_thin, z_sun[4]))
    
data = {'distance':distance,'CCSN_z0':CCSN_z0, 'CCSN_z25':CCSN_z25,'CCSN_z75':CCSN_z75, 'CCSN_z100':CCSN_z100, 'CCSN_z150':CCSN_z150 } 
df = pd.DataFrame(data) 
#print(df[['distance','CCSN_z0','CCSN_z25', 'CCSN_z75', 'CCSN_z100', 'CCSN_z150']]) 
df.to_csv('CCSN_zheights.csv')


# # Generating Dataset for Type Ia Supernvoae Solar Heights

# In[14]:


sn_z0 = []
sn_z25 = []
sn_z75 = []
sn_z100 = []
sn_z150 = []

distance = []


for i in r:
    distance.append(i)
    sn_z0.append((integral2_Adams(i, R_sn*0.5, H0_thin, R0_thin, z_sun[0])+integral2_Adams(i, R_sn*0.5, H0_thick, R0_thick, z_sun[0])))
    sn_z25.append((integral2_Adams(i, R_sn*0.5, H0_thin, R0_thin, z_sun[1])+integral2_Adams(i, R_sn*0.5, H0_thick, R0_thick, z_sun[1])))
    sn_z75.append((integral2_Adams(i, R_sn*0.5, H0_thin, R0_thin, z_sun[2])+integral2_Adams(i, R_sn*0.5, H0_thick, R0_thick, z_sun[2])))
    sn_z100.append((integral2_Adams(i, R_sn*0.5, H0_thin, R0_thin, z_sun[3])+integral2_Adams(i, R_sn*0.5, H0_thick, R0_thick, z_sun[3])))
    sn_z150.append((integral2_Adams(i, R_sn*0.5, H0_thin, R0_thin, z_sun[4])+integral2_Adams(i, R_sn*0.5, H0_thick, R0_thick, z_sun[4])))


    
data = {'distance':distance,'sn_z0':sn_z0, 'sn_z25':sn_z25,'sn_z75':sn_z75, 'sn_z100':sn_z100, 'sn_z150':sn_z150 } 
df = pd.DataFrame(data) 
print(df[['distance','sn_z0','sn_z25', 'sn_z75', 'sn_z100', 'sn_z150',]]) 
df.to_csv('SN_zheights.csv')


# In[ ]:




