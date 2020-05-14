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


def integral3_Adams(rhigh, H0, R0):
    
    def q_adams(R,z):
        return ((1)/(4*np.pi*H0*R0**2))*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)    
    
    def integrand(b,l,r): 
        R = np.sqrt(R_sun**2+(r*np.cos(b))**2-2*R_sun*r*np.cos(b)*np.cos(l))
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


# In[3]:


#integral Adams is the integral for the TRILEGAL model. 
def integral_Adams(R_sn, rhigh, H0, R0):
    
    def q_adams(R,z):
        return ((R_sn)/(4*np.pi*H0*R0**2))*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)    
    
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


# # Dataset for Table 6 in the paper

# In[4]:


r_t = [10, 20, 30, 40, 50, 100, 125, 150, 200, 500, 1000, 1500, 2000, 2500, 3000, 4000, 10000]


# In[5]:


Thin = []
Thick = []
distance = []

for i in r_t:
    distance.append(i)
    Thin.append(integral2_Adams(i, H0_thin,R0_thin))
    Thick.append(integral2_Adams(i,H0_thick,R0_thick))
    
data = {'distance':distance,'Thin':Thin, 'Thick':Thick} 
df = pd.DataFrame(data) 
print(df[['distance','Thin','Thick']]) 
df.to_csv('scale_table6.csv')


# # Exporting Results for the Rate vs Distance Plot

# In[82]:


r1 = np.logspace(0, 4.4759, num=100, base = 10)
r = r1.tolist()


# In[78]:


Thin = []
Thick = []
distance = []

for i in r:
    distance.append(i)
    Thin.append(integral2_Adams(i, H0_thin,R0_thin))
    Thick.append(integral2_Adams(i,H0_thick,R0_thick))
    
data = {'distance':distance,'Thin':Thin, 'Thick':Thick} 
df = pd.DataFrame(data) 
print(df[['distance','Thin','Thick']]) 
df.to_csv('scale_rates2.csv')


# In[3]:


#thin

R_lgrb_l = 1.59e-08 #* (1/u.yr)
R_lgrb = 1.22e-07 #* (1/u.yr)
R_lgrb_u=6.95e-07 #* (1/u.yr)

R_ccsn_l = (0.006) #*(1/u.yr)#Milky Way galactic CCSN rate lower limit
R_ccsn = (0.032) #*(1/u.yr)#Milky Way galactic CCSN rate 
R_ccsn_u = (0.105) #*(1/u.yr)#Milky Way galactic CCSN rate upper limit

lgrb_l=[]
lgrb=[]
lgrb_u =[]
ccsn_l=[]
ccsn=[]
ccsn_u=[]
sn_l = []
sn = []
sn_u= []


for i in Thin:
    lgrb_l.append(i*R_lgrb_l)
    lgrb.append(i*R_lgrb)
    lgrb_u.append(i*R_lgrb_u)
    ccsn_l.append(i*R_ccsn_l)
    ccsn.append(i*R_ccsn)
    ccsn_u.append(i*R_ccsn_u)
#thick

R_sgrb_l = 1.51e-07 #* (1/u.yr)
R_sgrb = 3.14e-06 #* (1/u.yr)
R_sgrb_u=0.000101 #* (1/u.yr)

R_sn_l = 0.006 #*(1/u.yr)
R_sn = 0.0139 #*(1/u.yr) #type Ia sn rate 
R_sn_u = 0.028 #*(1/u.yr)

sgrb_l = []
sgrb = []
sgrb_u = []

for i in Thick:
    sgrb_l.append(i*R_sgrb_l)
    sgrb.append(i*R_sgrb)
    sgrb_u.append(i*R_sgrb_u)
    sn_l.append(i*R_sn_l)
    sn.append(i*R_sn)
    sn_u.append(i*R_sn_u)



# In[86]:


data = {'distance':distance,'lgrb_l':lgrb_l, 'lgrb':lgrb, 'lgrb_u':lgrb_u , 'ccsn_l':ccsn_l, 'ccsn':ccsn, 'ccsn_u':ccsn_u, 'sgrb_l':sgrb_l, 'sgrb':sgrb, 'sgrb_u':sgrb_u, 'sn_l':sn_l, 'sn':sn , 'sn_u':sn_u} 
df = pd.DataFrame(data) 
#print(df[['distance','Thin','Thick']]) 
df.to_csv('final_rates2.csv')


# In[ ]:




