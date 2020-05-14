#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt 
import pylab
import scipy
import scipy.special as sps
from scipy import integrate


# In[2]:


H0_thin = 95 #*u.pc #the sacale height of the galaxy      LGRB&CCSN
R0_thin = 2.9e3 #*u.pc # the scale radius     LGRB&CCSN

H0_thick = 800#*u.pc #the sacale height of the galaxy        SGRB and Ia
R0_thick = 2.4e3#*u.pc # the scale radius     SGRB and Ia


z = 0#*u.pc #galaxy mid-plane height of the sun
R= 0#*u.pc #the placement of the sun 
R_sun = 0#*u.pc 

#Green Rate
a = 1.09
b = 3.87
h = 95#*u.pc #the scale height -> 
R_green = 8500#*u.pc 
z = 0 #*u.pc #galaxy mid-plane height of the sun


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

Rate_ccsn_l = (0.006) #*(1/u.yr)#Milky Way galactic CCSN rate lower limit
Rate_ccsn = (0.032) #*(1/u.yr)#Milky Way galactic CCSN rate 
Rate_ccsn_u = (0.105) #*(1/u.yr)#Milky Way galactic CCSN rate upper limit

R_sn_l = 0.006 #*(1/u.yr)
R_sn = 0.0139 #*(1/u.yr) #type Ia sn rate 
R_sn_u = 0.028 #*(1/u.yr)

z_max = 0 #*u.pc the average z_max from Bahcal


# In[3]:


def integral_Adams(rhigh, R_sn, H0, R0):
    
    def q_adams(R,z):
        return ((R_sn)/(4*np.pi*H0*R0**2))*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)    
    
    def integrand(b,l,r): 
        R = (R_sun**2+(r*np.cos(b))**2-2*R_sun*r*np.cos(b)*np.cos(l))**(1/2)
        z_sun = (2*z_max/(np.pi))
        z = r*np.sin(b)+z_sun
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


# In[5]:


xd = integral_Adams(10,Rate_ccsn_l,H0_thin,R0_thin)*10**9
x=integral_Adams(10,Rate_ccsn,H0_thin,R0_thin)*10**9
xu = integral_Adams(10,Rate_ccsn_u,H0_thin,R0_thin)*10**9

x1 = (integral_Adams(10,R_sn_u*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn_u*0.5,H0_thin, R0_thin))*10**9
y1 = (integral_Adams(10,R_sn*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn*0.5,H0_thin, R0_thin))*10**9
z1 = (integral_Adams(10,R_sn_l*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn_l*0.5,H0_thin, R0_thin))*10**9

upper = (xu+x1)-(y1+x)
mean = y1+x
lower = (y1+x)-(xd+z1)


print ("+ %s" %(upper))
print ("Trilegal SN rate at GC is %s" %(mean))
print ("- %s" %(lower))


# In[7]:


xd=integral_Adams(Rad_lgrb,R_lgrb_l,H0_thin,R0_thin)*10**9
x=integral_Adams(Rad_lgrb,R_lgrb,H0_thin,R0_thin)*10**9
xu=integral_Adams(Rad_lgrb,R_lgrb_u,H0_thin,R0_thin)*10**9

z1=integral_Adams(Rad_sgrb,R_sgrb_l,H0_thick,R0_thick)*10**9
y1=integral_Adams(Rad_sgrb, R_sgrb, H0_thick,R0_thick)*10**9
x1=integral_Adams(Rad_sgrb,R_sgrb_u, H0_thick,R0_thick)*10**9

upper = (xu+x1)-(y1+x)
mean = y1+x
lower = (y1+x)-(xd+z1)


print ("+ %s" %(upper))
print ("Trilegal GRB rate at GC is %s" %(mean))
print ("- %s" %(lower))

