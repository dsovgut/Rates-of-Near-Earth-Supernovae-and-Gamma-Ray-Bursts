#!/usr/bin/env python
# coding: utf-8

# In[18]:


import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt 
import pylab
import scipy
import scipy.special as sps
from scipy import integrate


# In[19]:


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

Rate_ccsn_l = (0.006) #*(1/u.yr)#Milky Way galactic CCSN rate lower limit
Rate_ccsn = (0.032) #*(1/u.yr)#Milky Way galactic CCSN rate 
Rate_ccsn_u = (0.105) #*(1/u.yr)#Milky Way galactic CCSN rate upper limit

R_sn_l = 0.006 #*(1/u.yr)
R_sn = 0.0139 #*(1/u.yr) #type Ia sn rate 
R_sn_u = 0.028 #*(1/u.yr)

z_max = 71 #*u.pc the average z_max from Bahcal


# In[20]:


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


# ______________________________________________________________________________________________________________________

# In[40]:


CCSN = []
r = np.linspace(1,500e6,1000)
c1= integral_Adams(10,Rate_ccsn,H0_thin,R0_thin)
CCSN = c1*r

CCSN_d = []
cd = integral_Adams(10,Rate_ccsn_l,H0_thin,R0_thin)
CCSN_d = cd*r

CCSN_up = []
cu = integral_Adams(10,Rate_ccsn_u,H0_thin,R0_thin)
CCSN_up = cu*r


# In[22]:


Ia = []
i1= integral_Adams(10,R_sn*0.5,H0_thick,R0_thick)+integral_Adams(10,R_sn*0.5,H0_thin,R0_thin)
Ia = i1*r

Ia_d = []
i2 = integral_Adams(10,R_sn_l*0.5,H0_thick,R0_thick)+integral_Adams(10,R_sn_l*0.5,H0_thin,R0_thin)
Ia_d = i2*r

Ia_up = []
i3 = integral_Adams(10,R_sn_u*0.5, H0_thick,R0_thick)+integral_Adams(10,R_sn_u*0.5, H0_thin,R0_thin)
Ia_up = i3*r


# In[46]:


print((1/(i1*10**9))) #actual mean time
print((1/(i2*10**9))) #lower mean time
print((1/(i3*10**9))) #upper mean time


# In[48]:


9.759065418708254-4.844678904287311


# In[23]:


LGRB = []
l1= integral_Adams(Rad_lgrb,R_lgrb,H0_thin,R0_thin)
LGRB = l1*r

LGRB_d = []
l2 = integral_Adams(Rad_lgrb,R_lgrb_l,H0_thin,R0_thin)
LGRB_d = l2*r

LGRB_up = []
l3 = integral_Adams(Rad_lgrb,R_lgrb_u,H0_thin,R0_thin)
LGRB_up = l3*r


# In[24]:


SGRB = []
s1= integral_Adams(Rad_sgrb,R_sgrb,H0_thick,R0_thick)
SGRB = s1*r

SGRB_d = []
s2 = integral_Adams(Rad_sgrb,R_sgrb_l,H0_thick,R0_thick)
SGRB_d = s2*r

SGRB_up = []
s3 = integral_Adams(Rad_sgrb,R_sgrb_u, H0_thick,R0_thick)
SGRB_up = s3*r


# In[32]:


y1 = [0,1,2,3,4,5] 
x1 = [66, 201.3, 252, 375,450,550]

a = np.linspace(0,500,1000)

plt.plot(a, CCSN, color = 'black', label = 'Expected CCSN')
plt.plot(a, CCSN_up, color = 'black',alpha=0.5)
plt.plot(a, CCSN_d, color = 'black',alpha=0.5)
plt.fill_between(a, CCSN_up,CCSN_d,color='gray',alpha=.2, label='CCSN Uncertainty')

plt.plot(a, LGRB, color = 'c', label = 'Expected LGRB')
plt.plot(a, LGRB_up, color = 'c',alpha=0.5)
plt.plot(a, LGRB_d, color = 'c',alpha=0.5)
plt.fill_between(a, LGRB_up,LGRB_d,color='c',alpha=.2, label='LGRB Uncertainty')

plt.xlabel ('Million Years Ago t [Myr]')
plt.ylabel (r'Cumilative Number of Extinctions $\Gamma(r)$ [$\frac{events}{Myr}$]')

plt.xlim(0,500)
plt.ylim(0,6)


plt.step(x1, y1, color = 'red',alpha=1, label='Mass Exctinction Events')
plt.legend()
plt.savefig('Exctinction_LGRB+CCSN.png', dpi=300)
plt.show()


# In[42]:


rate = c1+i1+l1+s1
rate_d = cd+i2+l2+s2
rate_u = cu+i3+l3+s3

R=[]
R=rate*r

R_d=[]
R_d =rate_d*r

R_u=[]
R_u=rate_u*r


# In[43]:


y = [0,1,2,3,4,5] 
x = [66, 201.3, 252, 375,450,550]

a = np.linspace(0,500,1000)

plt.plot(a, R, color = 'b', label = 'Expected Cosmic Explosions')
plt.plot(a, R_d, color = 'cornflowerblue',alpha=0.2)
plt.plot(a, R_u, color = 'cornflowerblue',alpha=0.2)
plt.fill_between(a, R_u,R_d,color='cornflowerblue',alpha=.2, label='Cosmic Explosions Uncertainty')


plt.xlabel ('Million Years Ago t [Myr]')
plt.ylabel (r'Cumilative Number of Extinctions $\Gamma(r)$ [$\frac{events}{Myr}$]')
plt.title("Extinction Events vs Expected Cosmic Explosion", fontsize = 13)

plt.xlim(0,500)
plt.ylim(0,6)


plt.step(x, y, color = 'red',alpha=1, label='Mass Exctinction Events')
plt.legend()
plt.savefig('Exctinction_combined.png', dpi=300)
pylab.show()


# In[11]:


Ext_rate = []
Ext_rate = r*(1e-08)


# In[12]:


y = [0,1,2,3,4,5] 
x = [66, 201.3, 252, 375,450,550]

a = np.linspace(0,500,1000)

plt.plot(a, R, color = 'b', label = 'Expected Cosmic Explosions')
plt.plot(a, R_d, color = 'cornflowerblue',alpha=0.2)
plt.plot(a, R_u, color = 'cornflowerblue',alpha=0.2)
plt.fill_between(a, R_u,R_d,color='cornflowerblue',alpha=.2, label='error')


plt.xlabel ('Million Years Ago')
plt.ylabel ('Cumilative Number of Extinctions')
plt.title("Extinction Events vs Expected Cosmic Explosion ", fontsize = 13)

plt.xlim(0,500)
plt.ylim(0,6)


plt.step(x, y, color = 'red',alpha=1, label='Mass Exctinction Events')
plt.legend()
#plt.savefig('Exctinction.png')
pylab.show()

