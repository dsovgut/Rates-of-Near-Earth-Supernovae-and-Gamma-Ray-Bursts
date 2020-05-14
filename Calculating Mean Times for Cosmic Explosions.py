#!/usr/bin/env python
# coding: utf-8

# In[15]:


import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt 
import pylab
import scipy
import scipy.special as sps
from scipy import integrate


# In[16]:


F_kill = 100e3 * u.J/u.m**2
E_Ia = 10**39 * u.J
E_sgrb = 10**43 * u.J
E_lgrb = 5*10**45 * u.J

def r_kill(E):
    return np.sqrt(E/(F_kill*4*np.pi))


# # Calculating Kill Distance from the Fluence

# In[17]:


print(('kill distance lgrb %s') %((r_kill(E_lgrb)).to(u.parsec)))
print(('kill distance sgrb %s') %((r_kill(E_sgrb)).to(u.parsec)))
print(('kill distance TypeIa %s') %((r_kill(E_Ia)).to(u.parsec)))


# In[5]:


H0_thin = 95 #*u.pc #the sacale height of the galaxy      LGRB&CCSN
R0_thin = 2.9e3 #*u.pc # the scale radius     LGRB&CCSN

H0_thick = 800#*u.pc #the sacale height of the galaxy        SGRB and Ia
R0_thick = 2.4e3#*u.pc # the scale radius     SGRB and Ia

H0_al = 800
R0_al = 7000
H0_fe = 300
R0_fe = 3500


z = 24#*u.pc #galaxy mid-plane height of the sun
R= 8700#*u.pc #the placement of the sun 
R_sun = 8700#*u.pc 

#Green Rate
a = 1.09
b = 3.87
h = 95#*u.pc #the scale height -> 
R_green = 8500#*u.pc 
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


# In[6]:


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


# In[7]:


xu = (1/(integral_Adams(1000,Rate_ccsn_l,H0_thin,R0_thin)*10**9))
x=(1/(integral_Adams(1000,Rate_ccsn,H0_thin,R0_thin)*10**9))
xd = (1/(integral_Adams(1000,Rate_ccsn_u,H0_thin,R0_thin)*10**9))

print ("+ %s" %(xu-x))
print ("Adams CCSN is %s" %(x))
print ("- %s" %(x-xd))


# In[24]:


xu = (1/(integral_Adams2(10,Rate_ccsn_l,H0_thin,R0_thin, 0,0)*10**9))
x=(1/(integral_Adams2(10,Rate_ccsn,H0_thin,R0_thin,0,0)*10**9))
xd = (1/(integral_Adams2(10,Rate_ccsn_u,H0_thin,R0_thin,0,0)*10**9))

print ("+ %s" %(xu-x))
print ("Adams CCSN is %s" %(x))
print ("- %s" %(x-xd))


# In[15]:


def integral_Green(r, R_sn,h,R0):
    
    def q_green(R,z):
        R0 = R_sun/b
        f1=(R_sn/(4*np.pi*sps.gamma(a+2)*h*(R0**2)))
        f2=((R_sun/R0)**a)
        f3=(np.exp(-R_sun/R0))
        f4=(np.exp(-np.abs(z)/h))
        return f1*f2*f3*f4
    
    #f = lambda b, l, r: q_adams(R,z)*r**2*np.cos(b) 
    
    def integrand(l,b,r): 
        R = (R_sun**2+(r*np.cos(b))**2-2*R_sun*np.cos(b)*np.cos(l*r))**(1/2)
        z = r*np.sin(b)+((2*z_max)/(np.pi))
        q=q_green(R,z)
        return q*r**2*np.cos(b) 

    value, error = scipy.integrate.tplquad(integrand, 
                                0,r, # r limits
                                lambda r: -np.pi/2, lambda r: np.pi/2, #b limits
                                lambda r, b: 0, lambda r, b: 2*np.pi)  #l limits   
    return (value)


# # 26Al Model Mean Times Calculations

# In[9]:


xu = (1/(integral_Adams(10,Rate_ccsn_l,H0_al,R0_al)*10**9))
x=(1/(integral_Adams(10,Rate_ccsn,H0_al,R0_al)*10**9))
xd = (1/(integral_Adams(10,Rate_ccsn_u,H0_al,R0_al)*10**9))

print ("+ %s" %(xu-x))
print ("Adams CCSN is %s" %(x))
print ("- %s" %(x-xd))


# In[10]:


xu=(1/(integral_Adams(Rad_lgrb,R_lgrb_l,H0_al,R0_al)*10**9))
x=(1/(integral_Adams(Rad_lgrb,R_lgrb,H0_al,R0_al)*10**9))
xd=(1/(integral_Adams(Rad_lgrb,R_lgrb_u,H0_al,R0_al)*10**9))

print ("+ %s" %(xu-x))
print ("Adams LGRB is %s" %(x))
print ("- %s" %(x-xd))


# In[11]:


x1 = integral_Adams(10,R_sn_u*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn_u*0.5,H0_al, R0_al)
y1 = integral_Adams(10,R_sn*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn*0.5,H0_al, R0_al)
z1 = integral_Adams(10,R_sn_l*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn_l*0.5,H0_al, R0_al)

ia_u = (1/((x1)*10**9))
ia = (1/((y1)*10**9))
ia_l = (1/((z1)*10**9))

ia_upper = ia-ia_u
ia_lower = ia_l-ia

print ("Adams TypeIa lower is + %s" %((ia_lower)))
print ("Adams TypeIa is %s" %ia)
print ("Adams TypeIa upper is - %s" %((ia_upper)))


# # Fe60 Model Mean Times Calculations

# In[12]:


xu = (1/(integral_Adams(10,Rate_ccsn_l,H0_fe,R0_fe)*10**9))
x=(1/(integral_Adams(10,Rate_ccsn,H0_fe,R0_fe)*10**9))
xd = (1/(integral_Adams(10,Rate_ccsn_u,H0_fe,R0_fe)*10**9))

print ("+ %s" %(xu-x))
print ("Adams CCSN is %s" %(x))
print ("- %s" %(x-xd))


# In[13]:


xu=(1/(integral_Adams(Rad_lgrb,R_lgrb_l,H0_fe,R0_fe)*10**9))
x=(1/(integral_Adams(Rad_lgrb,R_lgrb,H0_fe,R0_fe)*10**9))
xd=(1/(integral_Adams(Rad_lgrb,R_lgrb_u,H0_fe,R0_fe)*10**9))

print ("+ %s" %(xu-x))
print ("Adams LGRB is %s" %(x))
print ("- %s" %(x-xd))


# In[14]:


x1 = integral_Adams(10,R_sn_u*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn_u*0.5,H0_fe, R0_fe)
y1 = integral_Adams(10,R_sn*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn*0.5,H0_fe, R0_fe)
z1 = integral_Adams(10,R_sn_l*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn_l*0.5,H0_fe, R0_fe)

ia_u = (1/((x1)*10**9))
ia = (1/((y1)*10**9))
ia_l = (1/((z1)*10**9))

ia_upper = ia-ia_u
ia_lower = ia_l-ia

print ("Adams TypeIa lower is + %s" %((ia_lower)))
print ("Adams TypeIa is %s" %ia)
print ("Adams TypeIa upper is - %s" %((ia_upper)))


# # Rate of Explosions at the Galactic Center

# In[ ]:


x1 = integral_Adams(10,Rate_ccsn,H0_thick, R0_thick)*10**9

y1 = integral_Adams(10,Rate_ccsn,H0_thick, R0_thick)*10**9

z1 = integral_Adams(10,Rate_ccsn,H0_thick, R0_thick)*10**9

ia_u = x1
ia = y1
ia_l = z1

ia_upper = ia-ia_u
ia_lower = ia_l-ia


print ("Adams CCSN lower is + %s" %((ia_lower)))
print ("Adams CCSN is %s" %ia)
print ("Adams CCSN upper is - %s" %((ia_upper)))


# # TRILEGAL Model Mean Times Calculations [Table 5]

# In[8]:


xu = (1/(integral_Adams(10,Rate_ccsn_l,H0_thin,R0_thin)*10**9))
x=(1/(integral_Adams(10,Rate_ccsn,H0_thin,R0_thin)*10**9))
xd = (1/(integral_Adams(10,Rate_ccsn_u,H0_thin,R0_thin)*10**9))

print ("+ %s" %(xu-x))
print ("Adams CCSN is %s" %(x))
print ("- %s" %(x-xd))


# In[33]:


xu=(1/(integral_Adams(Rad_lgrb,R_lgrb_l,H0_thin,R0_thin)*10**9))
x=(1/(integral_Adams(Rad_lgrb,R_lgrb,H0_thin,R0_thin)*10**9))
xd=(1/(integral_Adams(Rad_lgrb,R_lgrb_u,H0_thin,R0_thin)*10**9))

print ("+ %s" %(xu-x))
print ("Adams LGRB is %s" %(x))
print ("- %s" %(x-xd))


# In[34]:


xu=(1/(integral_Adams(Rad_sgrb,R_sgrb_l,H0_thick,R0_thick)*10**9))
x=(1/(integral_Adams(Rad_sgrb, R_sgrb, H0_thick,R0_thick)*10**9))
xd=(1/(integral_Adams(Rad_sgrb,R_sgrb_u, H0_thick,R0_thick)*10**9))

print ("+ %s" %(xu-x))
print ("Adams SGRB is %s" %(x))
print ("- %s" %(x-xd))


# In[20]:


x1 = integral_Adams(10,R_sn_u*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn_u*0.5,H0_thin, R0_thin)
y1 = integral_Adams(10,R_sn*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn*0.5,H0_thin, R0_thin)
z1 = integral_Adams(10,R_sn_l*0.5,H0_thick, R0_thick)+integral_Adams(10,R_sn_l*0.5,H0_thin, R0_thin)

ia_u = (1/((x1)*10**9))
ia = (1/((y1)*10**9))
ia_l = (1/((z1)*10**9))

ia_upper = ia-ia_u
ia_lower = ia_l-ia

print ("Adams TypeIa lower is + %s" %((ia_lower)))
print ("Adams TypeIa is %s" %ia)
print ("Adams TypeIa upper is - %s" %((ia_upper)))


# # SNR Model Mean Times Calculations [Table 5]

# In[13]:


x1 = integral_Green(10,Rate_ccsn_l,H0_thin, R_green)
y1 = integral_Green(10,Rate_ccsn, H0_thin, R_green)
z1 = integral_Green(10,Rate_ccsn_u, H0_thin, R_green)


ia_l = (1/((x1)*10**9))
ia = (1/((y1)*10**9))
ia_u = (1/((z1)*10**9))

ia_upper = ia-ia_u
ia_lower = ia_l-ia


print ("Green CCSN lower is + %s" %((ia_lower)))
print ("Green CCSN is %s" %ia)
print ("Green CCSN upper is - %s" %((ia_upper)))


# In[20]:


x1 = integral_Green(Rad_lgrb,R_lgrb_l,H0_thin, R_green)
y1 = integral_Green(Rad_lgrb,R_lgrb, H0_thin, R_green)
z1 = integral_Green(Rad_lgrb,R_lgrb_u, H0_thin, R_green)


ia_l = (1/((x1)*10**9))
ia = (1/((y1)*10**9))
ia_u = (1/((z1)*10**9))

ia_upper = ia-ia_u
ia_lower = ia_l-ia


print ("Green LGRB lower is + %s" %((ia_lower)))
print ("Green LGRB is %s" %ia)
print ("Green LGRB upper is - %s" %((ia_upper)))


# In[22]:


x1 = integral_Green(Rad_sgrb,R_sgrb_l,H0_thick, R_green)
y1 = integral_Green(Rad_sgrb,R_sgrb, H0_thick, R_green)
z1 = integral_Green(Rad_sgrb,R_sgrb_u, H0_thick, R_green)


ia_l = (1/((x1)*10**9))
ia = (1/((y1)*10**9))
ia_u = (1/((z1)*10**9))

ia_upper = ia-ia_u
ia_lower = ia_l-ia


print ("Green SGRB lower is + %s" %((ia_lower)))
print ("Green SGRB is %s" %ia)
print ("Green SGRB upper is - %s" %((ia_upper)))


# In[67]:


x1 = integral_Green(10,R_sn_u*0.5,H0_thick, R_green)+integral_Green(10,R_sn_u*0.5,H0_thin, R_green)

y1 = integral_Green(10,R_sn*0.5,H0_thick, R_green)+integral_Green(10,R_sn*0.5,H0_thin, R_green)

z1 = integral_Green(10,R_sn_l*0.5,H0_thick, R_green)+integral_Green(10,R_sn_l*0.5,H0_thin, R_green)


ia_u = (1/((x1)*10**9))
ia = (1/((y1)*10**9))
ia_l = (1/((z1)*10**9))

ia_upper = ia-ia_u
ia_lower = ia_l-ia


print ("Green TypeIa lower is + %s" %((ia_lower)))
print ("Green TypeIa is %s" %ia)
print ("Green TypeIa upper is - %s" %((ia_upper)))

