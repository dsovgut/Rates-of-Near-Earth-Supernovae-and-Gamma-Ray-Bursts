#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# In[39]:


#loading the Trilegal data to get the thick component for the Type Ia
thin, thick = np.loadtxt(fname = 'scale_rates2.csv', dtype = 'float', unpack= 'True', usecols = (2,3), skiprows = 1, delimiter=',')
#loading the Wang scale data to plot and get the thin compoenent for the Type Ia
Al_scale, Fe_scale = np.loadtxt(fname = 'scale_rates_aluminum_ferum.csv', dtype = 'float', unpack= 'True', usecols = (2,3), skiprows = 1, delimiter=',')
#loading the Wang actual rates
lgrb_al, ccsn_al,sn_al,lgrb_fe,ccsn_fe,sn_fe = np.loadtxt(fname = 'rates-Al-Fe.csv', dtype = 'float', unpack= 'True', usecols = (2,3,4,5,6,7), skiprows = 1, delimiter=',')
#loading the Trilegal data
distance, LGRB_l, LGRB,LGRB_u,CCSN_l, CCSN,CCSN_u,SGRB_l, SGRB,SGRB_u,TypeIa_l,TypeIa,TypeIa_u = np.loadtxt(fname = 'final_rates2.csv', dtype = 'float', unpack= 'True', usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13), skiprows = 1, delimiter=',')


# In[37]:


#making Type Ia rates for the Trilegal Model
R_sn_l = 0.006 #*(1/u.yr)
R_sn = 0.0139 #*(1/u.yr) #type Ia sn rate 
R_sn_u = 0.028 #*(1/u.yr)

sn = []
for thin1,thick1 in zip(thin, thick):
    sn.append((thin1*R_sn*0.5)+(thick1*R_sn*0.5))


# In[41]:


#making Type Ia rates for the Wang model
sn_al2 = []
sn_fe2=[]

for al,fe,thin1,thick1 in zip(Al_scale,Fe_scale,thin,thick):
    sn_al2.append((al*R_sn*0.5)+(thick1*R_sn*0.5))
    sn_fe2.append((fe*R_sn*0.5)+(thick1*R_sn*0.5))


# In[125]:


#trilegal model
plt.plot(distance, CCSN, color = 'blue') 
plt.plot(distance, sn, color = 'red') 
plt.plot(distance, LGRB, color = 'm')
#iron 60 lines
plt.plot(distance, ccsn_fe, color = 'blue', linestyle=':') 
plt.plot(distance, sn_al2, color = 'red',linestyle=':')  
plt.plot(distance, lgrb_fe, color = 'm',linestyle=':')
#aluminum lines
plt.plot(distance, ccsn_al, color = 'blue', linestyle='--') 
plt.plot(distance, sn_al2, color = 'red', linestyle='--') 
plt.plot(distance, lgrb_fe, color = 'm', linestyle='--')

#dummy lines for generating legend
plt.plot(0, 0.0000001, color = 'blue', label = 'CCSN') 
plt.plot(0, 0.0000001, color = 'red', label = 'Type Ia')  
plt.plot(0, 0.0000001, color = 'm', label = 'Long GRB')

plt.plot(0, 0.0000001, color = 'black', label = 'Trilegal Model') 
plt.plot(0, 0.0000001, color = 'black', linestyle=':', label = r'$\rm ^{60}Fe$') #
plt.plot(0, 0.0000001, color = 'black', linestyle='--', label = r'$\rm ^{26}Al$') 

plt.ylim(1e-10,1)
plt.xlim(10,29915.75721530859)

plt.title("Rate vs Distance (Model Comparison)", fontsize = 15)

plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
#plt.title("Rate vs Distance (Comparing the Data)", fontsize = 15)

plt.legend()
plt.savefig('Rate_vs_Distance_c.png', dpi=300)

plt.show()


# # Generating Subplots

# In[121]:


#trilegal model
plt.plot(distance, CCSN, color = 'blue', label = 'Trilegal Model') 
#iron 60 lines
plt.plot(distance, ccsn_fe, color = 'blue', linestyle=':',label = r'$\rm ^{60}Fe$') 
#aluminum lines
plt.plot(distance, ccsn_al, color = 'blue', linestyle='--',label = r'$\rm ^{26}Al$') 


plt.ylim(1e-10,1)
plt.xlim(10,29915.75721530859)

plt.title("CCSN", fontsize = 15)

plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
#plt.title("Rate vs Distance (Comparing the Data)", fontsize = 15)

plt.legend()
plt.savefig('CCSN_c.png', dpi=300)

plt.show()


# In[122]:


#trilegal model
plt.plot(distance, sn, color = 'red',label = 'Trilegal Model') 
#iron 60 lines
plt.plot(distance, sn_al2, color = 'red',linestyle=':', label = r'$\rm ^{60}Fe$')  
#aluminum lines
plt.plot(distance, sn_al2, color = 'red', linestyle='--', label = r'$\rm ^{26}Al$') 


plt.ylim(1e-10,0.04)
plt.xlim(10,29915.75721530859)

plt.title("Type Ia", fontsize = 15)

plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
#plt.title("Rate vs Distance (Comparing the Data)", fontsize = 15)

plt.legend()
plt.savefig('SN_c.png', dpi=300)

plt.show()


# In[123]:


#trilegal model
plt.plot(distance, LGRB, color = 'm',label = 'Trilegal Model')
#iron 60 lines
plt.plot(distance, lgrb_fe, color = 'm',linestyle=':',label = r'$\rm ^{60}Fe$')
#aluminum lines
plt.plot(distance, lgrb_fe, color = 'm', linestyle='--',label = r'$\rm ^{26}Al$')


plt.ylim(1e-10,1*10**(-6.5))
plt.xlim(10,29915.75721530859)

plt.title("LGRB", fontsize = 15)

plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
#plt.title("Rate vs Distance (Comparing the Data)", fontsize = 15)

plt.legend()
plt.savefig('LGRB_c.png', dpi=300)

plt.show()


# In[ ]:




