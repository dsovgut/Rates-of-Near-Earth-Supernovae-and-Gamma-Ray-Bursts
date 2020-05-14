#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches


# In[2]:


#reading from the csv file
distance, LGRB_l, LGRB,LGRB_u,CCSN_l, CCSN,CCSN_u,SGRB_l, SGRB,SGRB_u,TypeIa_l,TypeIa,TypeIa_u = np.loadtxt(fname = 'final_rates2.csv', dtype = 'float', unpack= 'True', usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13), skiprows = 1, delimiter=',')

#splitting TypeIa supernovae into thick and thin disk
thin, thick = np.loadtxt(fname = 'scale_rates2.csv', dtype = 'float', unpack= 'True', usecols = (2,3), skiprows = 1, delimiter=',')
R_sn_l = 0.006 #*(1/u.yr)
R_sn = 0.0139 #*(1/u.yr) #type Ia sn rate 
R_sn_u = 0.028 #*(1/u.yr)
sn_l = []
sn = []
sn_u= []
for thin,thick in zip(thin, thick):
    sn_l.append((thin*R_sn_l*0.5)+(thick*R_sn_l*0.5))
    sn.append((thin*R_sn*0.5)+(thick*R_sn*0.5))
    sn_u.append((thin*R_sn_u*0.5)+(thick*R_sn_u*0.5))


# In[14]:


print ('Fe60 Rate: %s [1/yr] - %s [1/yr]' %((1/3e6),(1/2e6)))
print ('AD 1006 Rate: %s [1/yr] - %s [1/yr]' %((1/1014),(1/1011)))
print ('AD 1054 Rate: %s [1/yr] - %s [1/yr]' %((1/966),(1/964.25)))
print ('AD 1181 Rate: %s [1/yr] - %s [1/yr]' %((1/839),(1/838.5)))
print ('AD 1572 Rate: %s [1/yr] - %s [1/yr]' %((1/448),(1/446.5)))
print ('AD 1604 Rate: %s [1/yr] - %s [1/yr]' %((1/416),(1/415)))
print ('AD 1700 Rate: %s [1/yr] - %s [1/yr]' %((1/320),(1/358)))


# In[13]:


print ('Fe60 Distance: %s [pc] - %s [pc]' %((20),(150)))
print ('AD 1006 Distance: %s [pc] - %s [pc]' %((1000),(2000)))
print ('AD 1054 Distance: %s [pc] - %s [pc]' %((1820),(2040)))
print ('AD 1181 Distance: %s [pc] - %s [pc]' %((1700),(2300)))
print ('AD 1572 Distance: %s [pc] - %s [pc]' %((2000),(3000)))
print ('AD 1604 Distance: %s [pc] - %s [pc]' %((2500),(7000)))
print ('AD 1700 Distance: %s [pc] - %s [pc]' %(((3.4-0.1)*1000),((3.4+0.3)*1000)))


# In[16]:


#all the disances and rates in the list for each datapoint
Fe60 = [20,150,(1/3e6),(1/2e6)]
AD1006 = [1000,2000,0.0009861932938856016,0.0009891196834817012]
AD1054 = [1820,2040,0.0010351966873706005, 0.0010370754472387865]
AD1181 = [1700,2300, 0.0011918951132300357,0.0011926058437686344]
AD1572 = [2000,3000,0.002232142857142857,0.0022396416573348264]
AD1604 = [2500,7000,0.002403846153846154,0.0024096385542168677]
AD1700 = [3300,3700,0.002793296089385475,0.003125]


# In[39]:


#function that makes a rectangle, without a label
def rec(d1,d2,r1,r2,x):
    plt.plot([d1,d2], [r1,r1], linewidth=1, color = x)
    plt.plot([d1,d2], [r2,r2], linewidth=1, color = x)
    plt.plot([d1,d1], [r1,r2], linewidth=1, color = x)
    plt.plot([d2,d2], [r1,r2], linewidth=1, color = x)


# In[42]:


#function that makes a rectangle, with a label
def rec2(d1,d2,r1,r2,x, l):
    plt.plot([d1,d2], [r1,r1], linewidth=1, color = x, label = l)
    plt.plot([d1,d2], [r2,r2], linewidth=1, color = x)
    plt.plot([d1,d1], [r1,r2], linewidth=1, color = x)
    plt.plot([d2,d2], [r1,r2], linewidth=1, color = x)


# # All of the Data Points: Fe60 and 5 historic supernovae

# In[48]:


plt.plot(distance, CCSN, color = 'blue', label = 'CCSN') #red and blue, m, green 
plt.plot(distance, sn, color = 'red', label = 'Type Ia') #red and blue, m, green 

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, sn_u,sn_l,color='red',alpha=.09, linestyle = ':')


rec2(Fe60[0],Fe60[1],Fe60[2],Fe60[3], 'black','Fe60 Data')
rec2(AD1181[0],AD1181[1],AD1181[2],AD1181[3], 'cornflowerblue', 'Historic CCSN')
rec2(AD1006[0],AD1006[1],AD1006[2],AD1006[3], 'tomato', 'Historic Type Ia')
rec(AD1054[0],AD1054[1],AD1054[2],AD1054[3], 'cornflowerblue')
rec(AD1572[0],AD1572[1],AD1572[2],AD1572[3], 'tomato')
rec(AD1604[0],AD1604[1],AD1604[2],AD1604[3], 'tomato')
rec(AD1700[0],AD1700[1],AD1700[2],AD1700[3], 'cornflowerblue')


plt.ylim(1e-7,0.0035)
plt.xlim(10,7500)



plt.title("Rate vs Distance (Trilegal Model)", fontsize = 15)

plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance (Trilegal Model)", fontsize = 15)



plt.legend()
plt.savefig('Historic-data.png', dpi=300)

#using different shade, add Cas A, show plots side by side and change colors. 

plt.show()


# # Plotting Historic Supernovae

# In[53]:


plt.plot(distance, CCSN, color = 'blue', label = 'CCSN') #red and blue, m, green 
plt.plot(distance, sn, color = 'red', label = 'Type Ia') #red and blue, m, green 

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, sn_u,sn_l,color='red',alpha=.09, linestyle = ':')


rec2(AD1181[0],AD1181[1],AD1181[2],AD1181[3], 'cornflowerblue', 'Historic CCSN')
rec2(AD1006[0],AD1006[1],AD1006[2],AD1006[3], 'tomato', 'Historic Type Ia')
rec(AD1054[0],AD1054[1],AD1054[2],AD1054[3], 'cornflowerblue')
rec(AD1572[0],AD1572[1],AD1572[2],AD1572[3], 'tomato')
rec(AD1604[0],AD1604[1],AD1604[2],AD1604[3], 'tomato')
rec(AD1700[0],AD1700[1],AD1700[2],AD1700[3], 'cornflowerblue')

plt.ylim(0.00096,0.0035)
plt.xlim(800,7300)



#plt.title("Rate vs Distance (Trilegal Model)", fontsize = 15)

plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
#plt.title("Rate vs Distance (Trilegal Model)", fontsize = 15)


plt.tight_layout()

plt.legend()
#plt.savefig('Historic-SN.png', dpi=300)

plt.show()


# # Plotting individual Supernovae 

# In[131]:


plt.plot(distance, CCSN, color = 'blue', label = 'CCSN') #red and blue, m, green 
plt.plot(distance, sn, color = 'red', label = 'Type Ia') #red and blue, m, green 

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, sn_u,sn_l,color='red',alpha=.09, linestyle = ':')


rec(AD1006[0],AD1006[1],AD1006[2],AD1006[3], 'black')

plt.ylim(0.00097,0.00099)
plt.xlim(800,7300)


plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance (Trilegal Model), AD1006", fontsize = 15)



plt.legend()
#plt.savefig('Rate_vs_Distance_witherrors.png', dpi=300)

plt.show()


# In[132]:


plt.plot(distance, CCSN, color = 'blue', label = 'CCSN') #red and blue, m, green 
plt.plot(distance, sn, color = 'red', label = 'Type Ia') #red and blue, m, green 

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, sn_u,sn_l,color='red',alpha=.09, linestyle = ':')


rec(AD1054[0],AD1054[1],AD1054[2],AD1054[3], 'black')

plt.ylim(0.0009,0.0011)
plt.xlim(1800,3000)


plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance (Trilegal Model), AD1054 (Crab)", fontsize = 15)



plt.legend()
#plt.savefig('Rate_vs_Distance_witherrors.png', dpi=300)

plt.show()


# In[133]:


plt.plot(distance, CCSN, color = 'blue', label = 'CCSN') #red and blue, m, green 
plt.plot(distance, sn, color = 'red', label = 'Type Ia') #red and blue, m, green 

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, sn_u,sn_l,color='red',alpha=.09, linestyle = ':')


rec(AD1181[0],AD1181[1],AD1181[2],AD1181[3], 'black')

plt.ylim(0.0011918951132300357,0.0011926058437686344)
plt.xlim(1600,2400)


plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance (Trilegal Model), AD1181", fontsize = 15)



plt.legend()
#plt.savefig('Rate_vs_Distance_witherrors.png', dpi=300)

plt.show()


# In[134]:


plt.plot(distance, CCSN, color = 'blue', label = 'CCSN') #red and blue, m, green 
plt.plot(distance, sn, color = 'red', label = 'Type Ia') #red and blue, m, green 

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, sn_u,sn_l,color='red',alpha=.09, linestyle = ':')


rec(AD1572[0],AD1572[1],AD1572[2],AD1572[3], 'black')

plt.ylim(0.0021,0.0023)
plt.xlim(1900,3200)


plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance (Trilegal Model), AD1572 (Tycho)", fontsize = 15)



plt.legend()
#plt.savefig('Rate_vs_Distance_witherrors.png', dpi=300)

plt.show()


# In[135]:


plt.plot(distance, CCSN, color = 'blue', label = 'CCSN') #red and blue, m, green 
plt.plot(distance, sn, color = 'red', label = 'Type Ia') #red and blue, m, green 

plt.fill_between(distance, CCSN_l,CCSN_u,color='blue',alpha=.09,linestyle = '--')
plt.fill_between(distance, sn_u,sn_l,color='red',alpha=.09, linestyle = ':')


rec(AD1604[0],AD1604[1],AD1604[2],AD1604[3], 'black')

plt.ylim(0.0023,0.0025)
plt.xlim(2300,8000)


plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance (Trilegal Model), AD1604", fontsize = 15)



plt.legend()
#plt.savefig('Rate_vs_Distance_witherrors.png', dpi=300)

plt.show()


# In[ ]:




