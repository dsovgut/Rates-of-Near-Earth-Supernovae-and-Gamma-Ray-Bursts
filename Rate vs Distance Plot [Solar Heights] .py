#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# In[2]:


distance,CCSN_z0,CCSN_z25,CCSN_z75,CCSN_z100,CCSN_z150  = np.loadtxt(fname = 'CCSN_zheights.csv', dtype = 'float', unpack= 'True', usecols = (1,2,3,4,5,6), skiprows = 1, delimiter=',')
distance,sn_z0,sn_z25,sn_z75,sn_z100,sn_z150  = np.loadtxt(fname = 'SN_zheights.csv', dtype = 'float', unpack= 'True', usecols = (1,2,3,4,5,6), skiprows = 1, delimiter=',')


# In[7]:


from pylab import *

line1, = plt.plot(distance, CCSN_z0, color = 'blue', label = 'CCSN, z = 0 pc', linewidth=0.8) #red and blue, m, green 
line2, = plt.plot(distance, CCSN_z25, color = 'red', label = 'CCSN, z = 25 pc', linewidth=0.8) #red and blue, m, green 
line3, = plt.plot(distance, CCSN_z75, color = 'green', label = 'CCSN, z = 75 pc', linewidth=0.8)
line4, = plt.plot(distance, CCSN_z100, color = 'm', label = 'CCSN, z = 100 pc', linewidth=0.8)
line5, = plt.plot(distance, CCSN_z150,  color = 'tomato', label = 'CCSN, z = 150 pc', linewidth=0.8)


line11, = plt.plot(distance, sn_z0, color = 'blue',  label = 'Type Ia, z = 0 pc', linestyle = '--', linewidth=0.8) #red and blue, m, green 
line22, = plt.plot(distance, sn_z25, color = 'red', label = 'Type Ia, z = 25 pc', linestyle = '--', linewidth=0.8) #red and blue, m, green 
line33, = plt.plot(distance, sn_z75, color = 'green', label = 'Type Ia, z = 75 pc', linestyle = '--', linewidth=0.8)
line44, = plt.plot(distance, sn_z100, color = 'm', label = 'Type Ia, z = 100 pc', linestyle = '--', linewidth=0.8)
line55, = plt.plot(distance, sn_z150, color = 'tomato', label = 'Type Ia, z = 150 pc', linestyle = '--', linewidth=0.8)


plt.ylim(1e-10,10e-6)
plt.xlim(10,1000)


plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance for different solar heights z", fontsize = 15)


# Create a legend for the first line.
first_legend = plt.legend(handles=[line1, line2, line3, line4, line5], loc='upper left')

# Add the legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)

# Create another legend for the second line.
plt.legend(handles=[line11, line22, line33, line44, line55], loc='lower right')



plt.savefig('ZHeight.png', dpi=300)

plt.show()


# In[24]:


from pylab import *

line1, = plt.plot(distance, CCSN_z0, color = 'blue', label = 'CCSN', linewidth=0.8, linestyle='-') #red and blue, m, green 
plt.plot(distance, CCSN_z25, color = 'blue', linewidth=0.8, linestyle=':') #red and blue, m, green 
plt.plot(distance, CCSN_z75, color = 'blue', linewidth=0.8,linestyle='--')
plt.plot(distance, CCSN_z100, color = 'blue', linewidth=0.8, linestyle='-.')
plt.plot(distance, CCSN_z150,  color = 'blue', linewidth=0.8, linestyle='--', dashes=[10, 5, 20, 5], label='66')


line11, = plt.plot(distance, sn_z0, color = 'red',  label = 'Type Ia',  linewidth=0.8, linestyle='-') #red and blue, m, green 
plt.plot(distance, sn_z25, color = 'red',  linewidth=0.8, linestyle=':') #red and blue, m, green 
plt.plot(distance, sn_z75, color = 'red',  linewidth=0.8, linestyle='--')
plt.plot(distance, sn_z100, color = 'red', linewidth=0.8, linestyle='-.')
plt.plot(distance, sn_z150, color = 'red', linewidth=0.8, linestyle='--', dashes=[10, 5, 20, 5])


line2, = plt.plot(0,0.00000000001, color = 'black',  linewidth=0.8, linestyle='-', label = 'z=0 pc')
line3, = plt.plot(0,0.00000000001, color = 'black',  linewidth=0.8, linestyle=':', label = 'z=25 pc')
line4, = plt.plot(0,0.00000000001, color = 'black',  linewidth=0.8, linestyle='--', label = 'z=75 pc')
line5, = plt.plot(0,0.00000000001, color = 'black',  linewidth=0.8, linestyle='-.', label = 'z=100 pc')
line6, = plt.plot(0,100, color = 'black',  linewidth=0.8, dashes=[10, 5, 20, 5], linestyle='--',   label = 'z=150 pc')




plt.ylim(1e-10,10e-6)
plt.xlim(10,1000)


plt.xscale('log')
plt.yscale('log')


plt.xlabel (r'Distance from Earth $r$ [pc]')
plt.ylabel (r'Enclosed Explosion Rate $\Gamma(r)$ [$\frac{events}{yr}$]')
plt.title("Rate vs Distance for different solar heights z", fontsize = 15)


# Create a legend for the first line.
first_legend = plt.legend(handles=[line2, line3, line4, line5, line6], loc='upper left')

# Add the legend manually to the current Axes.
ax = plt.gca().add_artist(first_legend)

# Create another legend for the second line.
plt.legend(handles=[line11, line1], loc='lower right')



plt.savefig('ZHeight.png', dpi=300)

plt.show()


# In[ ]:




