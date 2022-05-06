#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 19:37:11 2021

@author: dealderod
"""


import matplotlib.pyplot as plt
import numpy as np

changeV = [-1.288
,0.960
,6.400
,5.724
,4.852
,5.966
,0.051
,-0.179
,-3.497
,-0.597
,-0.652] #* (u.km)/(u.s)

vrest = [18.88465
,21.28542
,19.75805
,19.21888
,20.80521
,18.80889
,22.92494
,25.71514
,25.71518
,23.65747
,21.07067] #* u.GHz

changeVerror= [0.043
,0.020
,0.010
,0.050
,0.021
,0.021
,0.028
,0.110
,0.112
,0.080
,0.012] #* (u.km)/(u.s)

changeA = [0.186
,-0.044
,-0.361
,0.079
,1.719
,0.047
,1.950
,1.700
,2.480
,0.090
,0.888]

changeAerror = [0.0019
,0.0031
,0.0054
,0.0001
,0.0160
,0.0022
,0.0030
,0.0063
,0.0048
,0.0033
,0.0056]

V12 = [0.18
,-0.80
,-0.83
,2.29
,-0.10
,-0.79
,-0.55
,-0.57
,0.16
,-0.33
,-0.58]

V12error = [0.09
,0.01
,0.01
,0.48
,0.01
,0.02
,0.04
,0.02
,0.05
,0.02
,0.02]
plt.suptitle("D.Deal and C.Henkel ammonia masers comparisons", fontname="Times New Roman")
plt.show()
    
plt.subplot(1,3,1)

plt.errorbar(vrest, changeV, yerr=changeVerror, fmt="o", capsize=4, color="black")
plt.title('W51North', fontname="Times New Roman")
plt.xlabel('Rest Frequency (GHz)',fontname="Times New Roman")
plt.ylabel('∆V_LSR (km s-1)',fontname="Times New Roman")
plt.xticks(np.arange(19,26))

plt.show()

plt.subplot(1,3,2)
plt.errorbar(changeA, changeV, yerr=changeVerror, xerr=changeAerror, fmt="o", capsize=4, color="black")
plt.title('W51North',fontname="Times New Roman")
plt.xlabel('∆S (Jy)',fontname="Times New Roman")
#plt.ylabel('∆V_LSR (km s-1)')
plt.show()

plt.subplot(1,3,3)
plt.errorbar(V12, changeV, yerr=changeVerror, xerr=V12error, fmt="o", capsize=4, color="black")
plt.title('W51North',fontname="Times New Roman")
plt.xlabel('∆V_1/2 (km s-1)',fontname="Times New Roman")
#plt.ylabel('∆V_LSR (km s-1)')
plt.show()
