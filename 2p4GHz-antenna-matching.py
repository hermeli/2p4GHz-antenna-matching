#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Script calculates optimal S11 (reflection) values for a 2.4 GHz
# antenna PI-matching circuit with 0.1pf/0.1nH resolution. 
#
# 2017, wyss@superspider.net
#
# Matching circuit is PI-filter:
# ---------Ls1-------------
#      |          |       |
#     Cp1        Cp2    Zant
#      |          |       |
# -------------------------
#-------------------------------------------------------------------------------
import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath
import pylab as pl
# import smithplot from local git path
sys.path.append("./pySmithPlot")
from smithplot import SmithAxes
pi = math.pi

# initialize value range of components (in Farad/Henry)
Cp1s = pl.frange(0.1e-12,5.0e-12,0.1e-12)
Ls1s = pl.frange(0.1e-9,10e-9,0.1e-9)
Cp2s = pl.frange(0.1e-12,5.0e-12,0.1e-12)

# initialize Zant antenna impedances (measured with a VNA)
Zant = {2*pi*2.4e9:complex(79,67),2*pi*2.44e9:complex(116,-70),2*pi*2.48e9:complex(31,-55)}

# do not change anything below this line!
#-------------------------------------------------------------------------------

# some definitions
Zw = complex(50,0)

# plot Zant antenna impedance in smith chart
plt.figure(figsize=(10, 10))
ax = plt.subplot(1, 1, 1, projection='smith')    
plt.plot(list(Zant.values()),label="Zant: antenna impedances (2.40/2.44/2.48 GHz)",datatype=SmithAxes.Z_PARAMETER)

# helper functions for complex impedances/admitances and Goal function
def Zres(R,w): return complex(R,0)
def Zind(L,w): return complex(0,w * L)
def Zcap(C,w): return 1.0 / (complex(0,w) * C)
def par(Z1,Z2): return Z1*Z2/(Z1+Z2)
def ser(Z1,Z2): return Z1+Z2
def S11(Z): return abs((Z-Zw)/(Z+Zw))
def Goal(S): return S**2

# the following dictionaries and lists will be used in the iteration loops
S11_temp = {list(Zant.keys())[0]:complex(0,0),list(Zant.keys())[1]:complex(0,0),list(Zant.keys())[2]:complex(0,0)}
Z_temp = {}                 # dictionary w:Z 
Z_opt = {}                  # dictionary w:Z at optimum tuning point
Goal_opt = float("inf")     # scalar variable for best goal value
Parts_opt = {}              # dictionary 'Cp1':valCp1

# calculate and display the total number of iterations
iterations = len(Cp1s)*len(Ls1s)*len(Cp2s)*len(Zant)
print('Total of {:d} iterations to calculate'.format(iterations))

# iterate over all Cp2s, Ls1s and Cp1s
for Cp2 in Cp2s:  
    print("#",end="")
    for Ls1 in Ls1s:
        for Cp1 in Cp1s:
            Goal_new = 0
            for w in Zant:
                Z_temp[w] = par(Zant[w],Zcap(Cp2,w))
                Z_temp[w] = ser(Z_temp[w],Zind(Ls1,w))
                Z_temp[w] = par(Z_temp[w],Zcap(Cp1,w))
                S11_temp[w] = S11(Z_temp[w])
                Goal_new += Goal(S11_temp[w])

            if Goal_new < Goal_opt:
                Goal_opt = Goal_new
                Z_opt = Z_temp.copy()
                Parts_opt = {'Cp1':Cp1/1e-12,'Ls1':Ls1/1e-9,'Cp2':Cp2/1e-12}

# print results to console
print("")
print('Optimum Part Values: Cp1={:.1f}pF Ls1={:.1f}nH Cp2={:.2f}pF S11_min={:.5f}'.format(Parts_opt['Cp1'],Parts_opt['Ls1'],Parts_opt['Cp2'],Goal_opt))                
print('Zopt[low/mid/high]={:.1f},{:.1f},{:.1f}'.format(list(Z_opt.values())[0],list(Z_opt.values())[1],list(Z_opt.values())[2]))

# show results in smith chart
plt.plot(list(Z_opt.values()),label="Zopt: PI-filter input impedances with min. sum of power reflections",datatype=SmithAxes.Z_PARAMETER)
plt.legend(loc="lower right", fontsize=12)
plt.title("Zant to Zopt with Cp1={:.1f}pF Ls1={:.1f}nH Cp2={:.1f}pF".format(Parts_opt['Cp1'],Parts_opt['Ls1'],Parts_opt['Cp2']))
plt.show()

