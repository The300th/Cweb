#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:12:17 2020

@author: aknebe
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import struct
from   read_Cweb import readCweb  # this is Weiguang's most flexible read routine
from   read_Cweb import read_Cweb # AK's routine to Cweb

def Randomsample_data(data, N):

    rand_pos     = np.random.randint(0,len(data),N)
    data_sampled = data[rand_pos]

    return data_sampled


#==============================================================================
#                            MAIN
#==============================================================================

# just look at a single file
#--------------------------------------------------
Cwebfile = '//Users/aknebe/Office/Analysis/COLA_CW/lcdm_z0p000.00128.Rs=5.00.Cweb'
Cweb     = readCweb(Cwebfile,PWEB=True)

x     = Cweb[:,0]
y     = Cweb[:,1]
z     = Cweb[:,2]
delta = Cweb[:,3]
l1    = Cweb[:,10]
l2    = Cweb[:,11]
l3    = Cweb[:,12]

# for large grids it makes more sense to use a subset of the data for visualization
# rand_pos = np.random.randint(0,len(x),128**3)
# x        = x[rand_pos]
# y        = y[rand_pos]
# z        = z[rand_pos]
# delta    = delta[rand_pos]
# l1       = l1[rand_pos]

fig = plt.figure(figsize=(20,16))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c=delta/max(delta))
plt.show()

'''
# here we compare various files against each other
#--------------------------------------------------
Cwebfiles = ['/Users/aknebe/Office/Source/Analysis/Cweb/GadgetX_NewMDCLUSTER_0001.snap_128.00064.Rs=0.00.Cweb-ascii',
             '/Users/aknebe/Office/Source/Analysis/Cweb/GadgetX_NewMDCLUSTER_0001.snap_128.00064.Rs=35.00.Cweb-ascii',
             '/Users/aknebe/Office/Source/Analysis/Cweb/GadgetX_NewMDCLUSTER_0001.snap_128.00064.Rs=50.00.Cweb-ascii',
             '/Users/aknebe/Office/Source/Analysis/Cweb/GadgetX_NewMDCLUSTER_0001.snap_128.00064.Rs=100.00.Cweb-ascii']

Plotfile = 'Cweb.pdf'
fig = plt.figure(figsize=(20,16))

for i,Cwebfile in enumerate(Cwebfiles):
    Cweb = read_CwebASCII(Cwebfile)
    x     = Cweb[:,0]
    y     = Cweb[:,1]
    z     = Cweb[:,2]
    delta = Cweb[:,3]
    l1    = Cweb[:,10]
    l2    = Cweb[:,11]
    l3    = Cweb[:,12]
    subplot = 221+i
    ax  = fig.add_subplot(subplot, projection='3d')
    ax.scatter(x, y, z, c=l1/max(l1))
    # ax.scatter(x, y, z, c=np.log10(delta)/np.log10(max(delta)))
    ax.set_xlabel('Mpc/h')
    ax.set_ylabel('Mpc/h')
    ax.set_zlabel('Mpc/h')


plt.show()

#plt.savefig(Plotfile)
'''
