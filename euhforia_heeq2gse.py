import glob 
import datetime
import os.path
import argparse
import ConfigParser
import matplotlib.pyplot as plt

import scipy as sp
import scipy.constants

import numpy as np 

import csv
import math
from math import *
import random

import os
import sys
print sys.argv

# -----------------------------------------------------------------------------------------------------
# Program that converts EUHFORIA outputs (v1.0) from HEEQ to GSE coordinates
# Input: EUHFORIA (v2.0) data at Earth
# -----------------------------------------------------------------------------------------------------

print "=====================================================================================:"
print 'Program that converts EUHFORIA output vectors (v0.1) from HEEQ to GSE coordinates'
print 'The converted file format is good for OpenGGCM runs'
print "=====================================================================================:"

event='Event1'

# directory where euh_10min.dsv data are stored
euh_dir = './'+event+'/Event1-euh/' #'./20170904/20170904_spr2/'#

# input file
if not os.path.exists(euh_dir):
    os.makedirs(euh_dir)

input_file = './'+event+'/hsphere_Earth.dsv'
# output file
output_file = euh_dir + 'euh_10min_gse.txt'

# ========================================================
# INPUT ROUTINE
# ========================================================

print 'Reading in EUHFORIA (v1.0) data file...'

#date r[AU] clt[rad] lon[rad] n[1/cm^3] P[Pa] vr[km/s] vclt[km/s] vlon[km/s] Br[nT] Bclt[nT] Blon[nT]
# column format in euhforia data file
data_earth = np.loadtxt(input_file, delimiter=' ',skiprows=1, usecols=[1,2,3,4,5,6,7,8,9,10,11])
date_earth = np.loadtxt(input_file, delimiter=' ', dtype='str',skiprows=1, usecols=[0])

# assigning date to euhforia data
final_dates = []
for dates in date_earth:
	final_dates.append( datetime.datetime.strptime(str(dates), "%Y-%m-%dT%H:%M:%S") )


# ========================================================
# CONVERTING FROM HEEQ TO GSE
# ========================================================

print 'Converting from HEEQ to GSE...'

file = open(output_file,"w+") 
file.write( '# date btot bxgse[nT] bygse[nT] bzgse[nT] vtot vxgse[km/s] vygse[km/s] vzgse[km/s] n[1/cm^3] P[pPa] x[gse] y[gse] z[gse]' + os.linesep)

euhforia_earth_xgse = data_earth[:,0] # xgse=r_heeq   #modified in satfilter.py
euhforia_earth_ygse = data_earth[:,1] # to be checked #modified in satfilter.py
euhforia_earth_zgse = data_earth[:,2] # to be checked #modified in satfilter.py

euhforia_earth_np = data_earth[:,3] # = np [part/cc]
euhforia_earth_pp = data_earth[:,4]*10**12 # = pp [pPa]

euhforia_earth_vx = - data_earth[:,5] # = -v_r
euhforia_earth_vy = - data_earth[:,7] # = -v_long
euhforia_earth_vz = - data_earth[:,6] # = -v_colat

euhforia_earth_bx=-data_earth[:,8]; euhforia_earth_by=-data_earth[:,10]; euhforia_earth_bz=-data_earth[:,9]

'''
# This is to add a random value between -5 and +5 to the whole time series.
# This can still have zero values in the flat part of the magnetic field time series.
for bb in range(len(euhforia_earth_vx)):
	euhforia_earth_bx.append(- data_earth[bb,8] )#+ random.randint(-50,51)/10) # = -b_r
	euhforia_earth_by.append(- data_earth[bb,10] )#+ random.randint(-50,51)/10) # = -b_long
	euhforia_earth_bz.append(- data_earth[bb,9] - 2)#+ random.randint(-50,51)/10) # = -b_colat
'''
# This is to add a positive (negative) random value when the magnetic field values are positive (negative) and zero or close to zero.
# This is to ensure no zero values in the flat part of the magnetic field time series.
# However, this still might give the same results as the previous processing
'''	
threshold=0.005
for bb in range(len(euhforia_earth_vx)):
	val = euhforia_earth_bz[bb]
        if val<threshold and val>(-1*threshold) and final_dates[bb]<datetime.datetime(2012,07,15):
                #euhforia_earth_bx.append(- data_earth[bb,8] + -1)
                #euhforia_earth_by.append(- data_earth[bb,10] + -1)
                euhforia_earth_bz+= -0.5
'''

'''
            if val<0:
                euhforia_earth_bx.append(- data_earth[bb,8]+random.randint(-50,-1)/10)
                euhforia_earth_by.append(- data_earth[bb,10]+random.randint(-50,-1)/10)
                euhforia_earth_bz.append(- data_earth[bb,9]+random.randint(-50,-1)/10)

            if val==0:
                euhforia_earth_bx.append(- data_earth[bb,8]+random.randint(-20,-1)/10)
                euhforia_earth_by.append(- data_earth[bb,10]+random.randint(-20,-1)/10)
                euhforia_earth_bz.append(- data_earth[bb,9]+random.randint(-20,-1)/10)
'''

euhforia_earth_vtot = np.sqrt( data_earth[:,5]**2 + data_earth[:,6]**2 + data_earth[:,7]**2 )
euhforia_earth_btot = np.sqrt( np.array(euhforia_earth_bx)**2 + np.array(euhforia_earth_by)**2 + np.array(euhforia_earth_bz)**2 )

# writing to a file
for i in range(0, len(final_dates)):
	file.write( str(final_dates[i].strftime("%Y-%m-%dT%H:%M:%S"))  + '   ' + format(euhforia_earth_btot[i], '.4f') + '   ' + format(euhforia_earth_bx[i], '.4f') + '   ' +  format(euhforia_earth_by[i], '.4f') + '   ' + format(euhforia_earth_bz[i], '.4f') + '   ' + format(euhforia_earth_vtot[i], '.3f') + '   ' + format(euhforia_earth_vx[i], '.3f')  + '   ' +  format(euhforia_earth_vy[i], '.3f') + '   ' +  format(euhforia_earth_vz[i], '.3f') + '   ' + format(euhforia_earth_np[i], '.3f') + '   ' + format(euhforia_earth_pp[i], '.4f')  + '   ' + format(euhforia_earth_xgse[i], '.4f') + '   ' + format(euhforia_earth_ygse[i], '.4f') + '   ' + format(euhforia_earth_zgse[i], '.4f') + os.linesep)
file.close()

