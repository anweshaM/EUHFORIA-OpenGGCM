#!/usr/bin/env python
# fix the multiple variable for loop. Currently plots for lower number of iterations out of the two variables. Need to make it independent

import calendar
import glob
import os.path
import argparse
#import configparser

import numpy as np
import scipy as sp
import scipy.constants

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, AutoLocator
#from matplotlib.dates import DayLocator, HourLocator, DateFormatter, strpdate2num
import matplotlib.dates as mdates
from scipy.interpolate import interp1d
import cPickle as pickle #for loading data from npz

import csv
from math import *
import random
import pandas as pd

import os
import sys
print(sys.argv)

#from empirical_ae_functions import *
from compute_dst import *

from datetime import * #datetime

def interpol(pressure, t_pressure, t_dst):
    f = interp1d(np.asarray([calendar.timegm(x.timetuple()) for x in t_pressure]), pressure , kind="linear", bounds_error=False, fill_value=np.nan,)
    t_final = np.asarray([calendar.timegm(x.timetuple()) for x in t_dst])
    return f(t_final)

def interp_nan(A):
    ok = ~np.isnan(A)
    xp = ok.ravel().nonzero()[0]
    fp = A[~np.isnan(A)]
    x  = np.isnan(A).ravel().nonzero()[0]
    A[np.isnan(A)] = np.interp(x, xp, fp)
    return A

def get_og_date_data(name):    

    #Loading OpenGGCM solar wind data (*.bzgse)

    og_date=np.loadtxt(name, dtype='str', usecols=[0,1,2,3,4,5]); og_dates = []; 
    year = og_date[:,0]; month = og_date[:,1]; day = og_date[:,2]; hour = og_date[:,3]; minute = og_date[:,4]; second = og_date[:,5]; 
    for ii in range(len(second)):
        second[ii] = int(float(second[ii]))

    for ii in range(len(year)):    
        res = str(year[ii] + "-" + month[ii] + "-" + day[ii] + "T" + hour[ii] + ":" + minute[ii] + ":" + str(second[ii]))#[:-5])
        og_dates.append(datetime.strptime(str(res), "%Y-%m-%dT%H:%M:%S"))

    og_data=np.loadtxt(name, usecols=[6])

    return og_dates, og_data

mu_0 = 4*np.pi*(10**(-7.))  # vacuum magnetic permeability [SI: N / A^2]
mu = 1              # mean molecular weight of the plasma
m_p = 1.6726 * (10**(-27.))     # proton mass [kg]

# directory where data are stored
sys_dir = '/media/u0141347/T7/openggcm/inputs/'#+spr+'/heliosphere/'
ext_dir = '/media/u0141347/T7/openggcm/Paper/' 
empirical_dir = '/media/u0141347/T7/openggcm/target/empirical_indices/'

#Relevant data files
nn=0 #Change this while changing event
save=np.bool(True)
kyoto=np.bool(True)
supermag=np.bool(False)
cut_dst=np.bool(True)
event=['20120712','20170904']
year_arr=['2012','2017']
year = year_arr[nn]; loc = 'Earth'

sw = '20120712_wi_lowcad_1'#'20170904_omni_lowcad3'#'20170904_lowcad_11'#'20120712_wi_lowres'#'20120712_wi_lowres'#'20170904_omni_lowcad3' is final for 20170904 #'20120712_wi_lowres' #'20170904_omni_lowres' 
spr= '20120712_wi_lowres' #'20120712_initial_-1nT' #'20170904_spr_med_highres' #'20170904_spr_med_highres' #'20120712_spr_135_med_lowres' #'20120712_spr_135_med_highres_pm5nT' 
spr1='20120712_wi_lowcad_init15h' #'20170904_lowcad_11' #'20120712_lowcad' # is final for 20170904 #'20170904_omni_lowcad'
spr2= '20120712_lowcad_init15h-1nT_1'#'20170904_spr_med_highres-2nT'


#For paper:
if nn==1: sw = '20170904_omni_lowcad3'; spr = '20170904_lowcad_11'
if nn==0: sw = '20120712_wi_lowcad_init14h_re'; spr = '20120712_lowcad_init15h_-1nT_1' #20120712_wi_lowcad_init14h
dsv = [sys_dir+event[nn]+'/'+sw, sys_dir+event[nn]+'/'+spr]#, sys_dir+event[nn]+'/'+spr1]#sys_dir+event[nn]+'/'+'wind0', sys_dir+event[nn]+'/'+sw]#,sys_dir+event[nn]+'/'+spr, sys_dir+event[nn]+'/'+spr1]#, sys_dir+event[nn]+'/'+spr2]
dst_files = [sw,spr]#, spr1]#, spr2]
png_name = 'dst_' + sw + '_' + spr

colors = ['k-','b-','m-','g-','c-','y-']
#labels = ['Wind init12h bz=4.06','Wind init21h bz=4.9','Wind init18h bz=-2.56','OMNI','Low cadence','High cadence', 'Bz +/-5nT', 'Bz -2nT']
labels = ['Observations+OpenGGCM','EUHFORIA+OpenGGCM']
plabels = ['Observations','EUHFORIA']
suptit = sw + '+' + spr 

# For loading Dst data from OMNI database
file_arr=['OMNI2_H0_MRG1HR_20120712_20120725', 'OMNI2_H0_MRG1HR_20170905_20170915']
file=file_arr[nn]+".txt" #"OMNI2_H0_MRG1HR_54276.txt" #"OMNI_HRO2_5MIN_201035.txt" #"OMNI_HRO_1MIN_239932_14_to_18.txt"

if nn==1:
    mc_t0 = '2017-09-07T09:13:00'
    mc1_t0 = '2017-09-08T11:00:00'
    mc_t1 = '2017-09-08T04:00:00'
    mc1_t1 = '2017-09-08T20:00:00'

date_beg_arr=["2012-07-14T14:00:00","2017-09-06T12:00:00"]
date_end_arr=["2012-07-16T12:00:00","2017-09-09T04:00:00"]
date_beg = date_beg_arr[nn]
date_end = date_end_arr[nn]

#------------------------------------------------------------------------#
#------------------------------- WIND Data ------------------------------#
#------------------------------------------------------------------------#
'''
  FORMAT OF THE SUBSETTED FILE
    
    ITEMS                      FORMAT   
     
 1 Year                          I4        
 2 Day                           I4        
 3 Hour                          I3        
 4 Minute                        I3        
 5 BX, GSE, nT                   F8.2      
 6 BY, GSE, nT                   F8.2      
 7 BZ, GSE, nT                   F8.2      
 8 KP_Speed, km/s                F8.1      
 9 KP_Vx,km/s                    F8.1      
10 Kp_Vy, km/s                   F8.1      
11 KP_Vz, km/s                   F8.1      
12 Kp_proton Density, n/cc       F7.2      
13 Kp_temperature, K             F9.0
'''
# For plotting solar wind data
omni_vB_dir = sys_dir+event[nn]+'/'
omni_vB_file_arr = ['wind_min_20120712_20120722.lst','wind_min_merge_20170905_20170915.lst']
omni_vB_file=omni_vB_file_arr[nn] 

#Loading OMNI solar wind data
omni_date=np.loadtxt(omni_vB_dir+omni_vB_file, dtype='str', usecols=[0,1,2,3])
omni_vB_dates = [] #omni _dates is a new array

day_num = omni_date[:,1]
hour = omni_date[:,2]
minute = omni_date[:,3]

for ii in range(len(day_num)):
    day_num[ii].rjust(3 + len(day_num[ii]), '0')
    res = datetime.strptime(year + "-" + day_num[ii] + "T" + hour[ii] + ":" + minute[ii] + ":" + "00", "%Y-%jT%H:%M:%S").strftime("%Y-%m-%dT%H:%M:%S")
    omni_vB_dates.append(datetime.strptime(str(res), "%Y-%m-%dT%H:%M:%S"))
    
omni_data0=np.loadtxt(omni_vB_dir+omni_vB_file, usecols=[4,5,6,7,8,9,10,11,12])
print(len(omni_data0))

#prepare for masking arrays - 'conventional' arrays won't do it
omni_data0 = np.ma.array(omni_data0)
omni_data0 = omni_data0.astype('float')
omni_data0[omni_data0 == 99.99] = np.nan
omni_data0[omni_data0 == 999.99] = np.nan
omni_data0[omni_data0 == 9999.99] = np.nan
omni_data0[omni_data0 == 99999.9] = np.nan
omni_data = omni_data0

omni_n = interp_nan(omni_data[:,7])
omni_v = interp_nan(omni_data[:,3]) #km/s
omni_vx = interp_nan(omni_data[:,4]); omni_vy = interp_nan(omni_data[:,5]); omni_vz = interp_nan(omni_data[:,6]) #km/s
omni_rho = mu*m_p*interp_nan(omni_data[:,7])*(10**6.)   # [kg/m3]
omni_bx = interp_nan(omni_data[:,0]); omni_by = interp_nan(omni_data[:,1]); omni_bz = interp_nan(omni_data[:,2]) #nT
omni_b = np.sqrt(interp_nan(omni_data[:,0])**2 + interp_nan(omni_data[:,1])**2 + interp_nan(omni_data[:,2])**2)
ap_ratio = 0.0494
omni_pmag = 10**9*(interp_nan(omni_data[:,0])*10**-9.0)**2.0/(2.0*10**-7) # [nPa]
omni_pram_doug = (1+4*ap_ratio)*omni_rho*omni_v*omni_v*1000*1000*10**9 # [nPa]
omni_ptherm = mu*interp_nan(omni_data[:,7])*scipy.constants.k*0.8e6*1e6*1e9 #multiplied by 1e6 to convert cc to m3 *scipy.constants.k*interp_nan(omni_data[:,9])

#---------------------------------------------------------
#Getting F10 data
#---------------------------------------------------------

F10_file_arr = ['omni2_hr_20120710-20120724_F10.lst','omni2_hr_20170901-20170915_F10.lst']
F10_file = F10_file_arr[nn]
#Loading OMNI solar wind data
F10_date=np.loadtxt(omni_vB_dir+F10_file, dtype='str', usecols=[0,1,2])
F10_data0 = np.loadtxt(omni_vB_dir+F10_file, usecols=[3])
F10_dates = [] #omni _dates is a new array

year = F10_date[:,0]; day_num = F10_date[:,1]; hour = F10_date[:,2]

for ii in range(len(day_num)):
    day_num[ii].rjust(3 + len(day_num[ii]), '0')
    res = datetime.strptime(year[ii] + "-" + day_num[ii] + "T" + hour[ii] + ":" + "00" + ":" + "00", "%Y-%jT%H:%M:%S").strftime("%Y-%m-%dT%H:%M:%S")
    F10_dates.append(datetime.strptime(str(res), "%Y-%m-%dT%H:%M:%S"))

# construct an interpolator for the F10.7 data
f = interp1d(np.asarray([calendar.timegm(x.timetuple()) for x in F10_dates]), F10_data0, kind="linear", bounds_error=False, fill_value=np.nan,)
t_final = np.asarray([calendar.timegm(x.timetuple()) for x in omni_vB_dates])
# approximate plasma data for the same timestamps as magnetic field
F10 = f(t_final)
if F10[F10 == np.nan]: raise ValueError("NaN found in F10")


#omni_pram = interp_nan(omni_data[:,10]) # [nPa] #Flow pressure = (2*10**-6)*Np*Vp**2 nPa (Np in cm**-3, Vp in km/s, subscript "p"

'''
  FORMAT OF THE SUBSETTED FILE

    ITEMS                      FORMAT

 1 YEAR                          I4
 2 DOY                           I4
 3 Hour                          I3
 4 Kp index                      I3
 5 Dst-index, nT                 I6
 6 AE-index, nT                  I5
 7 AL-index, nT                  I6
 8 AU-index, nT                  I6
'''

#------------------------------- OMNI Data ------------------------------#
# For plotting geomagnetic indices
#Description: https://cdaweb.gsfc.nasa.gov/misc/NotesO.html
#Data source: https://cdaweb.gsfc.nasa.gov/cgi-bin/eval2.cgi
indices_dir = '/media/u0141347/T7/openggcm/target/cdaweb/'
indices_file_arr = ['OMNI_HRO2_5MIN_222347_20120710-20120720_CDAWeb.txt','OMNI_HRO2_5MIN_254394_20170901_20170915.txt']
indices_file=indices_file_arr[nn] 

skprws=93

omni_date = np.loadtxt(indices_dir+indices_file, dtype='str', skiprows=skprws, usecols=[0])
omni_time = np.loadtxt(indices_dir+indices_file, dtype='str', skiprows=skprws, usecols=[1]) 
omni_gm_data = np.loadtxt(indices_dir+indices_file, skiprows=skprws, usecols=[2,3,4,5,6,7]) #[2,3,4,5,6,7,8]
omni_gm_data[omni_gm_data == 99999] = np.nan

omni_dates = [] #omni _dates is a new array

for i in range(len(omni_date)):
    time = omni_time[i]
    dt = omni_date[i] + " " + omni_time[i] 
    dt1 = datetime.strptime(str(dt),"%d-%m-%Y %H:%M:%S.%f").strftime("%Y-%m-%dT%H:%M:%S")
    omni_dates.append( datetime.strptime(str(dt1), "%Y-%m-%dT%H:%M:%S") )
#Lesson learnt: always ensure the format of datetime object before plotting

print (omni_dates[0], type(omni_dates[0]),omni_dates[-1],type(omni_dates[-1]))

dst = interp_nan(omni_gm_data[:,4]) #SYM-H index
print ("Checking Dst observations file: ",dst[0:5])
ae_index = interp_nan(omni_gm_data[:,0])
al_index = interp_nan(omni_gm_data[:,1])
au_index = interp_nan(omni_gm_data[:,2])

print ('AE index', ae_index[0], ae_index[-1])


def compute_r2(it,i,yt,y):
    # construct an interpolator 
    t_final = np.asarray([calendar.timegm(x.timetuple()) for x in yt])
    f = interp1d(np.asarray([calendar.timegm(x.timetuple()) for x in it]), i, kind="linear", bounds_error=False, fill_value=np.nan,)
    i = f(t_final)
    f = interp1d(np.asarray([calendar.timegm(x.timetuple()) for x in yt]), y, kind="linear", bounds_error=False, fill_value=np.nan,)
    y = f(t_final)
    correlation_matrix = np.corrcoef(i, y)
    correlation_xy = correlation_matrix[0,1]
    return correlation_xy**2

#print ('wind empirical dst: ', type(wind_emp_dst), len(wind_emp_dst))# , (wind_emp_dst[0]))

#------------------------------ OpenGGCM FR Data ----------------------------#

def get_og_dst(spr): #name = dsv[ii] == ext_dir+spr+'/'+'dst_'+spr.'txt'
    file=ext_dir+spr+'/'+'dst_'+spr+'.txt'
    datetime_og_dst = np.loadtxt(file, dtype='str', usecols=[1])
    og_dps_mhd_rcm= np.loadtxt(file, usecols=[21]) 
    og_delb_rcm= np.loadtxt(file, usecols=[22]) #delb_rcm
    og_delb_rcm2= np.loadtxt(file, usecols=[23]) #delb_rcm2
    og_dps_rcm= np.loadtxt(file, usecols=[2])  

    final_dates = []
    final_dates_dst = []

    for ii in range(len(datetime_og_dst)):
        dt = datetime_og_dst[ii] # + "T" + datetime_og_dst[ii]
        try:
            datetime.strptime(str(dt),"%Y:%m:%d:%H:%M:%S.%f")
            final_dates_dst.append(datetime.strptime(str(dt),"%Y:%m:%d:%H:%M:%S.%f"))
        except:
            final_dates_dst.append(datetime.strptime(str(dt),"%Y:%m:%d:%H:%M:%S"))

    #for ii in range(len(final_dates_dst)):
    #    dt1 = datetime.strptime(str(dt),"%Y:%m:%d:%H:%M:%S.%f").strftime("%Y-%m-%dT%H:%M:%S")
    #    final_dates_dst[ii] = datetime.strptime(str(dt1), "%Y-%m-%dT%H:%M:%S")

    return final_dates_dst, og_delb_rcm

def get_og_ae(spr): #name = dsv[ii] == ext_dir+spr+'/'+'dst_'+spr.'txt'
    file=ext_dir+spr+'/'
    date_og = np.loadtxt(file+'ae.txt', dtype='str', skiprows=7, usecols=[0])
    time_og = np.loadtxt(file+'ae.txt', dtype='str', skiprows=7, usecols=[1])
    og_ae = np.loadtxt(file+'ae.txt',  skiprows=7, usecols=[2])
    og_al = np.loadtxt(file+'al.txt',  skiprows=7, usecols=[2])
    og_au = np.loadtxt(file+'au.txt',  skiprows=7, usecols=[2])

    final_dates = []

    for ii in range(len(date_og)):
        dt = date_og[ii] + "T" + time_og[ii]
        try:
            datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S.%f")
            final_dates.append(datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S.%f"))
        except: 
            final_dates.append(datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S"))

    return final_dates, og_ae, og_al, og_au

def get_emp_ae(spr): #name = dsv[ii] == ext_dir+spr+'/'+'dst_'+spr.'txt'
    file=sys_dir+event[nn]+'/AUAL_prd_med5.txt'
    date_og = np.loadtxt(file, dtype='str', skiprows=1, usecols=[0,1,2,3,4,5])
    data_og = np.loadtxt(file, dtype='str', skiprows=1, usecols=[6,7])    
    og_al = np.loadtxt(file,  skiprows=1, usecols=[7])
    og_au = np.loadtxt(file,  skiprows=1, usecols=[6])
    og_ae = []; og_ae = [og_au[ii] - og_al[ii] for ii in range(len(og_au))]

    final_dates = []
    yy = date_og[:,0]; month = date_og[:,1]; day = date_og[:,2]; hour = date_og[:,3]; minute=date_og[:,4]; second=date_og[:,5] 
    
    for ii in range(len(month)):
        dt = yy[ii] + '-' + month[ii] + '-' + day[ii] + "T" + hour[ii] + ':' + minute[ii] + ':' + second[ii]
        try:
            datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S")
            final_dates.append(datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S"))
        except: 
            final_dates.append(datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S"))

    return final_dates, og_ae

def closest_time(final_dates_dst,omni_dates):
    cloz_dict = [abs((final_dates_dst[0] - date).total_seconds()) for date in omni_dates]
    start_date = omni_dates[cloz_dict.index(min(cloz_dict))+1] #+1 to take the value just after the first one to ensure you are within the bounds

    cloz_dict = [abs((final_dates_dst[-1] - date).total_seconds()) for date in omni_dates]
    end_date = omni_dates[cloz_dict.index(min(cloz_dict))-1] #-1 to take the value just before the last one to ensure you are within the bounds

    return start_date, end_date


# For plotting SMC data
sme_dir = sys_dir+event[nn]+'/'

def extract_sme(sme_dir):
    data_dir=sme_dir+'smc_data.txt'
    sme_dates=np.loadtxt(data_dir,dtype='str',skiprows=105,usecols=[0,1,2,3,4,5])
    year=sme_dates[:,0];    month=sme_dates[:,1];    day=sme_dates[:,2]
    hour=sme_dates[:,3];    minute=sme_dates[:,4];

    sme_dates_arr=[]

    for ii in range(len(day)):
        #res = datetime.strptime(year[ii] + "-" + month[ii] + "-" + day[ii] + "T" + hour[ii] + ":" + minute[ii] + ":" + "00", "%Y-%m-%dT%H:%M:%S").strftime("%Y-%m-%dT%H:%M:%S")
        #sme_dates_arr.append(str(res))

        dt = year[ii] + '-' + month[ii] + '-' + day[ii] + "T" + hour[ii] + ':' + minute[ii] + ':00'
        try:
            datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S")
            sme_dates_arr.append(datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S"))
        except: 
            sme_dates_arr.append(datetime.strptime(str(dt),"%Y-%m-%dT%H:%M:%S"))

    sme_data=np.loadtxt(data_dir,skiprows=105,usecols=[6,7,8])
    sme=sme_data[:,0]; sml=sme_data[:,1]; smu=sme_data[:,2];

    return sme_dates_arr, sme, sml, smu

#sme_dates_arr, sme, sml, smu=extract_sme(sme_dir)
#plt.plot(sme_dates_arr, sme)
#plt.show()


if nn==0:
    shock_t = '2012-07-14 17:39:56'  
    mc_t0 = '2012-07-15 06:14:22'
    mc_t1 = '2012-07-17 03:21:40'

if nn==1:
    shock_t = '2017-09-06 23:13:00'  #2017 09/06 22:21 https://wind.nasa.gov/ICME_catalog/ICME_catalog_viewer.php
    shock1_t = '2017-09-07 22:38:00'
    mc_t0 = '2017-09-07 09:13:00'
    mc1_t0 = '2017-09-08 11:00:00'
    mc_t1 = '2017-09-08 04:00:00'
    mc1_t1 = '2017-09-08 20:00:00'


line=1
def draw_vertical(line,nn):
    lw_loc = 1.0
    if line==1:
            plt.axvline(x=shock_t,c='m',linestyle='-',linewidth=lw_loc)
            plt.axvline(x=mc_t0,c='g',linestyle='-',linewidth=lw_loc)
            plt.axvline(x=mc_t1,c='g',linestyle='-',linewidth=lw_loc)

            if nn==1:
                plt.axvline(x=shock1_t,c='m',linestyle='--',linewidth=lw_loc)
                plt.axvline(x=mc1_t0,c='g',linestyle='--',linewidth=lw_loc)
                plt.axvline(x=mc1_t1,c='g',linestyle='--',linewidth=lw_loc)



fig = plt.figure(figsize=(10, 12), dpi=80)

num = 0

ae_max = []; t_ae_max = []
dst_min = []; t_dst_min = []
lw = 1.0; marker_ind=1.0
ytic = 10; ytit = 12; leg_size = 11


for i in range(len(dsv)): 

    # MAKING PLOTS
    # =================================================================================================

    #-----------------------------
    start_date=datetime.strptime(date_beg, "%Y-%m-%dT%H:%M:%S")
    end_date=datetime.strptime(date_end, "%Y-%m-%dT%H:%M:%S") 
    
 
    #-----------------------------
    ax = plt.subplot(7,1,1)
    #-----------------------------
    print (dsv[i])
    for name in glob.glob(dsv[i]+'/*.vxgse'): print(name)
    og_dates, og_vx = get_og_date_data(name)

    for name in glob.glob(dsv[i]+'/*.vygse'): print(name)
    og_dates, og_vy = get_og_date_data(name)

    for name in glob.glob(dsv[i]+'/*.vzgse'): print(name)
    og_dates, og_vz = get_og_date_data(name)
    og_v = np.sqrt(og_vx**2 + og_vy**2 + og_vz**2)

    #---------------------------------------------------------------------------------------------------------------------------------
    if nn==0:
        plt.ylim(250,1000)
    else:
        plt.ylim(250,1000)
    ax.set_xlim((start_date, end_date))
    plt.ylim(200,1200)
    ax.yaxis.set_label_text("$\\bf{v [km/s]}$", fontsize=ytit, fontweight="bold")

    #if num==0: ax.plot(omni_vB_dates, omni_v, colors[i], linestyle='-', linewidth=lw, marker='.',markersize=4)

    ax.plot(og_dates, og_v, colors[i], linewidth=lw, label=plabels[i])
    ax.xaxis.set_minor_locator( AutoMinorLocator(6) )
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.legend(prop={'size':leg_size}, loc=1, ncol=4)
    plt.grid()
    draw_vertical(line,nn)

    #plt.text(0.1, 1.1, suptit, fontsize=12, fontweight="bold", ha='left', va='top', transform=ax.transAxes)  

    
    #-----------------------------
    ax = plt.subplot(7,1,2)
    #-----------------------------

    for name in glob.glob(dsv[i]+'/*.np'):
        print(name)
    og_dates, og_np = get_og_date_data(name)
    #---------------------------------------------------------------------------------------------------------------------------------
    og_rho = mu*m_p*og_np*(10**6.) # [kg/m3]
    ap_ratio = 0.0494 #ap_ratio is ...
    og_pram = (1+4*ap_ratio)*og_rho*og_v*og_v*1000*1000*10**9

    ax.set_xlim((start_date, end_date))
    ax.yaxis.set_label_text("$\\bf{n_p [/cc]}$", fontsize=ytit, fontweight="bold")
    if nn==0:
        plt.ylim(0,50)
    else:
        plt.ylim(0,50)
    #if num==0: ax.plot(omni_vB_dates, omni_n, colors[i], linestyle='-', linewidth=lw, marker='.',markersize=4) 
    ax.plot(og_dates, og_np, colors[i], linewidth=lw, label=plabels[i])
    ax.xaxis.set_minor_locator( AutoMinorLocator(6) )
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.grid()
    draw_vertical(line,nn)
    
    #-----------------------------
    ax = plt.subplot(7,1,3)
    #-----------------------------

    #for name in glob.glob(dsv[i]+'/*.bxgse'): print(name)
    #og_dates, og_bx = get_og_date_data(name)
    #for name in glob.glob(dsv[i]+'/*.bygse'): print(name)
    #og_dates, og_by = get_og_date_data(name)
    for name in glob.glob(dsv[i]+'/*.bzgse'): print(name)
    og_dates, og_bz = get_og_date_data(name)
    #og_b = np.sqrt(og_bx**2 + og_by**2 + og_bz**2)
    #---------------------------------------------------------------------------------------------------------------------------------

    ax.set_xlim((start_date, end_date))
    if nn==00:
        plt.ylim(-30,20)
    else:
        plt.ylim(-50,20)

    ax.yaxis.set_label_text("$\\bf{B_z [nT]}$", fontsize=ytit, fontweight="bold")

    ax.plot(og_dates, og_bz, colors[i], linewidth=lw, label=plabels[i])
    ax.xaxis.set_minor_locator( AutoMinorLocator(6) )
    #if num==0: ax.plot(omni_vB_dates, omni_bz, colors[i], linestyle='-', linewidth=lw, marker='.',markersize=4)
    plt.axhline(y=0,linewidth=lw/2)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.grid()
    draw_vertical(line,nn)

    

    #-----------------------------
    ax = plt.subplot(7,1,4)
    #-----------------------------

    emp_dst_bool=np.bool(False)
    plt.xlim((start_date, end_date))
    plt.ylim((-150,100))

    if emp_dst_bool:
            emp_dst_omni = compute_dst(omni_vB_dates, omni_n, omni_v, omni_b, omni_vx, omni_vy, omni_vz, omni_bx, omni_by, omni_bz, start_date, end_date)
            emp_dst_sim = compute_dst(og_dates, og_np, og_v, og_b, og_vx, og_vy, og_vz, og_bx, og_by, og_bz, start_date, end_date)
            ax.plot(og_dates, emp_dst_sim, colors[i], linestyle='--', linewidth=0.3, label='Emp_'+labels[i], marker='.',markersize=2)
            
    ax.set_xlim((start_date, end_date))
    ax.yaxis.set_label_text("$\\bf{Dst [nT]}$", fontsize=ytit, fontweight="bold")
    
    if num==0:
        ax.plot(omni_dates, dst, 'r', linestyle=None, linewidth=lw, label='Measured Dst', marker='.',markersize=marker_ind)
        #ax.plot(omni_vB_dates, emp_dst_omni, '-m', linestyle='--', linewidth=0.3, label='WIND_Emp', marker='.',markersize=2)
    #print (min(dst), np.array(dst).index(min(dst)))

    min_dst_obs = min(dst); min_dst_pos = np.where(dst == dst.min())[0][0]
    print ('Minimum Dst (observations)= ', min(dst), np.where(dst == dst.min())[0][0], omni_dates[min_dst_pos])

    if glob.glob(ext_dir+dst_files[i]+'/'+'dst*'):
        final_dates_dst, og_delb_rcm = get_og_dst(dst_files[i])

        ind_cut=len(final_dates_dst)

        if nn==1 and cut_dst==True: 
            date_cut = datetime.strptime('2017-09-09T05:00:00', "%Y-%m-%dT%H:%M:%S")

            for df in range(len(final_dates_dst)):
                if final_dates_dst[df] > date_cut: 
                    ind_cut = df
                    print (final_dates_dst[df],final_dates_dst[df-1])
                    break
            
        if nn==0 and cut_dst==True: 
            date_cut = datetime.strptime('2012-07-16T13:00:00', "%Y-%m-%dT%H:%M:%S")

            for df in range(len(final_dates_dst)):
                if final_dates_dst[df] > date_cut: 
                    ind_cut = df
                    print (final_dates_dst[df],final_dates_dst[df-1])
                    break

        print ('ind_cut: ', ind_cut, final_dates_dst[ind_cut-1], dsv[i])
        print('Dates formatting: ',final_dates_dst[0])#,'How to cut datetime array: ',np.where(final_dates_dst>datetime(2017,9,8,23,0))[0][0])

        #ax.plot(final_dates_dst[:ind_cut], og_delb_rcm[:ind_cut], colors[i], linewidth=lw, label=labels[i])

        if num==0: 
            #ax.plot(final_dates, og_delb_rcm, colors[num], linestyle='-',linewidth=0.6, label=labels[num])
            print ('Minimum Dst = ', min(og_delb_rcm[0:-10]), labels[num])
        else:
            #ax.plot(final_dates, og_delb_rcm, colors[num], linewidth=0.3, label=labels[num])
            print ('Minimum Dst = ', min(og_delb_rcm[0:-10]), labels[num])

        print('Checking Dst values: ', og_delb_rcm[5:10],labels[num])
        plt.legend(loc=1,prop={'size':leg_size}, ncol=3)

        P = interpol(og_pram, og_dates, final_dates_dst)
        rcm_corr_doug=np.empty(len(og_delb_rcm)); rcm_corr_ak2=np.empty(len(og_delb_rcm)); rcm_corr_fl=np.empty(len(og_delb_rcm))
        #Fenrich and Luhmaan 1998 params as used by UCB_Doug
        
        for dd in range(len(og_delb_rcm)):

            #Obrien&McPheron params
            a=1; b=7.26; c=11 #AK2
            rcm_corr_ak2[dd] = a*og_delb_rcm[dd] + b*np.sqrt(P[dd]) - c #P is in nPa

            #Fenrich and Luhmaan 1998 params
            a = 1.25; b = 15.8; c = 20.0
            rcm_corr_fl[dd] = a*og_delb_rcm[dd] + b*np.sqrt(P[dd]) - c #P is in nPa

        

        #ax.plot(final_dates_dst[:ind_cut], rcm_corr_ak2[:ind_cut], colors[num], linestyle='--')#,label=labels[i]+' - AK2', markersize=2)
        ax.plot(final_dates_dst[:ind_cut], rcm_corr_ak2[:ind_cut], colors[i], linewidth=lw, label=labels[i])
        plt.legend(ncol=3)
        '''
        I am choosing AK2 over FL as it has the least RMSE as per Table 1 in O'brien and McPherron, 2000b.
        '''

        #ax.plot(final_dates_dst[:ind_cut], rcm_corr_fl[:ind_cut], colors[num], linestyle=':')#,label=labels[i]+' - FL', markersize=2)

        print('#----------------'+labels[num]+'----------------#')
        print('Dst data and OpenGGCM: ',compute_r2(omni_dates, dst, final_dates_dst, og_delb_rcm))
    ymin=-250; ymax=150
    plt.ylim(ymin,ymax)
    ax.xaxis.set_minor_locator( AutoMinorLocator(6) )
    plt.setp(ax.get_xticklabels(), visible=False, fontsize=ytic, fontweight="bold")
    plt.setp(ax.get_yticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    ax.xaxis.set_minor_locator( AutoMinorLocator(5) )

    # Show grid
    plt.grid()
    draw_vertical(line,nn)

    og_delb_rcm = rcm_corr_ak2

    final_dates_dst = final_dates_dst[:ind_cut]; og_delb_rcm = og_delb_rcm[:ind_cut]
    
    #uncommenting ends here for plotting Dst panel

    #Modifying Dst and plotting....
    
    '''
    #Obrien&McPheron params
    #a=1; b=0.20; c=20
    a=1; b=7.26; c=11 #AK2
    for dd in range(len(og_delb_rcm)):
        rcm_corr_ak2[dd] = a*og_delb_rcm[dd] + b*np.sqrt(P[dd]) - c #P is in nPa
    ax.plot(final_dates_dst, rcm_corr_ak2, '.r', linewidth=1.0, label='Obrien&McPheron', markersize=2) #colors[num],labels[num]

    #Fenrich and Luhmaan 1998 params
    a=1; b=0.20; c=20
    for dd in range(len(og_delb_rcm)):
        rcm_corr_fl[dd] = a*og_delb_rcm[dd] + b*np.sqrt(P[dd]) - c #P is in nPa
    ax.plot(final_dates_dst, rcm_corr_fl, '.g', label='Fedrich&Luhmaan+1998', markersize=2)
    '''
    #ax.plot(final_dates_dst, og_delb_rcm, colors[num], linestyle='-', linewidth=2.0, label='Dst_RCM')
    

    
    #ax2 = ax.twinx()
    #ax2.plot(final_dates_dst,b*P,'-m',linewidth=1.0, label='$P_{dyn}$[nPa]')

    #print("min Dst %/100: ",(min(cone_delb_rcm)-min(fr_delb_rcm))/min(cone_delb_rcm))
    #ax.plot(final_dates_cone, cone_dst, color3, linewidth=1.0, label=label7)

    sme_dates, sme, sml, smu = extract_sme(sme_dir)
    sme = smu-sml
    print ('#######################################################')
    print ('SME: ',len(sme_dates),len(sme))
    print ('Checking max and min vals SML: ', max(sml), min(sml))
    print ('Checking max and min vals SMU: ', max(smu), min(smu))
    print ('Checking max and min vals SME: ', max(sme), min(sme))
    print ('#######################################################')

    #-----------------------------
    ax = plt.subplot(7,1,5)
    #-----------------------------
    ax.set_xlim((start_date, end_date))
    ymin=0; ymax=800
    #plt.ylim(ymin,ymax) #((min(ae_index)), max(ae_index))
    
    if num==0: 
        if kyoto: ax.plot(omni_dates, au_index, '-r', linestyle=None, linewidth=0.4, label='Measured AU', marker='.',markersize=marker_ind)
        if supermag: ax.plot(sme_dates, smu, '-m', linestyle=None, linewidth=0.4, label='Measured SMU', marker='.',markersize=marker_ind)
        print ('Maximum AU (Observations) = ', max(au_index), 'SMU: ', max(smu), labels[num])

     
    if glob.glob(ext_dir+dst_files[i]+'/'+'au.*'): 
        print (glob.glob(ext_dir+dst_files[i]+'/'+'au.*'))
        final_dates_ae, og_ae, og_al, og_au = get_og_ae(dst_files[i])

        ind_cut=len(final_dates_ae)

        if nn==1 and cut_dst==True: 
            date_cut = datetime.strptime('2017-09-09T05:00:00', "%Y-%m-%dT%H:%M:%S")

            for df in range(len(final_dates_ae)):
                if final_dates_ae[df] > date_cut: 
                    ind_cut = df
                    break

        if nn==0 and cut_dst==True: 
            date_cut = datetime.strptime('2012-07-16T13:00:00', "%Y-%m-%dT%H:%M:%S")

            for df in range(len(final_dates_ae)):
                if final_dates_ae[df] > date_cut: 
                    ind_cut = df
                    break

        ax.plot(final_dates_ae[:ind_cut], og_au[:ind_cut], colors[i], linestyle='-',linewidth=1.0)#, label=labels[i])

        print('#----------------'+labels[num]+'----------------#')
        print('AE data and OpenGGCM: ', compute_r2(omni_dates,ae_index,final_dates_ae,og_ae))

    plt.grid()
    draw_vertical(line,nn)

    #Show legend
    ax.yaxis.set_label_text("AU [nT]", fontweight="bold", fontsize=ytit)
    ax.yaxis.set_major_locator( MultipleLocator(400) )
    ax.yaxis.set_minor_locator( MultipleLocator(100) )
    '''
    l = np.arange(ymin,ymax, 400)
    x = np.arange(ymin,ymax, 400)
    ax.set_yticks(l)
    ax.set_yticklabels(x)
    '''
    plt.setp(ax.get_yticklabels(), visible=True)
    plt.yticks(fontsize=8)
    ax.xaxis.set_minor_locator( AutoMinorLocator(6) )
    plt.setp(ax.get_xticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.setp(ax.get_yticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.legend(loc=2,prop={'size':leg_size},ncol=1)

    #-----------------------------
    ax = plt.subplot(7,1,6)
    #-----------------------------
    ax.set_xlim((start_date, end_date))
    ymin=-2100; ymax=0
    #plt.ylim(ymin,ymax) #((min(ae_index)), max(ae_index))
    
    if num==0: 
        if kyoto: ax.plot(omni_dates, al_index, '-r', linestyle=None, linewidth=0.4, label='Measured AL', marker='.',markersize=marker_ind)
        if supermag: ax.plot(sme_dates, sml, '-m', linestyle=None, linewidth=0.4, label='Measured SML', marker='.',markersize=marker_ind)
        print ('Minimum AL (Observations) = ', min(al_index), 'SML: ', min(sml), labels[num])

     
    if glob.glob(ext_dir+dst_files[i]+'/'+'al.*'): 
        print (glob.glob(ext_dir+dst_files[i]+'/'+'al.*'))
        final_dates_ae, og_ae, og_al, og_au = get_og_ae(dst_files[i])
        ax.plot(final_dates_ae[:ind_cut], og_al[:ind_cut], colors[i], linestyle='-',linewidth=1.0)#, label=labels[i])

        print('#----------------'+labels[num]+'----------------#')
        print('AE data and OpenGGCM: ', compute_r2(omni_dates,ae_index,final_dates_ae,og_ae))

    #if nn==0:
        #emp_dates, emp_ae = get_emp_ae(dst_files[i])
        #ax.plot(emp_dates, emp_ae, colors[i], linestyle='--',linewidth=1.0, label='Emp AE')


    plt.setp(ax.get_xticklabels(), visible=False)
    print('Starting to make plots...')

    plt.grid()
    draw_vertical(line,nn)

    #Show legend
    ax.yaxis.set_label_text("AL [nT]", fontweight="bold", fontsize=ytit)
    ax.yaxis.set_major_locator( MultipleLocator(600) )
    ax.yaxis.set_minor_locator( MultipleLocator(100) )
    '''
    l = np.arange(ymin,ymax, 600)
    x = np.arange(ymin,ymax, 600)
    ax.set_yticks(l)
    ax.set_yticklabels(x)
    '''
    
    plt.yticks(fontsize=8)
    ax.xaxis.set_minor_locator( AutoMinorLocator(6) )
    plt.setp(ax.get_xticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.setp(ax.get_yticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.legend(loc=4,prop={'size':leg_size},ncol=1)
    

    #-----------------------------
    ax = plt.subplot(7,1,7)
    #-----------------------------
    emp_ae_bool = 0

    if emp_ae_bool==1:        
            emp_ae_flg_omni = compute_ae_flg(omni_vB_dates, omni_n, omni_v, omni_b, omni_bz, start_date, end_date)
            emp_ae_lou_omni = compute_ae_flg(omni_vB_dates, omni_n, omni_v, omni_b, omni_bz, F10, start_date, end_date)
            #emp_ae_sim = compute_ae_flg(final_dates1, euhforia_earth_n, euhforia_earth_v, euhforia_earth_b, -1.0*euhforia_earth_vx, -1.0*euhforia_earth_vy, euhforia_earth_vz, -1.0*euhforia_earth_bx, -1.0*euhforia_earth_by, euhforia_earth_bz, start_date, end_date)

    ax.set_xlim((start_date, end_date))
    ymin=0; ymax=2100
    #plt.ylim(ymin,ymax) #((min(ae_index)), max(ae_index))
    
    if num==0: 
        if kyoto: ax.plot(omni_dates, ae_index, '-r', linestyle=None, linewidth=0.4, label='Measured AE', marker='.',markersize=marker_ind)
        if supermag: ax.plot(sme_dates, sme, '-m', linestyle=None, linewidth=0.4, label='Measured SME', marker='.',markersize=marker_ind)
        
        if emp_ae_bool==1: ax.plot(omni_vB_dates, emp_ae_flg_omni, '-g', linestyle=None, linewidth=0.4, label='OMNI_emp_FLG', marker='.',markersize=2)
        print ('Maximum AE (Observations) = ', max(ae_index), 'SME: ', max(sme), labels[num])

     
    if glob.glob(ext_dir+dst_files[i]+'/'+'al.*'): 
        print (glob.glob(ext_dir+dst_files[i]+'/'+'al.*'))
        final_dates_ae, og_ae, og_al, og_au = get_og_ae(dst_files[i])
        ax.plot(final_dates_ae[:ind_cut], og_ae[:ind_cut], colors[i], linestyle='-',linewidth=1.0)#, label=labels[i])

        print('#----------------'+labels[num]+'----------------#')
        print('AE data and OpenGGCM: ', compute_r2(omni_dates,ae_index,final_dates_ae,og_ae), 'max(au) = ', max(og_au), 'min(al) = ', min(og_al), 'max(au) = ', max(og_ae))

    #if nn==0:
        #emp_dates, emp_ae = get_emp_ae(dst_files[i])
        #ax.plot(emp_dates, emp_ae, colors[i], linestyle='--',linewidth=1.0, label='Emp AE')

    plt.setp(ax.get_xticklabels(), visible=False)
    print('Starting to make plots...')

    plt.grid()
    draw_vertical(line,nn)

    #Show legend
    ax.yaxis.set_label_text("AE [nT]", fontweight="bold", fontsize=ytit)
    ax.yaxis.set_major_locator( MultipleLocator(600) )
    ax.yaxis.set_minor_locator( MultipleLocator(100) )
    '''
    l = np.arange(ymin,ymax, 600)
    x = np.arange(ymin,ymax, 600)
    ax.set_yticks(l)
    ax.set_yticklabels(x)
    '''
    plt.yticks(fontsize=8)
    ax.xaxis.set_minor_locator( AutoMinorLocator(6) )
    plt.setp(ax.get_xticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.setp(ax.get_yticklabels(), visible=True, fontsize=ytic, fontweight="bold")
    plt.legend(loc=1,prop={'size':leg_size},ncol=1)
    


    num=num+1



    #---------------------------------------------------------------------------------
    #Saving data for Metric analysis
    #---------------------------------------------------------------------------------
    # construct an interpolator

    # interpolated to same timesteps as CDAWeb data


    #############################################################################################################################################################
    #############################################################################################################################################################   

    '''
    #This commented chunk is not that good. Use the chunk below to save data for DTW.
    #Interpolating the OMNI onto OpenGGCM: This works. Do not change
    print (final_dates_dst[0],final_dates_dst[-1],omni_dates[0],omni_dates[-1])
    print ('Before interpolation:  ', len(og_delb_rcm), len(final_dates_dst), len(dst), len(omni_dates))

    omni_dst_int = interpol(dst, omni_dates, final_dates_dst) 
    print ('After interpolation:  ', len(og_dst_int), len(final_dates_dst), len(omni_dates))
    #og_ae_int = interpol(og_ae, final_dates_ae, omni_dates)

    og_ref = np.zeros(len(final_dates_dst))
    
    fig1 = plt.figure(figsize=(2, 10), dpi=80)
    plt.plot(final_dates_dst,og_dst_int, '--k', label='OMNI interpolated to OG')
    #plt.plot(omni_dates,dst, '--k', label='OMNI')
    plt.plot(final_dates_dst,og_delb_rcm, '--r', label='OG')
    plt.legend()
    plt.show()
    
    ##################################
    #Saving interpolated Dst indices
    ################################## 
    data_dst = [final_dates_dst, omni_dst_int, og_delb_rcm, og_ref]
    with open(dsv[i]+'/t_dst.pkl', 'wb') as outfile:
        print ('Saving .pkl @',dsv[i]+'/t_dst.pkl')
        pickle.dump(data_dst, outfile, pickle.HIGHEST_PROTOCOL) 

    print ('Interpolated Dst written onto binary files successfully.')
    
    ##################################
    #Saving interpolated AE indices
    ##################################
    og_ae_int = interpol(ae_index, omni_dates, final_dates_ae)

    fig2 = plt.figure(figsize=(2, 10), dpi=80)
    plt.plot(final_dates_ae,og_ae_int, '--k', label='AE interpolated to OMNI')
    plt.plot(omni_dates,ae_index, '--k', label='OMNI')
    plt.legend()
    plt.show()

    start_date=omni_dates[0]#datetime.strptime(date_beg, "%Y-%m-%dT%H:%M:%S")
    end_date=omni_dates[-1]#datetime.strptime(date_end, "%Y-%m-%dT%H:%M:%S")  
       
    ind_beg = omni_dates.index(start_date) #np.where(start_date, final_dates)[0]#
    ind_end = omni_dates.index(end_date)

    data_ae = [omni_dates[ind_beg:ind_end], ae_index[ind_beg:ind_end], og_ae_int]

    with open(dsv[i]+'/t_ae.pkl', 'wb') as outfile:
        print ('Saving .pkl @',dsv[i]+'/t_ae.pkl')
        pickle.dump(data_ae, outfile, pickle.HIGHEST_PROTOCOL)  

    print ('Interpolated AE written onto binary files successfully.')
    '''


    ######################################################
    #Uncomment if want to save interpolated indices
    ######################################################
    
    #Interpolating the OpenGGCM onto OMNI (FINAL!!!)
    '''
    
    if glob.glob(ext_dir+dst_files[i]+'/'+'dst*'):
        og_dst_int = interpol(og_delb_rcm, final_dates_dst, omni_dates) 
        og_ref = np.zeros(len(omni_dates)) #Changed from omni_dates on Jul 15, 2017
        start_date, end_date = closest_time(final_dates_dst,omni_dates)
        ind_beg = omni_dates.index(start_date); ind_end = omni_dates.index(end_date)    
        data_dst = [omni_dates[ind_beg:ind_end], dst[ind_beg:ind_end], og_dst_int[ind_beg:ind_end], og_ref[ind_beg:ind_end]]
    
    print(omni_dates[ind_beg:ind_beg+2],len(omni_dates[ind_beg:ind_end]))
    
        with open(dsv[i]+'/t_dst.pkl', 'wb') as outfile:
            print ('Saving .pkl @',dsv[i]+'/t_dst.pkl')
            pickle.dump(data_dst, outfile, pickle.HIGHEST_PROTOCOL)
    
    if glob.glob(ext_dir+dst_files[i]+'/'+'ae*'):
        og_ae_int = interpol(og_ae, final_dates_ae, omni_dates)
        ae_ref = np.zeros(len(omni_dates))
        start_date, end_date = closest_time(final_dates_ae,omni_dates)    
        ind_beg = omni_dates.index(start_date); ind_end = omni_dates.index(end_date)
        data_ae = [omni_dates[ind_beg:ind_end], ae_index[ind_beg:ind_end], og_ae_int[ind_beg:ind_end], ae_ref[ind_beg:ind_end]]
        with open(dsv[i]+'/t_ae.pkl', 'wb') as outfile:
            print ('Saving .pkl @',dsv[i]+'/t_ae.pkl')
            pickle.dump(data_ae, outfile, pickle.HIGHEST_PROTOCOL)  

    print ('Interpolated data written onto binary files successfully.')
    '''
    
    
    


    #############################################################################################################################################################
    #############################################################################################################################################################


print('Starting to make plots...')
plt.setp(ax.get_xticklabels(), visible=True)

# x axis format
myFmt = mdates.DateFormatter('%m-%d %H:%M') #%d-%m-%y')#
ax.xaxis.set_major_formatter(myFmt)
#rotating labels
fig.autofmt_xdate() 

ax.xaxis.set_minor_locator( AutoMinorLocator(6) )
#ax.xaxis.set_minor_locator( MultipleLocator(10) )
ax.set_xlim((start_date, end_date))
plt.setp(ax.get_xticklabels(), fontsize=ytic, visible=True)
    
#plt.suptitle(suptit)  
ax.xaxis.set_label_text("$\\bf{Date}$", fontsize=ytit, fontweight="bold")
output_directory = sys_dir+'/'+event[nn]
plt.text(0.95, -0.5, year_arr[nn], fontsize=12, fontweight="bold", ha='left', va='top', transform=ax.transAxes)
'''
print (ae_max, 'AE baseline average error wrt Cone:', 100*(ae_max[3]-ae_max[2])/ae_max[0]) #0: data 1: OG_OMNI 2: OG_Cone 3: OG_Spr
print (t_ae_max)
print (dst_min, 'Dst baseline average error wrt Cone:', 100*(dst_min[3]-dst_min[2])/dst_min[0])
print (t_dst_min)

print (ae_max, 'AE baseline average error wrt ref:', 100*(ae_max[3]-ae_max[1])/ae_max[0]) #0: data 1: OG_OMNI 2: OG_Cone 3: OG_Spr
print (dst_min, 'Dst baseline average error wrt ref:', 100*(dst_min[3]-dst_min[1])/dst_min[0])

print ('Time shift between Spheromak and data (AE)', (t_ae_max[3] - t_ae_max[0]).seconds/3600) #gives time difference in hours
print ('Time shift between Spheromak and data (Dst)', (t_dst_min[3] - t_dst_min[0]).seconds/3600) #gives time difference in hours
'''

fig.align_ylabels()
plt.tight_layout()
 
#plt.savefig(os.path.join(output_directory, png_name))   #'dps'+cone_name+'_'+fr_name  
if save: 
    if cut_dst: 
        plt.savefig(os.path.join(output_directory, 'cut_'+png_name+'.pdf'), dpi=300) 
    else:
        plt.savefig(os.path.join(output_directory, 'uncut_'+png_name+'.pdf'), dpi=300) 

else: plt.show() 
#plt.close(fig)