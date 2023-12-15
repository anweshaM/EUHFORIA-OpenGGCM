import glob 
import datetime
import os.path
import ConfigParser

import numpy as np
import os
import sys
print sys.argv

from datetime import datetime
# https://openggcm.sr.unh.edu/?n=Main.Inputs

# =================================================================================================
# -----------------------------------------------------------------------------------------------------
# Program that processes a EUHFORIA 1.0 data file in GSE to be fed in OpenGGCM 3.0
# Developer: Camilla Scolini, Anwesha Maharana
# Last modified: May 2022
# Python vers. 2.7-18
# -----------------------------------------------------------------------------------------------------
event='Event1'
euh_dir = './'+event+'/Event1-euh/' #/Event2-euh

# file you want to process
# column entries:
# date btot bxgse[nT] bygse[nT] bzgse[nT] vtot vxgse[km/s] vygse[km/s] vzgse[km/s] n[1/cm^3] temp[K] P[pPa] x[gse] y[gse] z[gse]
input_file_euh = euh_dir + 'euh_10min_gse.txt'

# output files
output_file_1 = euh_dir +'euh.np'
output_file_2 = euh_dir +'euh.pp'
output_file_3 = euh_dir +'euh.temp'
output_file_4 = euh_dir +'euh.bxgse'
output_file_5 = euh_dir +'euh.bygse'
output_file_6 = euh_dir +'euh.bzgse'
output_file_7 = euh_dir +'euh.btot'
output_file_8 = euh_dir +'euh.vxgse'
output_file_9 = euh_dir +'euh.vygse'
output_file_10 = euh_dir +'euh.vzgse'
output_file_11 = euh_dir +'euh.vtot'
output_file_12 = euh_dir +'euh.xgse'
output_file_13 = euh_dir +'euh.ygse'
output_file_14 = euh_dir +'euh.zgse'

output_file_15 = euh_dir +'euh.year'
output_file_16 = euh_dir +'euh.month'
output_file_17 = euh_dir +'euh.day'
output_file_18 = euh_dir +'euh.hr'
output_file_19 = euh_dir +'euh.min'
output_file_20 = euh_dir +'euh.sec'

output_file_21 = euh_dir +'euh.rr'

#--------------------------------------------------------

euh_data = np.loadtxt(input_file_euh, skiprows=1, usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13])
euh_date = np.loadtxt(input_file_euh, dtype='str', skiprows=1, usecols=[0])

# assigning date to euhforia data
final_dates = []
for dates in euh_date:
	final_dates.append(datetime.strptime(str(dates), "%Y-%m-%dT%H:%M:%S") )

for i in range(len(euh_data[:,10])):
        euh_data[i,10] = 30.0 #R_E (unit of radius of Earth)
        euh_data[i,11] = 0.0
        euh_data[i,12] = 0.0
        euh_data[:,1] = 0.0 #At magnetopause boundary, Br = 0.0

euh_xgse = euh_data[:,10] 
euh_ygse = euh_data[:,11]
euh_zgse = euh_data[:,12]

euh_np = euh_data[:,8]*0.5 #In EUHFORIA density = total density, not proton [part/cc]
euh_rr = euh_np #[part/cc]
euh_pp = euh_data[:,9] #[pPa]
print ("Dumb: ",euh_pp[20:30]) 
euh_tp = euh_pp*10**(-12)/(euh_np*(10**6.0)*1.38065*10**(-23.0)) # = tp [K]

euh_vxgse = euh_data[:,5]
euh_vygse = euh_data[:,6]
euh_vzgse = euh_data[:,7]

euh_bxgse = euh_data[:,1]
euh_bygse = euh_data[:,2]
euh_bzgse = euh_data[:,3]

euh_btot = euh_data[:,0]
euh_vtot = euh_data[:,4]

euh_yr=np.empty(len(euh_xgse))
euh_mm=np.empty(len(euh_xgse))
euh_dd=np.empty(len(euh_xgse))
euh_hr=np.empty(len(euh_xgse))
euh_min=np.empty(len(euh_xgse))
euh_sec=np.empty(len(euh_xgse))

for i in range (0, len(final_dates)):
    euh_yr[i] = int(final_dates[i].strftime("%Y"))
    euh_mm[i] = int(final_dates[i].strftime("%m"))
    euh_dd[i] = int(final_dates[i].strftime("%d"))	
    euh_hr[i] = int(final_dates[i].strftime("%H"))
    euh_min[i] = int(final_dates[i].strftime("%M")) 
    euh_sec[i] = final_dates[i].strftime("%S") 

# =================================================================================================
# CREATING OUTPUT FILES
# =================================================================================================

print "Creating output files..."

file1 = open(output_file_1, "w")
file2 = open(output_file_2, "w")
file3 = open(output_file_3, "w")
file4 = open(output_file_4, "w")
file5 = open(output_file_5, "w")
file6 = open(output_file_6, "w")
file7 = open(output_file_7, "w")
file8 = open(output_file_8, "w")
file9 = open(output_file_9, "w")
file10 = open(output_file_10, "w")
file11 = open(output_file_11, "w")
file12 = open(output_file_12, "w")
file13 = open(output_file_13, "w")
file14 = open(output_file_14, "w")
file15 = open(output_file_15, "w")
file16 = open(output_file_16, "w")
file17 = open(output_file_17, "w")
file18 = open(output_file_18, "w")
file19 = open(output_file_19, "w")
file20 = open(output_file_20, "w")
file21 = open(output_file_21, "w")

for i in range(len(euh_yr)):
    file1.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.5f" % euh_np[i]) + "\n")
    file21.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.5f" % euh_rr[i]) + "\n")
    file2.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.8f" % euh_pp[i]) + "\n")
    file3.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.1f" % euh_tp[i]) + "\n")
    file4.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.6f" % euh_bxgse[i]) + "\n")
    file5.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.6f" % euh_bygse[i]) + "\n")
    file6.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.6f" % euh_bzgse[i]) + "\n")
    file7.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.6f" % euh_btot[i]) + "\n")
    file8.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.3f" % euh_vxgse[i]) + "\n")
    file9.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.3f" % euh_vygse[i]) + "\n")
    file10.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.3f" % euh_vzgse[i]) + "\n")
    file11.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.3f" % euh_vtot[i]) + "\n")
    file12.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.4f" % euh_xgse[i]) + "\n")
    file13.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.4f" % euh_ygse[i]) + "\n")
    file14.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.4f" % euh_zgse[i]) + "\n")
    file15.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%04d" % euh_yr[i]) + "\n")
    file16.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%02d" % euh_mm[i]) + "\n")
    file17.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%02d" % euh_dd[i]) + "\n")
    file18.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%02d" % euh_hr[i]) + "\n")
    file19.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%02d" % euh_min[i]) + "\n")
    file20.write(str("%.0f" % euh_yr[i]) + " " + str("%.0f" % euh_mm[i]) + " " + str("%.0f" % euh_dd[i]) + " " + str("%02d" % euh_hr[i]) + " " + str("%02d" % euh_min[i]) + " " + str("%02.3f" % euh_sec[i]) + " " + str("%.3f" % euh_sec[i]) + "\n")

file1.close()
file2.close()
file3.close()
file4.close()
file5.close()
file6.close()
file7.close()
file8.close()
file9.close()
file10.close()
file11.close()
file12.close()
file13.close()
file14.close()
file15.close()
file16.close()
file17.close()
file18.close()
file19.close()
file20.close()
file21.close()
