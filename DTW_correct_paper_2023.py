import os 
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DayLocator, HourLocator, DateFormatter
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator, AutoLocator
import scipy
from scipy.signal import savgol_filter
#%matplotlib inline
from scipy.ndimage.filters import uniform_filter1d
import cPickle as pickle #for loading data from npz

# real_stand and modeled_stand are the sequences being compared
# This is the only part that needs to be changed, based on the user's time series

sys_dir = './'

#Relevant data files
#Event1 is associated with the CME event on 20120712.
#Event2 is associated with the CME event on 20170904
event=['Event1','Event2']
year_arr=['2012','2017']


nn=1 #Change this while changing event

if nn==0:  #Also check line 125
        sw = 'Event1-obs'; spr = 'Event1-euh'
else: 
	sw = 'Event2-obs'; spr = 'Event2-euh'
	
euh_run = sw #Change this depending on which run you are analyzing
save=np.bool(True)
put_window=np.bool(True)

year = year_arr[nn]; loc = 'Earth'
if nn==0 and euh_run == spr:
	event_lab = 'Event1-euh'
elif nn==0 and euh_run == sw:
	event_lab = 'Event1-obs'
elif nn==1 and euh_run == sw:
	event_lab = 'Event2-obs'
else:
	event_lab = 'Event2-euh'
	

#real_stand=[1, 4, 7, 8, 10, 9, 6, 5, 2, 3, 3, 3]
#model_stand=[1, 1, 1, 4, 7, 8, 10, 9, 6, 5, 2, 3]
labsize=24; ticsize=16;

#The following locations are for modelled output
out_dir = sys_dir + event[nn]  + '/' + euh_run + '/' #'/home/u0141347/kul/paper_ongoing/Paper2_OpenGGCM/' #
outfile_dst = sys_dir + event[nn]  + '/' + euh_run + '/' + 't_dst.pkl' #'20120712_spr1_high_cad.pkl' #
outfile_ae = sys_dir + event[nn]  + '/' + euh_run + '/' + 't_ae.pkl'
index='dst' #Change this index while performing DTW on Dst or AE. Option: 'ae' or 'dst'

with open(outfile_dst, 'rb') as infile:
    infi_dst = pickle.load(infile)

    
with open(outfile_ae, 'rb') as infile:
    infi_ae = pickle.load(infile)	

#I was saving the modelled and observed dst in the same file after interpolating in the plotting files 
# infi_dst[dates, observed Dst, modelled Dst, Averaged observed Dst]   		
if index=='ae':
	real_stand = infi_ae[1]  
	model_stand = infi_ae[2]  
	avg_real_stand = infi_ae[3]
else:
	real_stand = infi_dst[1] 
	model_stand = infi_dst[2] 
	avg_real_stand = infi_dst[3]


'''
beg_ind = 0; end_ind = 0
for ii in range(1,len(model_stand)-1):
	print (isinstance(model_stand[ii], float), type(model_stand[ii]))
	print( model_stand[ii].isnumeric() )
	#if model_stand[ii]!=float(nan):# and isinstance(model_stand[ii-1], (float))==False:
	#beg_ind = np.where(model_stand.tolist().isnumeric())[0][0]
	if isinstance(model_stand[ii], (float))==True and isinstance(model_stand[ii+1], (float))==False:
		end_ind = ii
model_stand=model_stand[beg_ind:end_ind]
print (model_stand)
'''

#Smoothing both time series similarly
'''
Syntax:
savgol_filter(x, window_length, polyorder, deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
Parameter interpretation:
Low window length means the smooth profile is closer to the real data. 
High window length means the smooth profile is actually smoothed/averaged out. 
(Value of 7-11 works fine in our case)

yhat = savgol_filter(real_stand, 17, 2) worked well for 20170904_lowcad_1 dst
'''
yhat = savgol_filter(real_stand, 11, 2) #5 instead of 21 
plt.plot(real_stand, 'b-', label = 'Data', linewidth=1)
plt.plot(yhat, 'b--', label = 'Smoothed', linewidth=1)

xhat = savgol_filter(model_stand, 11, 2) #5 instead of 21
plt.plot(model_stand, 'r-' ,label='EUHFORIA+OpenGGCM', linewidth=3)
plt.plot(xhat, 'r--', label = 'Smoothed', linewidth=1)
plt.legend()
plt.xlabel('Time elements')
plt.ylabel('Amplitude')

if save: 
	if put_window:
            plt.savefig(out_dir + 'avg_win_'+ index + '_20120712_smooth.png', dpi=300, bbox_inches='tight')
	else:
	    plt.savefig(out_dir + index + '_20120712_smooth.png', dpi=300, bbox_inches='tight') #'/home/u0141347/kul/paper_ongoing/Paper2_OpenGGCM/dtw20120712/ae_20120712_smooth.png'
		
	
plt.tight_layout()
plt.close()
#plt.show()

real_stand = yhat
###########################################################
# Uncomment if want to perform DTW with the reference model
model_stand = xhat #xhat is used for OpenGGCM modelled; avg_real_stand is used for Dst=0 case.
###########################################################

cost = np.zeros((len(real_stand), len(model_stand))) # euclidean distance
DTW = np.ones((len(real_stand)+1, len(model_stand)+1)) # DTW cost array # Make the array 1 row and column larger than it has to be to allow for infty trick
DTW = DTW*np.infty # This is faster than iterating over the full matrix and replacing each value with infty

# The below-mentioned "w" is a time window of 300 elements which correspond to ~2days in real time
# This means that, for my purposes, I only permit the "matching" of the elements within a window of ~2days and not more than that. 
# If, in your case, you dont want to apply this window, please comment the following line and uncomment the line 
# "for j in range(1, len(model_stand)+1):" in the last for-loop.

w = np.max([24, abs(len(real_stand)-len(model_stand))]) #constraint the alignment at max 2days (~288 elements)
#30*5 = 150 minutes
 
DTW[0,0] = cost[0,0] # initialization of DTW
        
########## calculate all the other elements #############
for i in range(1, len(real_stand)+1):
    if put_window: 
        jmin=np.max([1, i-w]); jmax=np.min([len(model_stand)+1, i+w])
    else:
        jmin=1; jmax=len(model_stand)+1

    for j in range(jmin, jmax):   #### applying also a small punishment to the non-diagonal elements
        cost[i-1,j-1] = abs(model_stand[j-1]-real_stand[i-1]) #**2
        DTW[i,j] = min(DTW[i-1, j-1], DTW[i-1, j], DTW[i, j-1]) + cost[i-1,j-1] #is it the quantity that gives us the min cost to warp our timeseries?
        
DTW = DTW[1:,1:] # Remove the artificial shell of infty around the DTW matrix

print_cost = DTW[len(real_stand)-1, len(model_stand)-1]
print("The DTW cost is ", print_cost)

print(DTW)

# This is the cummulative cost matrix
def distance_distance_plot(DTW):
    im = plt.imshow(DTW, interpolation='nearest', cmap='Greens') 
    #plt.gca().invert_yaxis()
    plt.xlabel("sequence 1")
    plt.ylabel("sequence 2")
    plt.grid()
    plt.colorbar()
    #plt.clim(0,45)

###########################
# Figure 1
###########################
distance_distance_plot(DTW)


###### Backtracking technique to reveal the warpping path ############

path = [[len(model_stand)-1, len(real_stand)-1]]
i = len(real_stand)-1
j = len(model_stand)-1
while i>0 or j>0:
    if i==0:
        j = j - 1
    elif j==0:
        i = i - 1
    else:
        if DTW[i-1, j] == min(DTW[i-1, j-1], DTW[i-1, j], DTW[i, j-1]):
            i = i - 1
        elif DTW[i, j-1] == min(DTW[i-1, j-1], DTW[i-1, j], DTW[i, j-1]):
            j = j-1
        else:
            i = i - 1
            j = j - 1
    path.append([j, i])
path.append([0,0])    

path_x = [point[0] for point in path]
path_y = [point[1] for point in path]


distance_distance_plot(DTW)
plt.plot(path_x, path_y, linewidth=5);

# This is to reverse
def path_DTW(model_stand, real_stand, DTW, cost):
    path = [[len(model_stand)-1, len(real_stand)-1]]
    DTW_new = 0
    i = len(real_stand)-1
    j = len(model_stand)-1
    while i>0 or j>0:
        if i==0:
            j = j - 1
        elif j==0:
            i = i - 1
        else:
            if DTW[i-1, j] == min(DTW[i-1, j-1], DTW[i-1, j], DTW[i, j-1]):
                i = i - 1
            elif DTW[i, j-1] == min(DTW[i-1, j-1], DTW[i-1, j], DTW[i, j-1]):
                j = j-1
            else:
                i = i - 1
                j= j- 1
        path.append([j, i])
    path.append([0,0])
    for [real_stand, model_stand] in path:
        DTW_new = DTW_new + cost[model_stand, real_stand]
    return path, DTW_new

fig, ax = plt.subplots(figsize=(25, 10))
plt.plot(real_stand, 'b-', label = 'Observed data', linewidth=3)
if euh_run == spr: mod_label = 'EUHFORIA+OpenGGCM'
else: mod_label = 'Observations+OpenGGCM'
plt.plot(model_stand, 'r-' ,label=mod_label, linewidth=3)
plt.plot(avg_real_stand, 'k-' ,label='Reference model', linewidth=3)

paths = path_DTW(model_stand, real_stand, DTW, cost)[0]
for [map_x, map_y] in paths:
    #print map_x, x[map_x], ":", map_y, y[map_y]    
    plt.plot([map_x, map_y], [model_stand[map_x], real_stand[map_y]], 'g', linewidth=0.5)
    
import datetime
plt.xlim(0,len(model_stand))
plt.xticks(fontsize=ticsize)
plt.yticks(fontsize=ticsize)

plt.xlabel('Time elements', fontsize=labsize,fontweight='bold')
if (index == 'ae'):
	plt.ylabel(r'AE (nT)', fontsize=labsize,fontweight='bold')
if (index == 'dst'):
	plt.ylabel(r'Dst (nT)', fontsize=labsize,fontweight='bold')
	
plt.suptitle(event_lab,fontsize=labsize+8,fontweight='bold')
fig.subplots_adjust(top=0.94)
plt.legend(fontsize=labsize+1)

if save: 
	if put_window: 
	    plt.savefig(out_dir+'avg_win_'+index+'_' + euh_run + '.png', dpi=300, bbox_inches='tight')
	else: 
	    plt.savefig(out_dir+index+'_' + euh_run + '.png', dpi=300, bbox_inches='tight')
quit()
###########################
# Figure 2
###########################

print("======= see the alignments ========")
delta_t=[]; delta_amp=[]
for [map_x, map_y] in paths:
    #print(map_x,model_stand[map_x], ":", map_y, real_stand[map_y]) # x-coords, by subtracting them I estimate difference in time
    
    Time_diff = (map_y-map_x) # observed-modeled in days --- multiply by  10.0/1440.0 to convert to days or 10.0/24.0 to convert to hours if you have a 10-min resolution --> This was for Lia
    Amp_diff = real_stand[map_y] - model_stand[map_x]   #observed-modeled in nT
    delta_t.append(Time_diff*5.0); #Time_diff is in 5 minutes cadence --> Multiply by 5.0/60 to convert to hours
    delta_amp.append(Amp_diff)
    
print (type(map_y),map_y)
print ('range of delta_t', min(delta_t), max(delta_t))
print ('range of delta_amp', min(delta_amp), max(delta_amp))
#------------------ Plot histograms for Time_diff and Amp_diff ---------------------#

##########################################################################################################################
##########################################################################################################################
# Plotting the offset histograms in time and amplitude
##########################################################################################################################
##########################################################################################################################
from matplotlib.colors import Normalize
from matplotlib import cm

#bins_num = np.linspace(int(min(delta_amp))-1, int(max(delta_amp))+1, 25) #-60, 20, 9)
bins_num = np.arange(int(min(delta_amp)), int(max(delta_amp)), 10)
print ('bin size amps: ',int(min(delta_amp)),int(max(delta_amp)), bins_num)
wind_weights1 = np.ones_like(delta_amp)/float(len(delta_amp)) #Return an array of ones with the same shape and type as a given array. This implies all the amplitudes are given the same weight

data=[delta_amp]
weights_all=[wind_weights1]
colors=['blue']
labels=['WIND']

###########################
# Figure 3
###########################
fig, ax = plt.subplots(figsize=(12, 10))

n, bins, patches = plt.hist(data, bins_num, weights=weights_all, histtype='bar', edgecolor='black', color=colors, label=labels)
norm = Normalize(min(delta_t), max(delta_t)) #(-60, 20) #Normalization is not being applied
ticks = [(patch._x0 + patch._x1)/2 for patch in patches]
if (index=='ae'):
	ax.set_xlabel(r'$\Delta$AE$_{Observed- Predicted}$ (nT)', fontsize=labsize,fontweight='bold')
if (index=='dst'):
	ax.set_xlabel(r'$\Delta$Dst$_{Observed- Predicted}$ (nT)', fontsize=labsize,fontweight='bold')
	ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax.set_ylabel('frequency', fontsize=labsize,fontweight='bold')
ax.set_title('', fontsize=15)
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
#plt.xlim(-100,100)
plt.xticks(ticks,fontsize=ticsize)
plt.yticks(fontsize=ticsize)

if save: 
	print ('saving amp hists')
	if put_window:
	    plt.savefig(out_dir + 'win_hist_ampdiff_'+index+euh_run+'.png', dpi=300, bbox_inches='tight')    
	else: 
	    plt.savefig(out_dir + 'hist_ampdiff_'+index+euh_run+'.png', dpi=300, bbox_inches='tight') #/home/u0141347/kul/paper_ongoing/Paper2_OpenGGCM/dtw20120712/
#plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0, hspace=0)


#-----------------------------------------------------------------------------------------------

#bins_num = np.linspace(min(delta_t), max(delta_t),11) #-1, 1, 17) #this is when i have maximum +/- 1day difference
bins_num = np.arange(int(min(delta_t)), int(max(delta_t)), 10) #the integral number is the interval size
print ('bin size times: ',int(min(delta_t)),int(max(delta_t)), bins_num)
wind_weights1 = np.ones_like(delta_t)/float(len(delta_t)) #Return an array of ones with the same shape and type as a given array.

data=[delta_t]
weights_all=[wind_weights1]
colors=['red']
labels=['WIND']

###########################
# Figure 4
###########################
fig, ax = plt.subplots(figsize=(12, 10))

n, bins, patches=plt.hist(data, bins_num, weights=weights_all, histtype='bar', edgecolor='black', color=colors,label=labels)
norm = Normalize(min(delta_t), max(delta_t)) #(-1,1) #Normalization is not being applied
ticks = [(patch._x0 + patch._x1)/2 for patch in patches]
#ticklabels = [i for i in range(len(bins_num))]
plt.xticks(ticks,fontsize=ticsize)#, ticklabels)
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel('$\Delta$t$_{Observed - Predicted}$ (minutes)', fontsize=labsize,fontweight='bold')
ax.set_ylabel('frequency', fontsize=labsize,fontweight='bold')
ax.set_title('', fontsize=15)
plt.yticks(fontsize=ticsize)

if save:        
	if put_window: 
	    print ('saving time hists')
	    plt.savefig(out_dir+'win_hist_tdiff_'+index+euh_run+'.png', dpi=300, bbox_inches='tight')
	else: 
	    plt.savefig(out_dir+'hist_tdiff_'+index+euh_run+'.png', dpi=300, bbox_inches='tight') #'/home/u0141347/kul/paper_ongoing/Paper2_OpenGGCM/dtw20120712/
	    
#plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0, hspace=0)
else: plt.show()
#plt.savefig('Histos/CR2198_TimeDiff.png', dpi=100, bbox_inches='tight')
#plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0, hspace=0)
