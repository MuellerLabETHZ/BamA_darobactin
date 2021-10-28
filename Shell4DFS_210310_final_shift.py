#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 17:44:22 2020

@author: rnoah
"""
import os
from shutil import copy
import numpy as np
import pickle
import matplotlib.pyplot as plt

from scipy.stats import linregress

from FD import FD
from BLC import BLC
from PF import PF
from Transformation import Trafo
from Denoising import denoising
from Autoalign import autoalign
from Aligntoprev import aligntoprev
from Peak_Preview import peak_preview
from Move_Interface import move_interface
from Mergeandcopy import mergeandcopy
from Shift_preview import shift_preview
from Density_plot_Final_shift import density_plot
from Optics import optics
from Pathway_rgb_input import pathway

from SplitByBoarders_final_shift import splitbyboarders
from SplitByBoarders_final_shift import assign_single_by_boarders

from Import_colors import import_rgb_colors

#%% User Input

#### Location of inputfiles
inputfolder = 'dummy'
inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_Daro/1ums_D+'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_Daro/3ums_D+'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_Daro/4_5ums_D+'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_Daro/05ums_D+'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_Daro/6ums_D+'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_Daro/07ums_D+'

# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_No_Daro/1ums_D-'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_No_Daro/3ums_D-'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_No_Daro/4_5ums_D-'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_No_Daro/05ums_D-'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_No_Daro/6ums_D-'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Aligned_No_Daro/07ums_D-'

# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Darobactin_Initial_Unfolding/Aligned/Daro'
# inputfolder = '/Users/rnoah/Desktop/DFS_Aligned/Darobactin_Initial_Unfolding/Aligned/NoDaro'
# Name of final pickle-file
picklename = 'Alignment'

#### Only use this if new data should be added to a previous experiment
# Add path to pickle-file of previous dataset to combine, else False
prev_pickle_dir = False
# Add path to folder containing files of previous dataset, else False
prev_files_dir = False

#### Optional steps for the procedure, choose with True or False
# First_Align=True -> Data is aligned for the first time and pickle is created
# First_Align=False -> Data has already been aligned and respective pickle is loaded
First_Align = False
# If True auto-alignment will be performed
do_auto_align = False
# If True, selected curves will be copied into a new subdirectory in inputfolder
saveselected = False
# If True, aligned curves can in the end be shifted in x by the user
adjust_global_shift = True
# If True, clustering and pathway analysis will be performed
perform_clustering = True
# If True, density plot will be generated
do_density_plot = True

#### General settings
# Defines which segment of a force-curve should be analyzed
Segment = 2
# Enter Force-CL to align curves in Contourlength space, else enter Force-Dist
Space = 'Force-CL'
#Threshold used for Transformation - Vdef vals < trafo_tresh are considered as noise
trafo_tresh = 20
# Last x percent for linear fit to perform Baselinecorrection - enter 0 to skip this step
percentforBLC = 20
#Threshold used to discard tiny peaks, e.g. in Baseline (minimal peak-width [nm])
PF_Thr = 7


#### Denoising settings for peak-finding
#Denoising parameters (for more details check pywt documentation):
dwt_iterations=3 #Levels of decomposition
dwt_thr=0.1 #thr = multiplier of denoising threshold - is multiplied with ppnm
dwt_level_factor=4 #level_factor = scaling of denoising threshold per level
dwt_wave = 'db4' #type of wavelet
dwt_thr_mode = 'soft' #choose soft or hard threshold

#### Alignment settings
# Define Bin-sizes for alignment of curves
binsize_coarse = 10 # Binsize for Cross-alignment -> Less precise but fast
binsize_fine = 1 # Binsize for final alignment -> more precise but slow
Allow_y_shift = False # If True, user can shift curves in y

#### Define Region of interest
roi = 350 #  Define region of interest for Density Plot and Auto-alignment (x-max and y_max)
bins_per_nm = 2 # Defines resolution of density plot

#### Nr of curves to display to optimize peak-finding
nr_of_example_curves = 5 #Cannot be smaller than files in dataset, suggested value: 5-10
prom = 40 #Start-prominence for peak-finding -> will be re-adjusted by user

#### Parameters for optics-clustering
# min_samples: min nr of neighboring points for a point to become core-point
# This is defined as sqrt(nr_of_curves) - can be altered in code below
op_y_scaling_factor=50 # 50 according to Yosh's algo. In order to 'normalize' vdef and TSS.
op_xi=.04 # minimum steepness on reachability plot to define cluster boundary
op_min_cluster_size=.02 # min samples per cluster as fraction of total samples
#%% Import Files
directory = os.path.abspath(inputfolder)
pickle_dir = os.path.join(directory, 'Pickles')


if First_Align == True:
    allfiles = os.listdir(directory)
    inputfiles = []
    for filename in allfiles:
        if filename.endswith('.txt') and not filename.startswith('Analysis'):
            inputfiles.append(filename)


    ND = {} # Dict that pairs filename with its directory

    cutoff1 = 40000
    cutoff2 = 40100
    
    ### Pairs every filename with its directory
    int_filenames = []
    for file in inputfiles:
        drct = os.path.join(directory, file)
        filename = drct.split('/')[-1].rstrip('.txt')
        ### Use this to remove parts of the dataset:
        # int_filenames.append(int(filename))
        # if int(filename) >= cutoff1 and int(filename) <= cutoff2:
        #     ND[filename] = drct
        ND[filename] = drct

    TotalNrOfCurves = len(inputfiles)
    # TotalNrOfCurves = len(ND)



    FD_Curves = []

    Progresscounter = 0
    print('Load Curves:')
    for index, name in enumerate(ND):

        ### Initialize Progresscounter
        print(name)
        Progresscounter += 1
        print(str(Progresscounter) + '/' + str(TotalNrOfCurves) + ' Imported')

        ### Extract TSS (x_raw) and VD (y_raw) from txt-file:
        FD_Curves.append({})
        FD_Curves[index]['name'] = name
        with open(ND[name], 'r') as f:
            x_raw = []
            y_raw = []
            for line in f:
                if line.startswith('# Sensitivity [nm/V]'):
                    sensitivity = float(line.split()[-1])
                elif line.startswith('# Spring Constant [N/m]'):
                    spring_const = float(line.split()[-1])
                if not line.startswith('#'):
                    x_raw.append(float(line.split(',')[0])) # TSS Imported as nm
                    y_raw.append(float(line.split(',')[1])*1000) # Vdef Convert nN to pN
        print('Sensitivity = ', sensitivity, 'Spring Constant = ', spring_const)

        x_raw = np.asarray(x_raw)
        y_raw = np.asarray(y_raw)
        ppnm = 10

        if percentforBLC != 0:
            x_BLC, y_BLC = BLC(x_raw, y_raw, percentforBLC)
        else:
            x_BLC = x_raw
            y_BLC = y_raw

        # DWT needs even-length samples: adjust if necessary:
        if len(x_raw)%2 != 0:
            x_raw = x_raw[:-1]
            y_raw = y_raw[:-1]
        y_DWT = denoising(y_BLC,
                          ppnm,
                          dwt_iterations,
                          dwt_thr,
                          dwt_level_factor,
                          dwt_wave,
                          dwt_thr_mode)

        ###------------------------------------------
        x_PF,y_PF = PF(x_BLC, y_DWT, ppnm, PF_Thr)
        x_Trafo, y_Trafo = Trafo(x_BLC, -y_BLC, trafo_tresh)
        FD_Curves[index]['x_raw'] = x_BLC
        FD_Curves[index]['y_raw'] = -y_BLC
        FD_Curves[index]['x_PF'] = x_PF
        FD_Curves[index]['y_PF'] = -y_PF
        FD_Curves[index]['x_Trafo'] = x_Trafo
        FD_Curves[index]['y_Trafo'] = y_Trafo

#%% Auto-ALignment
    if do_auto_align == True:
        if prev_pickle_dir != False:
            FD_Curves_prev = pickle.load(open(prev_pickle_dir,"rb"))
            FD_Curves = aligntoprev(FD_Curves, FD_Curves_prev, binsize_fine, Space, trafo_tresh)
        else:
            FD_Curves_prev = False
            FD_Curves = autoalign(FD_Curves, binsize_coarse, binsize_fine, Space, trafo_tresh)
        # Sort curves based on their global alignment scores
        # 1. Create list of all alignment scores
        all_scores = []
        for i in range(TotalNrOfCurves):
            all_scores.append(FD_Curves[i]['Glob_Align_Score'])
        # 2. Perform argsort
        sorted_ind = np.argsort(all_scores)
        # 3. Change indeces of FD_Curves accordingly
        FD_Curves = [FD_Curves[i] for i in sorted_ind]
    else:
        for i in range(TotalNrOfCurves):
            FD_Curves[i]['Auto_Shift'] = 0
        FD_Curves_prev = False

#%% Manual section
    #First window: User can adjust peak-finding parameters to respective data
    prom = peak_preview(FD_Curves, nr_of_example_curves, prom)
    #Second window: Curves are aligned by the user
    FD_Curves = move_interface(FD_Curves, FD_Curves_prev, Space, Allow_y_shift, prom, roi)

#%% Copy selected Files to a new folder and merge FD_Curves and FD_Curves_prev
    if prev_files_dir != False:
        FD_Curves = mergeandcopy(FD_Curves,
                                 FD_Curves_prev,
                                 prev_files_dir,
                                 directory,
                                 saveselected)

    # If data is not added from 2 directories, just copy selected curves
    else:
        # Copy selected Files to a new folder
        if saveselected == True:
            destdir = os.path.join(directory, 'Selected_Curves')
            os.makedirs(destdir, exist_ok=True)
            for i in range(TotalNrOfCurves):
                if FD_Curves[i]['Selected'] == True:
                    filename = FD_Curves[i]['name'] + '.txt'
                    filedir = os.path.join(directory, filename)
                    copy(filedir, destdir)

    # Store FD_Curves in new pickle
    if not os.path.exists(pickle_dir):
        os.makedirs(pickle_dir)
    pickle.dump(FD_Curves, open(os.path.join(pickle_dir, picklename),"wb"))

#%% Pickle-Load of pre-aligned data
else:
    FD_Curves = pickle.load(open(os.path.join(pickle_dir, picklename),"rb"))
    
TotalNrOfCurves =len(FD_Curves)

count=0
#Sensitivity and Spring Cons were not stored initially. Extract here:
for i in range(TotalNrOfCurves):
    if FD_Curves[i]['Selected'] == True:
        count+=1
        filedir = os.path.join(inputfolder, 'Selected_Curves', FD_Curves[i]['name']+'.txt')
        with open(filedir, 'r') as f:
            for line in f:
                if line.startswith('# Sensitivity [nm/V]'):
                    FD_Curves[i]['sensitivity'] = float(line.split()[-1])
                elif line.startswith('# Spring Constant [N/m]'):
                    FD_Curves[i]['spring_const'] = float(line.split()[-1])   
        

#%% For DFS
#%% Shift all curves to introduce global shift
if adjust_global_shift == True:
    shift_dict = {
        '1ums_D+': -198.18969,
        '3ums_D+': -178.99964,
        '4_5ums_D+': -162.10903,
        '05ums_D+': -148.41,
        '6ums_D+': -170.848,
        '07ums_D+': -216.385,
        '1ums_D-': -103.073,
        '3ums_D-': -69.5375,
        '4_5ums_D-': -131.51,
        '05ums_D-': -138.862,
        '6ums_D-': -167.316,
        '07ums_D-': -135.783
        }

    # Add 406aa to shift (146.16nm) to move first HP from zero to correct dist
    for key in shift_dict:
        shift_dict[key] += 146.16

    condition = inputfolder.split('/')[-1]
    global_shift = shift_dict[condition]

else:
    global_shift = 0

# #%% For initial unfolding
# global_shift = 74
# # global_shift= shift_preview(FD_Curves, Space)
# print('final global shift:', global_shift)
#%%

# Shift peak positions by global shift
for i in range(TotalNrOfCurves):
    if FD_Curves[i]['Selected'] == True:
        FD_Curves[i]['x_peaks'] = [a+global_shift for a in FD_Curves[i]['x_peaks']]

# Create new entries in FD-Curves with shifted curves
for i in range(TotalNrOfCurves):
    if FD_Curves[i]['Selected'] == True:
        auto_shift = FD_Curves[i]['Auto_Shift']
        x_shift = FD_Curves[i]['x_shift']
        y_shift = FD_Curves[i]['y_shift']
        total_x_shift = global_shift + auto_shift + x_shift
        if Space == 'Force-CL':
            curve_x = [a+total_x_shift for a in FD_Curves[i]['x_Trafo']]
            curve_y = [a+y_shift for a in FD_Curves[i]['y_Trafo']]
        else:
            curve_x = [a+total_x_shift for a in FD_Curves[i]['x_raw']]
            curve_y = [a+y_shift for a in FD_Curves[i]['y_raw']]
        FD_Curves[i]['curve_x_shifted'] = curve_x
        FD_Curves[i]['curve_y_shifted'] = curve_y


#%% Import color scheme as rgb values
rgb_file_dir = '/Users/rnoah/Documents/Manuscripts/BamA_DFS_MS/Figures/Color_scheme.txt'
color_scheme, name_scheme = import_rgb_colors(rgb_file_dir)

#%% For DFS

# Extracts Infos from Title
from_title = inputfolder.split('/')[-1]

unfolding_speed = from_title.split('u')[0]
if unfolding_speed == '4_5':
    unfolding_speed = 4.5
elif unfolding_speed == '05':
    unfolding_speed = 0.5
elif unfolding_speed == '07':
    unfolding_speed = 0.7
else:
    unfolding_speed = int(unfolding_speed)
daro = from_title[-1]

Dir_Title = str(unfolding_speed) + ' um_s'+ ' Unfolding Speed, ' + daro + 'Darobactin'
if daro == '-':
    daro = '–'
Title = str(unfolding_speed) + ' $\mu$'+'$s^{-1}$'+ ' Pulling Speed, ' + daro + 'Darobactin'


plot_dir = os.path.join('/Users/rnoah/Desktop/DFS_Aligned/Python_Plot_Autosave_210720', Dir_Title)
font_size = 18
x_size = 7.4
y_size = x_size*0.8
fig_size = [x_size, y_size]
DPI = 300

#%% For Intial Unfolding
# Title = 'initial_' + inputfolder.split('/')[-1]
# plot_dir = os.path.join('/Users/rnoah/Desktop/DFS_Aligned/Python_Plot_Autosave_210720', Title)

#%%

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

#%% Clustering and Pathway analysis
# Pool all detected peaks of selected curves for clustering
Nr_of_Selected = 0
if perform_clustering == True:
    all_peaks_x = []
    all_peaks_y = []
    for i in range(TotalNrOfCurves):
        if FD_Curves[i]['Selected']:
            Nr_of_Selected += 1
            all_peaks_x = np.append(all_peaks_x, FD_Curves[i]['x_peaks'])
            all_peaks_y = np.append(all_peaks_y, FD_Curves[i]['y_peaks'])

#     # Do the optics clustering and plotting
#     op_min_samples=int(np.sqrt(TotalNrOfCurves))
#     clustering = optics(all_peaks_x,
#                 all_peaks_y,
#                 op_y_scaling_factor,
#                 op_min_samples,
#                 op_xi,
#                 op_min_cluster_size,
#                 Space)

#     # Run the pathway analysis on the clustered data
#     pathway(clustering, FD_Curves, Space)

#%% Use this to clusterpeaks by fixed boarders
boarderclusters, noise_points, boarders = splitbyboarders(all_peaks_x, all_peaks_y, 0)
fig, ax = plt.subplots(figsize=fig_size)
# plt.plot(noise_points['x'], noise_points['y'], '.', markerfacecolor="None", markeredgecolor='lightgrey', alpha=200/Nr_of_Selected)
# [plt.vlines(i, 0, 300) for i in boarders]
for index, cluster in enumerate(boarderclusters):
    fig.suptitle(Title, size=20)
    ax.plot(boarderclusters[cluster]['x'], boarderclusters[cluster]['y'],  '.', markersize=4, alpha=100/Nr_of_Selected, color=color_scheme[index])
    ax.set_xlabel('Contour Length [nm]', size=font_size, labelpad=2)
    ax.set_ylabel('Force [pN]', size=font_size)
    ax.set_xlim(xmin = 0, xmax=300)
    ax.set_ylim(ymin = 0, ymax=450)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(os.path.join(plot_dir, 'clustering.png'), dpi=DPI) # For paper, just add required dpi=x



#%% Do a pathway analysis and probability of peak-occurences on the boarderclusters
pathway(boarderclusters, FD_Curves, Space, color_scheme, plot_dir, font_size, fig_size, DPI, Title)
fig, ax = plt.subplots(figsize=fig_size)
fig.suptitle(Title, size=20)
ax.set_xlabel('Peak-Position', size=font_size)
ax.set_ylabel('Frequency [a.u.]', size=font_size)
ax.set_ylim(ymin = 0, ymax=1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
for cluster in boarderclusters:
    ax.bar(cluster, len(boarderclusters[cluster]['x'])/Nr_of_Selected, color=color_scheme[cluster], edgecolor='k', linewidth=0.5)

# extract avg and stdev of clusters in x and y
    clusterdata = np.zeros((len(boarderclusters), 5))
    for cluster in boarderclusters:
        x_avg = np.average(boarderclusters[cluster]['x'])
        x_std = np.std(boarderclusters[cluster]['x'])
        y_avg = np.average(boarderclusters[cluster]['y'])
        y_std = np.std(boarderclusters[cluster]['y'])
        clusterdata[cluster,:] = cluster, x_avg, x_std, y_avg, y_std
plt.xticks(range(len(boarderclusters)), name_scheme)
plt.xticks(fontsize=15)
plt.yticks(fontsize=font_size)
plt.savefig(os.path.join(plot_dir, 'Peak_Prob.png'), dpi=DPI)

#%% Write Forces and LRs to txt file
f = open("/Users/rnoah/Desktop/DFS_Aligned/Python_Plot_Autosave/initial_NoDaro/All_Forces.txt", "w")
for clust in boarderclusters:
    f.write('segment ' + str(clust) + ', ')
    for val in boarderclusters[clust]['y']:
        f.write(str(val) + ', ')
    f.write('\n')
f.close()


#%% Do a density plot of all selected curves
if do_density_plot == True:
    density_plot(FD_Curves, roi, bins_per_nm, Space, Title, color_scheme, name_scheme, clusterdata, plot_dir, font_size, DPI)



#     #Define boarders between clusters -> Midpoint between cluster x avg + stds
#     #First and last point are mirrored from 2nd from start and end
#     boarders = []
#     for i in range(len(clusterdata)-1):
#         left_avg = clusterdata[i,1]
#         left_std = clusterdata[i,2]
#         right_avg = clusterdata[i+1,1]
#         right_std = clusterdata[i+1,2]
#         boarder = ((right_avg-right_std)+(left_avg+left_std))/2
#         boarders.append(boarder)
#     last_distance = clusterdata[-1,1] - boarders[-1]
#     first_distance = boarders[0]-clusterdata[0,1]
#     boarders.append(clusterdata[-1,1]+last_distance)
#     boarders.insert(0,clusterdata[0,1]-first_distance)


#%%
# #### Loading Rate Stuff

# # Anal_Dir = os.path.join(inputfolder, 'Analysis_'+str(unfolding_speed)+'ums_5percNoise.txt')
# Anal_Dir = os.path.join(inputfolder, 'Dummy.txt')
# with open(Anal_Dir, 'w+') as f:
#     f.write('Filename, Cluster, Distance[nm], Force[pN], Loading-Rate[pN/s], Sensitivity [nm/V], Spring Constant [N/m]\n')
#     sampling_rate = unfolding_speed*10000
    
#     # Update i_peaks -> not updated in move_interface -> needs fix in the future
#     for entry in FD_Curves:
#         i_peaks = []
#         #1st remove duplicates from x- and y-peaks
#         entry['x_peaks'] = list(dict.fromkeys(entry['x_peaks']))
#         entry['y_peaks'] = list(dict.fromkeys(entry['y_peaks']))

#         #Find closest instead of exact match (exact match was buggy):
#         #If peak_index is larger than previous, add it to i_peaks
#         #If smaller, remove from x- and y-peaks
#         for peak_index, peakval in enumerate(entry['y_peaks']):
#             closest = 100000000
#             for index, raw_val in enumerate(entry['y_raw']):
#                 if abs(peakval-raw_val) < closest:
#                     closest = abs(peakval-raw_val)
#                     closest_ind = index
#             if len(i_peaks) == 0:
#                 i_peaks.append(closest_ind)
#             elif closest_ind > i_peaks[-1]:
#                 i_peaks.append(closest_ind)
#             else:
#                 entry['x_peaks'].pop(peak_index)
#                 entry['y_peaks'].pop(peak_index)
#         entry['i_peaks'] = i_peaks
        
        
#         # plt.figure()
#         # plt.title(entry['name'])
#         # plt.plot(entry['x_raw'], entry['y_raw'])
#         # plt.plot(entry['x_peaks'], entry['y_peaks'], '.')
#         # plt.plot(entry['x_raw'][i_peaks], entry['y_raw'][i_peaks], 'x')
#         # print(entry['name'], len(entry['i_peaks']), len(entry['x_peaks']), len(entry['y_peaks']))
    
#     # Extract Force and Time of each Curve
#     percent_for_fit = 20
#     min_length = 50 # minimal points for reliable fit
    
#     for entry in FD_Curves:
#         if entry['Selected'] == True:
#             force = entry['y_raw']
#             time = np.arange(0,len(force)*(1/sampling_rate), 1/sampling_rate)
#             peak_inds = entry['i_peaks']
#             sensitivity = str(entry['sensitivity'])
#             spring_const = str(entry['spring_const'])
            
#             # plt.figure()
#             # plt.plot(time, entry['y_raw'], color='k', alpha=0.5)
#             # plt.plot(time, entry['y_PF'], color='g', alpha=0.5)
#             # plt.plot(time[peak_inds], force[peak_inds], 'x', color='r')
#             # plt.title(entry['name'])
            
#             #Write to txt-file
        
#         # Find Valleys between and Ascends before Peaks
#             for i in range(0, len(peak_inds), 1):
#                 #Find Ascend of first peak:
#                 if i == 0:
#                     min_ind = next((i for i, x in enumerate(force) if x > 0), None)
#                     # plt.plot(time[min_ind], force[min_ind], 'x', color='green')
#                 else:
#                     #Define Interval between peak(i-1) and peak(i)
#                     interval = range(peak_inds[i-1],peak_inds[i],1)
#                     #Find index of lowest point
#                     # print('interval_Length:', len(interval), peak_inds[i-1], peak_inds[i])
#                     min_ind = interval[0] + np.argmin(force[interval])
#                     # plt.plot(time[min_ind], force[min_ind], 'x', color='b')
#                 #Identify Ascend-indeces
#                 ascend_ind = range(min_ind, peak_inds[i], 1)
#                 # If ascend is long enough, so that percent_for_fit >= min_length:
#                 if len(ascend_ind) > min_length/percent_for_fit*100:
#                     fit_length = int(float(len(ascend_ind))/100*percent_for_fit)
#                     # If ascend is too short, but longer than min_length:
#                 elif len(ascend_ind) < min_length/percent_for_fit*100 and len(ascend_ind) > min_length:
#                     fit_length = min_length
#                 else:
#                     #If ascend is shorter than min_length
#                     fit_length = len(ascend_ind)
#                 fit_ind = range(peak_inds[i]-fit_length, peak_inds[i], 1)
#                 slope, intercept, r_value, p_value, std_err = linregress(time[fit_ind], force[fit_ind])
#                 cluster_by_boarder = assign_single_by_boarders(entry['x_peaks'][i])
                
#                 f.write(entry['name']+','
#                         +str(cluster_by_boarder)+','
#                         +str(entry['x_peaks'][i])+','
#                         +str(entry['y_peaks'][i])+','
#                         +str(slope)+','
#                         +sensitivity+','
#                         +spring_const
#                         +'\n')
            
                
#             #     line = slope*time[fit_ind]+intercept
#             #     plt.plot(time[fit_ind], force[fit_ind], '.', color='b')
#             #     plt.plot(time[fit_ind], line, color='r')
#             # plt.ylim(ymin = 0, ymax = 400)