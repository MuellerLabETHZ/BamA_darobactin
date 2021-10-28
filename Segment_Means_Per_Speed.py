#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 11:04:24 2021

@author: rnoah
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import random

input_file = '/Users/rnoah/Desktop/DFS_Aligned/Prism_Analysis/Sorted_By_Segment.csv'
plot_dir = '/Users/rnoah/Desktop/DFS_Aligned/Segment_Means_Per_Speed'
df = pd.read_csv(input_file)
LR_cols = [col for col in df.columns if 'LR' in col]
df = df[df.columns.drop(LR_cols)]


Segments = ['S' + str(i) for i in range(15)]
Daro = ['-', '+']
Speeds = [500, 700, 1000, 3000, 4500, 6000]


markerdict = {
    500:'o',
    700:'v',
    1000:'s',
    3000:'D',
    4500:'^',
    6000:'h'
    }


bg_col1 = [0.9,0.9,0.9]
bg_col2 = [0.027, 0.94, 0.33]
mean_col1 = [0.68,0.68,0.68]
mean_col2 = [0.015, 0.63, 0.22]
x1 = 13.2
x2 = 33
bg_msize = 3
mean_msize = 12

# order = [0,1.5,3,4.5,6,7.5]
order = [0,2,4,6,8,10]
for seg in Segments:
    random.shuffle(order)
    fig, ax = plt.subplots(figsize=(2.45, 3.7))
    ax.set_xlim(xmin = 0, xmax=46)
    ax.set_ylim(ymin = 0, ymax=500)
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticklabels([])
    for cond in Daro:
        for index, speed in enumerate(Speeds):
            col_name = 'F_D' + cond + str(speed) + '_' + seg
            m = markerdict[speed]
            x_pos = x1 if cond == '-' else x2
            bg_col = bg_col1 if cond == '-' else bg_col2
            mean_col = mean_col1 if cond == '-' else mean_col2
            col_sum = 0
            col_count = 0
            for entry in df[col_name]:
                if not np.isnan(entry):
                    col_sum += entry
                    col_count += 1
                    plt.plot(x_pos+(np.random.normal(0,2)), entry, marker=m, color=bg_col, markersize=bg_msize, markeredgecolor='k', markeredgewidth=0.25, alpha=0.4, zorder=0)
            col_mean = col_sum/col_count
            plt.plot(x_pos-5+order[index], col_mean, marker=m, color=mean_col, markersize=mean_msize, markeredgecolor='k', markeredgewidth=1.5, zorder=100)
            #Don't always add random, but make symmetrical to both sides of x_pos
                
    plt.rcParams.update({
        "figure.facecolor":  (0.0, 0.0, 0.0, 0.0),
        "axes.facecolor":    (0.0, 0.0, 0.0, 0.0),
        "savefig.facecolor": (0.0, 0.0, 0.0, 0.0)
    })
    
    plt.savefig(os.path.join(plot_dir, 'L_shift'+seg+'_means.png'), dpi=600)
            

# input_file = '/Users/rnoah/Desktop/DFS_Aligned/Segment_Means_Per_Speed.csv'
# plot_dir = '/Users/rnoah/Desktop/DFS_Aligned/Segment_Means_Per_Speed'

# all_data = pd.read_csv(input_file)
# col1 = [0.68,0.68,0.68]
# col2 = [0.015, 0.63, 0.22]
# mew = 0.5
# x1 = 13.2
# x2 = 33
# msize=8

# for seg in all_data['Segment']:
#     fig, ax = plt.subplots(figsize=(0.48*5.1,0.725*5.1))
    
#     ax.set_xlim(xmin = 0, xmax=46)
#     ax.set_ylim(ymin = 0, ymax=500)
#     ax.spines['top'].set_visible(False)
#     ax.get_xaxis().set_ticks([])
#     ax.get_yaxis().set_ticklabels([])
    
#     plt.plot(x1, all_data['D-_1ums'][seg], 'o', color=col1, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x1, all_data['D-_3ums'][seg], 'v', color=col1, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x1, all_data['D-_4_5ums'][seg], 's', color=col1, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x1, all_data['D-_05ums'][seg], 'D', color=col1, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x1, all_data['D-_6ums'][seg], '^', color=col1, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x1, all_data['D-_07ums'][seg], 'h', color=col1, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
    
#     plt.plot(x2, all_data['D+_1ums'][seg], 'o', color=col2, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x2, all_data['D+_3ums'][seg], 'v', color=col2, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x2, all_data['D+_4_5ums'][seg], 's', color=col2, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x2, all_data['D+_05ums'][seg], 'D', color=col2, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x2, all_data['D+_6ums'][seg], '^', color=col2, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
#     plt.plot(x2, all_data['D+_07ums'][seg], 'h', color=col2, markersize=msize, markeredgecolor='k', markeredgewidth=mew)
    
#     plt.rcParams.update({
#         "figure.facecolor":  (0.0, 0.0, 0.0, 0.0),  # red   with alpha = 30%
#         "axes.facecolor":    (0.0, 0.0, 0.0, 0.0),  # green with alpha = 50%
#         "savefig.facecolor": (0.0, 0.0, 0.0, 0.0),  # blue  with alpha = 20%
#     })
    
#     plt.savefig(os.path.join(plot_dir, 'S'+str(seg)+'_means.png'), facecolor='None', dpi=600)