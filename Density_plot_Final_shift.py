#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 16:12:51 2020

@author: rnoah
"""
import matplotlib.pyplot as plt
import numpy as np
import os

def density_plot(FD_Curves, roi, bins_per_nm, Space, Title, color_scheme, name_scheme, clusterdata, plot_dir, font_size, DPI):
    x_shift = 0
    y_shift = 0
    bins_per_nm = 0.2
    resolution = int(300*bins_per_nm)
    
    fig, ax = plt.subplots(figsize=[7.5, 6.1])
    fig.suptitle(Title, size=20, y=1.0)
    if Space == 'Force-CL':
        ax.set_xlabel('Contour Length [nm]', size=font_size)
    else:
        ax.set_xlabel('Distance [nm]', size=font_size)
        ax.set_xlim(xmin = 0, xmax = 300)
        ax.set_ylim(ymin = 0, ymax = 350)
    ax.set_ylabel('Force [pN]', size=font_size)
    x = []
    y = []
    for i in range(len(FD_Curves)):
        if FD_Curves[i]['Selected'] == True:
            x_temp = FD_Curves[i]['curve_x_shifted']
            y_temp = FD_Curves[i]['curve_y_shifted']
            x_trunc = [i for i in x_temp if i < roi]
            y_trunc = y_temp[:len(x_trunc)]
            x += x_trunc
            y += y_trunc

    x = np.asarray(x)
    y = np.asarray(y)
    
    # Make the plot
    im = plt.hexbin(x, y, gridsize=(resolution*7, resolution*6), cmap=plt.cm.Greys)
    max_bin = int(max(im.get_array()))
    plt.text(310, 200, '0', color='k', size=18)
    plt.text(310, 350, str(max_bin), color='k', size=18)
    # plt.colorbar()
    plt.xlim(xmin = 0, xmax=300)
    plt.ylim(ymin = 20, ymax=350)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    # plt.title(Title, pad=30)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    n = 0
    for entry in FD_Curves:
        if entry['Selected'] == True:
            n += 1
    plt.text(0.95, 0.99, 'n='+str(n), fontsize=16, horizontalalignment='right',
              verticalalignment='top', transform=ax.transAxes)
    for i in range(len(clusterdata)):
        x_pos = clusterdata[i][1]
        x_std = clusterdata[i][2]
        y_pos = clusterdata[i][3] + clusterdata[i][4]
        plt.plot(x_pos, y_pos, 'v', color=color_scheme[i], markersize=9, markeredgewidth=0.35, markeredgecolor='k')
        plt.text(x_pos-3.5, y_pos+11, str(int(np.round(x_pos/0.36)))
                 + ' ± ' + str(int(np.round(x_std)))
                 , color='k', rotation='vertical', size=13)

    plt.savefig(os.path.join(plot_dir, 'Density_Plot.png'), dpi=DPI)