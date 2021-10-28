#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:57:18 2020

@author: rnoah
"""
import numpy as np
# Use this script to replace optics-clustering to assign clusters to peaks,
# Then proceed as normal

def assign_single_by_boarders(x_peak_pos):

    boarders = [-139.98015495938643,
                -114.02659506094548,
                -88.41534714029555,
                -63.16989894999022,
                -38.451182651128406,
                -16.446102731192042,
                -3.44891102454197,
                5.090348498573519,
                16.261381106364126,
                30.116869399237178,
                48.50258269275616,
                66.49272202005179,
                79.56143921888864,
                98.09892746370926,
                116.45820240167137,
                133.97092073043893]
    
    # Add 406aa to shift (146.16nm) to move first HP from zero to correct dist
    for i in range(len(boarders)):
        boarders[i] += 146.16

    for i in range(len(boarders)-1):
        if x_peak_pos > boarders[i] and x_peak_pos < boarders[i+1]:
            return i
    return -1


def splitbyboarders(all_peaks_x, all_peaks_y, noise_percents):

    boarders = [-139.98015495938643,
                -114.02659506094548,
                -88.41534714029555,
                -63.16989894999022,
                -38.451182651128406,
                -16.446102731192042,
                -3.44891102454197,
                5.090348498573519,
                16.261381106364126,
                30.116869399237178,
                48.50258269275616,
                66.49272202005179,
                79.56143921888864,
                98.09892746370926,
                116.45820240167137,
                133.97092073043893]

    # Add 406aa to shift (146.16nm) to move first HP from zero to correct dist
    for i in range(len(boarders)):
        boarders[i] += 146.16

    clusters = {}
    is_in_cluster = []
    for i in range(len(boarders)-1):
        clusters[i] = {}
        clusters[i]['x'] = []
        clusters[i]['y'] = []
        for index, x_pos in enumerate(all_peaks_x):
            # if x_pos > boarders[i] and x_pos < boarders[i+1]:
            # percentage = np.linspace(boarders[i], boarders[i+1], 100)
            # if x_pos > percentage[noise_percents] and x_pos > percentage[-(noise_percents+1)]:
            interval_length = boarders[i+1] - boarders[i]
            noise_length = interval_length/100*noise_percents
            lower_boarder = boarders[i] + noise_length
            upper_boarder = boarders[i+1] - noise_length
            if x_pos > lower_boarder and x_pos < upper_boarder:
                clusters[i]['x'].append(all_peaks_x[index])
                clusters[i]['y'].append(all_peaks_y[index])
                is_in_cluster.append(index)
    
    noise_points = {'x':[], 'y':[]}
    for index, x_pos in enumerate(all_peaks_x):
        if not index in is_in_cluster:
            noise_points['x'].append(all_peaks_x[index])
            noise_points['y'].append(all_peaks_y[index])
    return clusters, noise_points, boarders