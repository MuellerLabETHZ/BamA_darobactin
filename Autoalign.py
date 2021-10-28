#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 16:45:42 2020

@author: rnoah
"""
import numpy as np

from Align import align
from Binning import binning

def autoalign(FD_Curves, binsize_coarse, binsize_fine, Space, trafo_tresh):
    TotalNrOfCurves = len(FD_Curves)

    # Binning of each curve
    for index in range(TotalNrOfCurves):
        if Space == 'Force-CL':
            filler = trafo_tresh
            FD_Curves[index]['coarse_bins'] = binning(FD_Curves[index]['x_Trafo'], FD_Curves[index]['y_Trafo'], binsize_coarse, filler)
            FD_Curves[index]['fine_bins'] = binning(FD_Curves[index]['x_Trafo'], FD_Curves[index]['y_Trafo'], binsize_fine, filler)
        else:
            filler = 0
            FD_Curves[index]['coarse_bins'] = binning(FD_Curves[index]['x_PF'], FD_Curves[index]['y_PF'], binsize_coarse, filler)
            FD_Curves[index]['fine_bins'] = binning(FD_Curves[index]['x_PF'], FD_Curves[index]['y_PF'], binsize_fine, filler)
        FD_Curves[index]['Glob_Align_Score'] = 0

    # Align all curves against each other
    NrofAlignments = sum([a for a in range(TotalNrOfCurves)])
    counter = 0
    print('Cross-Alignment:')
    for i in range(TotalNrOfCurves):
        for j in range(i):
            counter += 1
            print(counter, '/', NrofAlignments)
            seq_i = FD_Curves[i]['coarse_bins']
            seq_j = FD_Curves[j]['coarse_bins']
            shift, score = align(seq_i, seq_j, binsize_coarse, Space, trafo_tresh)
            FD_Curves[i]['Glob_Align_Score'] += score
            FD_Curves[j]['Glob_Align_Score'] += score

    # Find curve that fits best to all others using coarse bins
    Glob_Align_Min = np.inf
    Best_aligned_Curve = 0
    for i in range(TotalNrOfCurves):
        if FD_Curves[i]['Glob_Align_Score'] < Glob_Align_Min:
            Best_aligned_Curve = i
            Glob_Align_Min = FD_Curves[i]['Glob_Align_Score']
    FD_Curves[Best_aligned_Curve]['Auto_Shift'] = 0

    # Align all curves against optimal aligned curve
    print('Alignment to curve Nr', Best_aligned_Curve, FD_Curves[Best_aligned_Curve]['name'])
    for i in range(TotalNrOfCurves):
        print(i+1, '/', TotalNrOfCurves)
        if i != Best_aligned_Curve:
            seq_i = FD_Curves[i]['fine_bins']
            seq_j = FD_Curves[Best_aligned_Curve]['fine_bins']
            shift, score = align(seq_i, seq_j, binsize_fine, Space, trafo_tresh)
            FD_Curves[i]['Auto_Shift'] = shift*binsize_fine
    return FD_Curves