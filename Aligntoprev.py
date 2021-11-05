#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 16:45:42 2020

@author: rnoah
"""
import numpy as np

from Align import align
from Binning import binning

def aligntoprev(FD_Curves, FD_Curves_prev, binsize_fine, Space, trafo_tresh):
    TotalNrOfCurves = len(FD_Curves)

    # Binning of each curve
    for index in range(TotalNrOfCurves):
        if Space == 'Force-CL':
            filler = trafo_tresh
            FD_Curves[index]['fine_bins'] = binning(FD_Curves[index]['x_Trafo'], FD_Curves[index]['y_Trafo'], binsize_fine, filler)
        else:
            filler = 0
            FD_Curves[index]['fine_bins'] = binning(FD_Curves[index]['x_PF'], FD_Curves[index]['y_PF'], binsize_fine, filler)
        FD_Curves[index]['Glob_Align_Score'] = 0


    # Find curve with lowest Global Alignment score in prev_set
    Glob_Align_Min = np.inf
    Best_aligned_Curve = 0
    for i in range(len(FD_Curves_prev)):
        if FD_Curves_prev[i]['Glob_Align_Score'] < Glob_Align_Min:
            Best_aligned_Curve = i
            Glob_Align_Min = FD_Curves_prev[i]['Glob_Align_Score']


    # Align all curves against optimal Global Aligned curve in prev_set
    print('Alignment to curve Nr', Best_aligned_Curve, FD_Curves_prev[Best_aligned_Curve]['name'])
    for i in range(TotalNrOfCurves):
        print(i+1, '/', TotalNrOfCurves)
        seq_i = FD_Curves[i]['fine_bins']
        seq_j = FD_Curves_prev[Best_aligned_Curve]['fine_bins']
        shift, score = align(seq_i, seq_j, binsize_fine, Space, trafo_tresh)
        FD_Curves[i]['Auto_Shift'] = shift*binsize_fine
    return FD_Curves