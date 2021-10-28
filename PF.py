#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import copy


def PF(x_input, y_input, ppnm, Thr = 5):
    """Filters Peaks based on their width.
    Thin peaks and positive values are set to zero.
    Arg0 = x_val
    Arg1 = y_val
    Arg2 = ppnm (points per nm)
    Arg3 = Thr - is multiplied with ppnm,
    to define peakwidth-threshold (optional) Default = 5,
    everything below is not considered.
    """

    #Make independent copies of input
    x_vals = copy.deepcopy(x_input)
    y_vals = copy.deepcopy(y_input)
    
    Threshold = Thr * ppnm
                
    ### Sets all positive values to zero.
    for index, value in enumerate(y_vals):
        if value > 0:
            y_vals[index] = 0
    
    ### Variables used to filter peaks based on their width.
    InPeak = False # True if last value was part of a peak. Else it is False.
    count = 0 # Counts peaks of each rectraction curve.
    peaks = {} # Dictionary that contains All peaks (Unfiltered).
    PeakList = [] # Pooled indeces of broad peaks.
            
    ### Goes through every value in y_vals and detects if it is at the start, end or within a peak.
    ### Saves the corresponding indeces in lists (peakX).
    for index, value in enumerate(y_vals):
        
        ###Start of a Peak
        if value < 0 and InPeak == False:
            InPeak = True
            peaks['peak' + str(count)] = []
            peaks['peak' + str(count)].append(index)
        
        ###During Peak
        elif value < 0 and InPeak == True:
            peaks['peak' + str(count)].append(index)
            
        ###End of Peak
        elif value == 0 and InPeak == True:
            InPeak = False
            count += 1
    
    ### If a peak is broader than the Threshold it is added to PeakList
    for number in peaks:
        if len(peaks[number]) > Threshold:
            PeakList += peaks[number]

    ### Sets all values at indeces not contained in PeakList to zero. Only broad peaks stay.
    for index, value in enumerate(y_vals):
        if index not in PeakList:
            y_vals[index] = 0

    return x_vals, y_vals
