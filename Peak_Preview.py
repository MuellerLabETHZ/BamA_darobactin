#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 11:51:55 2020

@author: rnoah
"""
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


def keypress(event):
    global key
    key = event.key

def draw_peaks(FD_Curves, prom, iteration, nr_of_example_curves, axs):
    """
    Called by peak_preview function
    in first iteration:
        plots raw and PF-data, then detects peaks with current settings and plots them
    in all other iterations:
        removes previous peaks from plot and plots newly detected peaks instead 
    """
    for i in range(nr_of_example_curves):
        curve_x = FD_Curves[i]['x_raw']
        curve_y = FD_Curves[i]['y_raw']
        if iteration != 0:
            axs[i].lines.remove(axs[i].lines[-1])
        if iteration == 0:
            axs[i].set_ylim([-10, max(FD_Curves[i]['y_raw'])+30])
            axs[i].plot(curve_x, curve_y, color='k')
            axs[i].plot(FD_Curves[i]['x_PF'], FD_Curves[i]['y_PF'], color='b')
        peaks_PF, _ = find_peaks(FD_Curves[i]['y_PF'], prominence = prom)
        x_peaks_PF = [curve_x[a] for a in peaks_PF]
        y_peaks_PF = [curve_y[a] for a in peaks_PF]
        axs[i].plot(x_peaks_PF, y_peaks_PF, 'x', color='r')

def peak_preview(FD_Curves, nr_of_example_curves = 5, prom = 40):
    """
    Plots a preview of some curves to adjust peakfinder-parameters
    Sensitivity for peaks can be adjusted by changing required prominence:
        left: -1, z: -10 (lowers necessary prominence -> detects more peaks)
        right: +1, x: +10 (increases necessary prominence -> detects less peaks)
        up: Use current settings and proceed to curve-alignment
    draw-peaks function is called to update plot
    """
    global key
    print('Adjust the peak-detection sensitivity:')
    print('Left arrow key: a bit more sensitive')
    print('Right arrow key: a bit less sensitive')
    print('z: much more sensitive')
    print('x: much less sensitive')
    fig, axs = plt.subplots(nr_of_example_curves, sharex=True, gridspec_kw={'hspace': 0})
    fig.suptitle('Peak-Prominence = '+str(prom)+' -lower for more peaks, increase for less peaks')
    axs[-1].set_xlabel('Distance [nm]')
    axs[-1].set_ylabel('Force [pN]')
    stay_in_loop = True
    iteration = 0
    draw_peaks(FD_Curves, prom, iteration, nr_of_example_curves, axs)
    while stay_in_loop == True:
        key = '' # Reinitialize key in every loop
        fig.canvas.mpl_connect('key_press_event', keypress)
        print(key, prom)
        while key not in ['left', 'right', 'up', 'z', 'x']:
            while not plt.waitforbuttonpress(0): pass
        if key == 'left':
            prom -= 1
        elif key == 'z':
            prom -= 10
        elif key == 'right':
            prom += 1
        elif key == 'x':
            prom += 10
        elif key == 'up':
            stay_in_loop = False
        fig.suptitle('Peak-Prominence = '+str(prom)+' -lower for more peaks, increase for less peaks')
        iteration+=1
        draw_peaks(FD_Curves, prom, iteration, nr_of_example_curves, axs)
        fig.canvas.draw()
    return prom