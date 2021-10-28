#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 17:19:28 2020

@author: rnoah
"""
import matplotlib.pyplot as plt

def keypress(event):
    global key
    key = event.key

def draw_shift(FD_Curves, global_shift, iteration, ax, nr_of_selected, Space):
    """
    Called by shift_preview function
    removes all curves from plot and redraws them.
    """

    if iteration != 0:
        ax.clear()
        if Space == 'Force-CL':
            ax.set_xlabel('Contour Length [nm]')
        else:
            ax.set_xlabel('Distance [nm]')
            ax.set_ylim(ymin = -50, ymax=300)
        ax.set_ylabel('Force [pN]')
        # plt.axvline(x=146.16, ymin=0, ymax=1) # Can be used to draw ref v_line
    for i in range(len(FD_Curves)):
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
            ax.plot(curve_x, curve_y, '.', markersize = 0.5, color='r', alpha = 3/(nr_of_selected+1))

def shift_preview(FD_Curves, Space):
    """
    Plots a preview of all shifted curves
        left: -1, z: -10 (lowers global shift)
        right: +1, x: +10 (increases global shift)
        up: Use current settings and proceed to clustering
    draw-shift function is called to update plot
    """
    global key
    TotalNrOfCurves = len(FD_Curves)
    global_shift = 0
    fig, ax = plt.subplots()
    fig.suptitle('Global Shift = '+str(global_shift) + ' nm')
    if Space == 'Force-CL':
        ax.set_xlabel('Contour Length [nm]')
    else:
        ax.set_xlabel('Distance [nm]')
        ax.set_ylim(ymin = -50, ymax=300)
    ax.set_ylabel('Force [pN]')
    # plt.axvline(x=146.16, ymin=0, ymax=1) # Can be used to draw ref v_line
    stay_in_loop = True
    iteration = 0
    nr_of_selected = 0
    for i in range(TotalNrOfCurves):
        if FD_Curves[i]['Selected'] == True:
            nr_of_selected += 1
    draw_shift(FD_Curves, global_shift, iteration, ax, nr_of_selected, Space)
    while stay_in_loop == True:
        key = '' # Reinitialize key in every loop
        fig.canvas.mpl_connect('key_press_event', keypress)
        while key not in ['left', 'right', 'up', 'z', 'x']:
            while not plt.waitforbuttonpress(0): pass
        if key == 'left':
            global_shift -= 1
        elif key == 'z':
            global_shift -= 10
        elif key == 'right':
            global_shift += 1
        elif key == 'x':
            global_shift += 10
        elif key == 'up':
            stay_in_loop = False
        fig.suptitle('Global Shift = '+str(global_shift) + ' nm')
        iteration+=1
        draw_shift(FD_Curves, global_shift, iteration, ax, nr_of_selected, Space)
        fig.canvas.draw()
    return global_shift