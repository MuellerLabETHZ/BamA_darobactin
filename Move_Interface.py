#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 11:55:22 2020

@author: rnoah
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from MoveCurve_200520 import plot_curve
from MovePeaks_200908 import plot_peaks
from Transformation import Trafo




def keypress(event):
    global key
    key = event.key

def move_interface(FD_Curves, FD_Curves_prev, Space, Allow_y_shift, prom, roi):
    global key
    if Space == 'Force-CL':
        fig, [ax, ax2] = plt.subplots(2,1, gridspec_kw={'hspace': 0.4})
        ax.set_xlabel('Contour Length [nm]')
        ax.set_ylabel('Force [pN]')
        ax2.set_xlabel('Distance [nm]')
        ax2.set_ylabel('Force [pN]')
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Distance [nm]')
        ax.set_ylabel('Force [pN]')
    TotalNrOfCurves = len(FD_Curves)
    event_hist = []
    BG_Curves = []
    i = 0
    nr_of_selected = 0

    #Plot selected curves of prev_set to BG
    if FD_Curves_prev != False:
        for index in range(len(FD_Curves_prev)):
            if FD_Curves_prev[index]['Selected'] == True:
                auto_shift = FD_Curves_prev[index]['Auto_Shift']
                x_shift = FD_Curves_prev[index]['x_shift']
                y_shift = FD_Curves_prev[index]['y_shift']
                if Space == 'Force-CL':
                    curve_x = [a+auto_shift+x_shift for a in FD_Curves_prev[index]['x_Trafo']]
                    curve_y = [a+y_shift for a in FD_Curves_prev[index]['y_Trafo']]
                else:
                    curve_x = [a+auto_shift+x_shift for a in FD_Curves_prev[index]['x_raw']]
                    curve_y = [a+y_shift for a in FD_Curves_prev[index]['y_raw']]
                prev_curve, = ax.plot(curve_x, curve_y, '.', markersize = 0.5, color='red')
                BG_Curves.append(prev_curve)
                nr_of_selected += 1
        for element in BG_Curves:
            element.set_alpha(3/(nr_of_selected+1))


    while i < TotalNrOfCurves+1:
        print(i, '/',  TotalNrOfCurves, ' Aligned')
        key = '' # Reinitialize key in every loop
        if i < TotalNrOfCurves:
            #Define data to align, based on selected space
            auto_shift = FD_Curves[i]['Auto_Shift']
            if Space == 'Force-CL':
                curve_x = [a+auto_shift for a in FD_Curves[i]['x_Trafo']]
                curve_y = FD_Curves[i]['y_Trafo']
            else:
                curve_x = [a+auto_shift for a in FD_Curves[i]['x_raw']]
                curve_y = FD_Curves[i]['y_raw']
            # Find starting x and y values, to later define shift
            x0 = curve_x[0]
            y0 = curve_y[0]
            # Plot title for each Curve
            fig. suptitle(str(i+1)+'/'+ str(TotalNrOfCurves)+' '+FD_Curves[i]['name'])
            # Runs Program to manually drag curves
            dr, curve = plot_curve(ax, curve_x, curve_y) # Runs DraggableCurve
            # Find Peak indeces and x and y positions in both spaces
            peaks, _ = find_peaks(FD_Curves[i]['y_PF'], prominence = prom)
            if Space == 'Force-CL':
                ref_curve, = ax2.plot(FD_Curves[i]['x_raw'], FD_Curves[i]['y_raw'], 'k')
                ax2.set_ylim(ymin = -50, ymax = roi)
            if len(peaks) != 0:
                x0_peaks = [FD_Curves[i]['x_raw'][a] for a in peaks]
                y0_peaks = [FD_Curves[i]['y_raw'][a] for a in peaks]
            else:
                # If no peaks are detected, a peak is added to prevent errors
                peaks = np.array([100])
                x0_peaks = np.array([100])
                y0_peaks = np.array([100])
                print(peaks.shape, x0_peaks.shape, y0_peaks.shape)
            if Space == 'Force-CL':
                x0_peaks, y0_peaks = Trafo(np.asarray(x0_peaks), np.asarray(y0_peaks), 0)
            else:
                ax.set_ylim(ymin = -50)
            ax.set_xlim(xmin = 0, xmax = roi)
            ax.set_ylim(ymax = roi)
            # Plots peaks onto graph and add autoshift to x_peaks
            x0_peaks = [a+auto_shift for a in x0_peaks]
            dr_peaks, peakdraw = plot_peaks(ax, x0_peaks, y0_peaks, curve)
            # Update Plot and wait for user input
            fig.canvas.draw() # Updates Plot
            fig.canvas.mpl_connect('key_press_event', keypress)
            #remove current peaks from graph to plot new ones later
            while key not in ['left', 'right', 'down']:
                while not plt.waitforbuttonpress(0): pass
            print(key)
            ax.lines.remove(peakdraw)
            if Space == 'Force-CL':
                ax2.lines.remove(ref_curve)
            # Get current shift-value and shift curves and peaks by dx and dy
            x1 = curve.get_xdata(orig=True)[0]
            y1 = curve.get_ydata(orig=True)[0]
            dx = x1-x0
            dy = y1-y0
            curve_x_shifted = [a+dx for a in curve_x]
            curve_y_shifted = [a+dy for a in curve_y]

            # Now get all peaks from graph:
            x1_peaks = np.asarray(dr_peaks[0].peaks.get_xdata())
            y1_peaks = np.asarray(dr_peaks[0].peaks.get_ydata())
            if Allow_y_shift == False:
                y1_peaks = y1_peaks - dy
            # Sort peaks based on their x-positions
            sorted_peak_ind = np.argsort(x1_peaks)
            x1_peaks = x1_peaks[sorted_peak_ind]
            y1_peaks = y1_peaks[sorted_peak_ind]

            # Write data to FD-Curves dictionary
            FD_Curves[i]['x_shift'] = dx
            FD_Curves[i]['y_shift'] = dy
            FD_Curves[i]['i_peaks'] = peaks
            FD_Curves[i]['x_peaks'] = x1_peaks
            FD_Curves[i]['y_peaks'] = y1_peaks
        # Add new Curve and display current in BG
            if key == 'right':
                nr_of_selected += 1
                ax.lines.remove(curve)
                if Allow_y_shift == True:
                    BG_Curve, = ax.plot(curve_x_shifted, curve_y_shifted, '.', markersize = 0.5, color='red')
                else:
                    BG_Curve, = ax.plot(curve_x_shifted, curve_y, '.', markersize = 0.5, color='red')
                BG_Curves.append(BG_Curve)
                FD_Curves[i]['Selected'] = True
                for element in BG_Curves:
                    element.set_alpha(3/(nr_of_selected+1))
                event_hist.append('right')
                i += 1
        # If at the start nothing changes (curve from dotheplot is removed)
            elif key == 'left' and i == 0:
                ax.lines.remove(curve)
        # Go Back to previous Curve and remove curve in front and BG
            elif key == 'left' and event_hist[-1]!='down':
                nr_of_selected -= 1
                ax.lines.remove(curve)
                ax.lines.remove(BG_Curves[-1])
                BG_Curves = BG_Curves[:-1]
                del event_hist[-1]
                i -= 1
            # Go Back to previous Curve and remove only curve in front
            elif key == 'left' and event_hist[-1]=='down':
                ax.lines.remove(curve)
                del event_hist[-1]
                i -= 1
        # Skip this curve -> Curve is ignored/deleted
            elif key == 'down':
                FD_Curves[i]['Selected'] = False
                ax.lines.remove(curve)
                event_hist.append('down')
                i += 1
        else:
            fig.suptitle('finished! Happy? Yes -> press UP, otherwhise -> press left')
            fig.canvas.draw() # Updates Plot
            fig.canvas.mpl_connect('key_press_event', keypress)
            while key not in ['up', 'left']:
                while not plt.waitforbuttonpress(0): pass
            if key == 'left':
                ax.lines.remove(BG_Curves[-1])
                BG_Curves = BG_Curves[:-1]
                i -= 1
            elif key == 'up':
                i+= 1
    return FD_Curves