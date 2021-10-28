#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy import stats
import copy

def BLC(x_input, y_input, PercentForBLC = 10, name = False):
    """Estimates tilt in last x percent (3rd arg, or 10%),
    fits a line and corrects all y-vals to the line.
    Arg0 = x_val
    Arg1 = y_val
    Arg2 = Last x percent for fit (optional) - default = 10%
    Arg3 = Name (optional) - only required to plot
    """
    
    #Make independent copies of input
    x_val = copy.deepcopy(x_input)
    y_val = copy.deepcopy(y_input)

    ### Defines Start-Position for linear-regression
    StartInArray = int(len(x_val)* float((100 - PercentForBLC)/100))
    
    ### Linear regression of last X percent
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_val[StartInArray:],y_val[StartInArray:])
    line = slope*x_val+intercept
    
    ### Baseline-corrected y_vals
    y_BLC = y_val - line
    
    if name != False:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.title(name + ' - BLC')
        plt.xlabel('Tip-Sample Separation [nm]')
        plt.ylabel('Vertical Deflection [pN]')
        plt.plot(x_val,y_val)
        plt.plot(x_val, line)
        plt.plot(x_val[:StartInArray], y_BLC[:StartInArray])
        plt.plot(x_val[StartInArray:],y_BLC[StartInArray:])
        plt.axis([min(x_val)*1.1,max(x_val)*1.1,min(y_val)*1.1,max(y_val)*1.1])
        plt.show()
    
    return x_val, y_BLC