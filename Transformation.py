#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 11:10:55 2018

@author: rnoah
"""

def Trafo(x_vals, y_vals, y_threshold):
    P = 0.4 # Persistence Length = 0.4 nm
    T=4.1
    CL_list = []
    indeces = []
    for index, val in enumerate(x_vals):
        if y_vals[index] > y_threshold:
            z = x_vals[index]
            F = y_vals[index]
            k = F*P/T
        
            # Fractions can be calculated and filled into formula to increase speed
            A = (16*k**3-72*k**2+27*k-54)
            B = (-72*k**3-162*k+54*k**2)
            C = 1/186624*A**2/k**6+1/272097792*B**3/k**9
            D = -1/432*A/k**3
            E = complex(C)**(1/2)
            
            CL = (((D+E)**(1/3)+(D-E)**(1/3)+1/2/k+2/3)*z).real
    #            if CL < 450:
            CL_list.append(CL)
            indeces.append(index)
    return CL_list, y_vals[indeces]
