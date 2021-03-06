#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 18:49:08 2021

@author: rnoah
"""
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
import os


# def create_col_map(color_resolution):
#     # Color map red (low), to white, to blue (high)
#     R = 1
#     G = 0
#     B = 0
#     color_map = []
#     for i in range(color_resolution):
#         B = i/color_resolution
#         G = i/color_resolution
#         color_map.append([R, G, B])
#     B = 1
#     for i in range(color_resolution, 0, -1):
#         R = i/color_resolution
#         G = i/color_resolution
#         color_map.append([R, G, B])
#     B = 0
#     G = 0
#     color_map.append([R, G, B])
#     return color_map

def create_col_map(color_resolution):
    # Color map white (low) to red (high)
    R = 1
    G = 1
    B = 1
    color_map = []
    for i in range(color_resolution, 0, -1):
        B = i/color_resolution
        G = i/color_resolution
        color_map.append([R, G, B])
    return color_map


def interpolation(segments, params, start, protein_length, resolution):
    x = np.asarray(segments)
    y = params.to_numpy()
    
    #Extend because of edge-problem
    x = np.insert(x, 0, start)
    x = np.append(x, protein_length)
    y = np.insert(y, 0, y[0])
    y = np.append(y, y[-1])
    
    y_range = max(y) - min(y)
    
    f = interpolate.interp1d(x, y, kind='quadratic')
    
    xnew = np.arange(min(x), max(x)+1, 1)
    ynew = np.arange(min(y), max(y), y_range/((resolution*2)+1))
    ynew = f(xnew)   # use interpolation function returned by `interp1d`

    return xnew, ynew


def map_to_rgb(params_cont, color_map, map_range, save_to_dir, poi):
    fig, ax = plt.subplots()
    ax.set_ylim(map_range[0], map_range[1])
    ax.set_xlabel('Residue')
    ax.set_ylabel('Parameter Value')
    mapped_colors = []
    # mapping = interpolate.interp1d([min(params_cont), max(params_cont)], [0, len(color_map)-1])
    mapping = interpolate.interp1d([map_range[0], map_range[1]], [0, len(color_map)-1])
    for index, val in enumerate(params_cont):
        mapped = mapping(val)
        mapped = int(np.round(mapped))
        rgb_col = color_map[mapped]
        mapped_colors.append(rgb_col)
        plt.plot(residues[index], params_cont[index], '.', color=rgb_col)
    plt.plot(segments_aa, params, 'o', color='k')
    plt.savefig(os.path.join('/Users/rnoah', save_to_dir, 'Plots', poi.replace(' ', '')))
    return mapped_colors

def lists_for_pymol(residues, mapped_colors):
    set_list = []
    coloration_list = []
    for i in range(len(residues)):
        set_list.append('set_color r' + str(residues[i]) + ', ' + str(mapped_colors[i]))
        coloration_list.append('color r' + str(residues[i]) + ', ' + ' Resi ' + str(residues[i]))    
    return set_list, coloration_list

def add_render_settings(coloration_list, save_to_dir, poi):
    #Optional if view should also be standardized
    # For unbound Structure
    # view = 'set_view (\
    #    -0.304142803,   -0.438496441,   -0.845667481,\
    #    -0.901852608,    0.418401033,    0.107412912,\
    #     0.306732208,    0.795361519,   -0.522753775,\
    #     0.002427377,   -0.001643389, -420.548797607,\
    #   -11.691201210,   47.120197296,    4.419474125,\
    #  -15818.869140625, 16659.865234375,  -20.000000000 )'
    
    # Use this for unfolding view
#     'set_view (\
#        -0.287582964,   -0.457459092,   -0.841405272,\
#     -0.907554150,    0.410770535,    0.086875401,\
#      0.305887133,    0.788629651,   -0.533340693,\
#      0.002130874,   -0.004011184, -784.260864258,\
#    -28.578849792,   63.210720062,   34.982677460,\
#   -14050.592773438, 15619.072265625,  -20.000000000 )'

    # For Daro-Bound Structure
    view = 'set_view (\
    -0.098090716,   -0.529967666,    0.842271328,\
    -0.571739435,   -0.662703097,   -0.483583689,\
     0.814499080,   -0.529034376,   -0.238000110,\
     0.002358153,    0.023387626, -497.144622803,\
   172.756408691,  178.151794434,  176.638854980,\
  -1857.567260742, 2851.541259766,  -20.000000000 )'

    edge_col = 'ray_trace_color, black'
    individual_path = os.path.join(save_to_dir, poi.replace(' ', ''))
    render_to_file = 'png ~/' + individual_path + '.png, width=4cm, height=7cm, dpi=300, ray=1'
    
    coloration_list.append(view)
    coloration_list.append(edge_col)
    coloration_list.append(render_to_file)
    return coloration_list

    
#%%
params_dir = '/Users/rnoah/Desktop/DFS_Aligned/Prism_Analysis/params_and_errors.csv'
save_to_dir = 'Documents/Manuscripts/BamA_DFS_MS/Pymol_Stuff/HeatMaps_Autogenerated'

resolution = 200
start = 21
protein_length = 810 # Length in aa
df_of_all_params = pd.read_csv(params_dir)
poi = 'D+ kappa'
params = df_of_all_params[poi]

# Only for k0:
# params = np.log10(params)


# Enter CL in aa of rupture events
segments_aa = [80, 148, 219, 288, 355, 404, 430, 459, 491, 535, 593, 633, 674, 728, 774]
segments_aa = [i + 20 for i in segments_aa] # Compensate for 20 aa introduced by signal seq in pymol

color_map = create_col_map(resolution)

residues, params_cont = interpolation(segments_aa, params, start, protein_length, resolution)

map_range = [0, 6] # Manual Range
# map_range = [min(params_cont), max(params_cont)] # Auto Range
# print(map_range)

mapped_colors = map_to_rgb(params_cont, color_map, map_range, save_to_dir, poi)

a_set_list, b_coloration_list = lists_for_pymol(residues, mapped_colors)
b_coloration_list = add_render_settings(b_coloration_list, save_to_dir, poi)