#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 15:45:56 2020

@author: rnoah
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os


def pathway(clustering, FD_Curves, Space, color_scheme, plot_dir, font_size, fig_size, DPI, Title):
    TotalNrOfCurves = len(FD_Curves)

    # Assigns each peak of each FD_Curve to as cluster
    for cluster in clustering:
        for i in range(TotalNrOfCurves):
            if 'clusters' not in FD_Curves[i]:
                FD_Curves[i]['clusters'] = {}
            FD_Curves[i]['clusters'][cluster] = []
            x_peaks = FD_Curves[i]['x_peaks']
            y_peaks = FD_Curves[i]['y_peaks']
            for index, entry in enumerate(x_peaks):
                if entry in clustering[cluster]['x']:
                    FD_Curves[i]['clusters'][cluster] = [entry, y_peaks[index]]


    # Clustering Matrix: cols=Clusters, rows=FD-Curves, 1 if curve goes through cluster, else zero
    m = TotalNrOfCurves
    n = len(clustering) # Subtract 1 because noise-points with cluster -1 are not included
    Clustering_Matrix = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            if len(FD_Curves[i]['clusters'][j]) != 0:
                Clustering_Matrix[i, j] = 1
    
    # Transition Matrix: rows and cols=clusters, each entry is nr of transition from one cluster to another
    Transition_Matrix = np.zeros((n, n))
    for i in range(m):
        transitions = []
        for j in range(n):
            if Clustering_Matrix[i, j] == 1:
                transitions.append(j)
        for index, entry in enumerate(transitions):
            if index != 0:
                Transition_Matrix[transitions[index-1], entry] += 1
    
    # Plot pathways: x,-axis: avg dist, y-axis: avg forces
    # markersize=nr of points in cluster, linewidth=nr of transitions
    min_transition = 10000
    max_transition = np.amax(Transition_Matrix)
    total_transitions = np.sum(Transition_Matrix)
    total_clusterpoints = np.sum(Clustering_Matrix)
    min_cluster = min([np.sum(Clustering_Matrix[:,j]) for j in range(n)])
    max_cluster = max([np.sum(Clustering_Matrix[:,j]) for j in range(n)])
    fig, ax = plt.subplots(figsize=fig_size)
    fig.suptitle(Title, size=20)
    if Space == 'Force-CL':
        ax.set_xlabel('Contour Length [nm]', size=font_size)
    else:
        ax.set_xlabel('Distance [nm]', size=font_size)
    ax.set_ylabel('Force [pN]', size=font_size)
    
    ax.set_xlim(xmin = 0, xmax = 300)
    ax.set_ylim(ymin = 0, ymax = 300)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    for i in range(n):
        for j in range(i+1, n):
            if Transition_Matrix[i, j] > total_transitions/100:
                if Transition_Matrix[i, j] < min_transition:
                    min_transition = Transition_Matrix[i, j]
                x_from = np.average(clustering[i]['x'])
                x_to = np.average(clustering[j]['x'])
                y_from = np.average(clustering[i]['y'])
                y_to = np.average(clustering[j]['y'])
                # ax.plot([x_from, x_to], [y_from, y_to], linewidth=5*Transition_Matrix[i, j]/max_transition, color='lightgrey', label='lines')
                ax.plot([x_from, x_to], [y_from, y_to], linewidth=(Transition_Matrix[i, j]/total_transitions)*100, color='lightgrey', label='lines')
    
    
    for j in range(n):
        x_avg = np.average(clustering[j]['x'])
        y_avg = np.average(clustering[j]['y'])
        # ax.plot(x_avg, y_avg, 'o', markersize=12*np.sum(Clustering_Matrix[:,j])/max_cluster, label='marks')
        ax.plot(x_avg, y_avg, 'o', markersize=(2*np.sum(Clustering_Matrix[:,j])/total_clusterpoints)*100, label='marks', color=color_scheme[j], markeredgecolor='k', linewidth=0.5)
    
    legend_elements = [Line2D([0], [0], color='lightgrey', lw=(min_transition/total_transitions)*100, label=int(min_transition)),
                       Line2D([0], [0], color='lightgrey', lw=(max_transition/total_transitions)*100, label=int(max_transition)),
                       Line2D([0], [0], marker='o', color='w', label=int(min_cluster), markerfacecolor='lightgrey', markersize=(2*min_cluster/total_clusterpoints)*100, markeredgecolor='k', markeredgewidth=.8),
                       Line2D([0], [0], marker='o', color='w', label=int(max_cluster), markerfacecolor='lightgrey', markersize=(2*max_cluster/total_clusterpoints)*100, markeredgecolor='k', markeredgewidth=.8)]
    ax.legend(handles=legend_elements, handletextpad=2, labelspacing=1, borderpad=1, fontsize=15, loc='upper left')
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.savefig(os.path.join(plot_dir, 'Pathway.png'), dpi=DPI)