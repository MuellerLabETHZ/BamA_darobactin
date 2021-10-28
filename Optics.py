#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 15:44:59 2020

@author: rnoah
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import OPTICS
import matplotlib.gridspec as gridspec

def optics(all_peaks_x,
           all_peaks_y,
           y_scaling_factor,
           min_samples_input,
           xi_input,
           min_cluster_size_input,
           Space):

    scaled_peaks_y = [a/y_scaling_factor for a in all_peaks_y]
    X = np.transpose(np.asarray([all_peaks_x, scaled_peaks_y]))

    clust = OPTICS(min_samples=min_samples_input,
                    xi=xi_input,
                    min_cluster_size=min_cluster_size_input)

    clust.fit(X)
    clusters = clust.fit_predict(X)

    space = np.arange(len(X))
    reachability = clust.reachability_[clust.ordering_]
    labels = clust.labels_[clust.ordering_]

    plt.figure(figsize=(10, 7))
    G = gridspec.GridSpec(2, 3)
    ax1 = plt.subplot(G[0, :])
    ax2 = plt.subplot(G[1, :])

    # Reachability plot
    for klass in range(0, max(labels)+1):
        Xk = space[labels == klass]
        Rk = reachability[labels == klass]
        # ax1.plot(Xk, Rk, alpha=0.3) # Added to plot line
        ax1.plot(Xk, Rk, '.', alpha=0.3)
    ax1.plot(space[labels == -1], reachability[labels == -1], 'k.', alpha=0.3)
    ax1.set_ylabel('Reachability (epsilon distance)')
    ax1.set_title('Reachability Plot')

    # OPTICS
    ax2.plot(X[clust.labels_ == -1, 0], X[clust.labels_ == -1, 1]*y_scaling_factor, '.', markersize=4, color='lightgrey', alpha=0.3)
    for klass in range(0, max(labels)+1):
        Xk = X[clust.labels_ == klass]
        ax2.plot(Xk[:, 0], Xk[:, 1]*y_scaling_factor, '.', alpha=0.3)
        ax2.set_ylim([-10, 600])
        if Space == 'Force-CL':
            ax2.set_xlabel('Contourlength [nm]')
        else:
            ax2.set_xlabel('Distance [nm]')
        ax2.set_ylabel('Force [pN]')


    # Dict clustering keys=cluster_nr, val = x and y vals of all peaks in this cluster
    clustering = {}
    for index, point in enumerate(clusters):
        if point not in clustering:
            clustering[point] = {'x':[], 'y':[]}
        clustering[point]['x'].append(all_peaks_x[index])
        clustering[point]['y'].append(all_peaks_y[index])

    # Sorting of clusters (without noise-points) based on their x-avg
    cluster_x_avgs = []
    for i in range(len(clustering)-1):
        cluster_x_avgs.append(np.average(clustering[i]['x']))
        sorted_clust_i = np.argsort(cluster_x_avgs)

    clustering_sorted = {}
    clustering_sorted[-1] = clustering[-1]
    for index, i in enumerate(sorted_clust_i):
        clustering_sorted[index] = {}
        clustering_sorted[index]['x'] = clustering[i]['x']
        clustering_sorted[index]['y'] = clustering[i]['y']

    return clustering_sorted