#! /usr/bin/env py -3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 11:13:19 2018

@author: Simon

Plotting module for the BEM model
"""

import matplotlib.pyplot as plt

def plot_data(data):
    """
    Plots data in a graph with LaTeX support
    """
    fig = plt.figure()
    ax=fig.add_axes(111)
    return

"""
def plot_lines_x(data, names, xlabel, ylabel, filename, extension = '.eps', fill=False):
    '''
    plot each line in the data matrix
    '''
    num_plots = len(names)-1
    print 'Number of lines: ', num_plots
    
    
    
    
    names_text = defchar.replace(names,'_', ' ')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    for i in range(num_plots):
        ax1.plot(data[names[0]], data[names[i+1]], label=names_text[i+1])
    if fill:
        #Find lowest values until the first vertical line, defined by nan in any other than the vline itself (cannot be the first)
        x_end = data[names[0]][np.where(np.isnan(data[names[1]]))][0]
        region = find_lowest(data, names, x_end=x_end)
        x = data[names[0]][np.where(data[names[0]]<=x_end)]
        ax1.fill_between(x,0, region, facecolor='#00FF00', alpha = 0.75)
        
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_ylim([-0.01,0.21])
    ax1.set_xlim([1.,1.3])
    ax1.legend()
    ax1.grid(b=True,alpha=0.5,linestyle='--')
    plt.show()
    
    filename = filename + extension
    
    plt.savefig(filename)
    return
"""