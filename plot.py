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