#!/usr/bin/env py -3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 17:30:34 2018

@author: SiMa

About: Main file for a simple BEM code for AE4135
"""

import numpy as np

def load_polar(filename, appendCSV=True):
    if appendCSV:
        filename += 'csv' 
    polar = np.genfromtxt(filename, delimiter=',', skip_header=2)
    
    return polar

