#!/usr/bin/env py -3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 17:30:34 2018

@author: SiMa

About: Main file for a simple BEM code for AE4135
"""

import numpy as np
import matplotlib.pyplot as plt

def load_polar(filename, appendCSV=True):
    """
    Loads a polar from a file with title and headers
    
    Input:
            filename  = str, name of the file which contains the polar
            AppendCSV = Bool, defines if '.csv' should be added to filename, defaults to True
          
    Output:
            polar     =  numpy array containing AoA, Cl, Cd, Cm
    """
    if appendCSV:
        filename += '.csv' 
    polar = np.genfromtxt(filename, delimiter=',', skip_header=2)
    
    return polar

