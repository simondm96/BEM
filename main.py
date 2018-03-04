#!/usr/bin/env py -3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 17:30:34 2018

@author: SiMa

About: Main file for a simple BEM code for AE4135
"""

import numpy as np
import matplotlib.pyplot as plt

class rotor:
    """
    A class which contains properties and a discretisation for a rotor for the 
    BEM model.
    
    Attributes:
            elements    = ndarray, contains the r-coordinates of the middle 
                          sections of the blade
            ends        = ndarray, contains the ends of each section, has length 
                          len(blades)+1
            twist       = ndarray, contains the twist of each section in radians
            pitch       = float, global blade pitch in radians
            chord       = ndarray, chord distribution of each section in meter
            num_blades  = float, number of blades
            
    """
    def __init___(self, num_blades, r_in, r_out, twist, chord, N=20):
        """
        Initialises the rotor class
        
        Input:
            num_blades  = int, number of blades
            r_in        = float, starting radius of the blade
        """
        

def twist(section, r_start, r_end):
    """
    Generates a twist distrubution
    
    Input:
        section     = ndarray, the section(s) for which the twist should be
                      calculated
                      
    Output:
        twist       = ndarray, the twist for the section(s). If sections is a
                      float, returns a float
    """
    section_norm = map_values(section,r_start, r_end, 0, 1)
    twist = 14*(1-section_norm)*np.pi/180.
    return twist


def chord(section, r_start, r_end):
    """
    Generates a chord distribution
    
    Input:
        section     = ndarray, the section(s) for which the chord should be
                      calculated
                      
    Output:
        twist       = ndarray, the chord for the section(s). If section is a
                      float, returns a float
    """
    section_norm = map_values(section,r_start, r_end, 0, 1)
    chord = (3*(1-section_norm)+1)
    return chord


def map_values(data, x_start1, x_end1, x_start2, x_end2):
    """
    Maps data with boundaries x_start1 and x_end1 to x_start2 and x_start2
    """
    return x_start2 + (data-x_start1)*(x_end1-x_start1)/(x_end2-x_start1)


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

