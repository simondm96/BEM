#!/usr/bin/env py -3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 17:30:34 2018

@author: SiMa

About: Main file for a simple BEM code for AE4135
"""

import numpy as np
import matplotlib.pyplot as plt
import xlrd
import xlwt

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
            pitch       = float, global blade pitch in radians (should this be included?)
            chord       = ndarray, chord distribution of each section in meter
            num_blades  = float, number of blades
            
    """
    def __init__(self, num_blades, r_in, r_out, twist, chord, N=20):
        """
        Initialises the rotor class
        
        Input:
            num_blades  = int, number of blades
            r_in        = float, inner radius of the blade
            r_out       = float, outer radius of the blade
            chord       = function, defines the chord distribution
            twist       = function, defines the twist distribution
            N           = int, number of ends. The number of elements is N-1
            
        """
        self.ends = np.linspace(r_in, r_out, num=N)
        self.elements = middle_vals(self.ends)
        self.mu = map_values(self.elements,r_in, r_out, 0.2, 1)
        self.twist = twist(self.mu)
        self.chord = chord(self.mu)
        self.num_blades = num_blades
        
        
def middle_vals(data):
    return np.delete((np.roll(data,1)-data)/2+data,0)


def twist(section):
    """
    Generates a twist distrubution
    
    Input:
        section     = ndarray, the section(s) for which the twist should be
                      calculated
                      
    Output:
        twist       = ndarray, the twist for the section(s). If sections is a
                      float, returns a float
    """
    return 14*(1-section)*np.pi/180.


def chord(section):
    """
    Generates a chord distribution
    
    Input:
        section     = ndarray, the section(s) for which the chord should be
                      calculated
                      
    Output:
        twist       = ndarray, the chord for the section(s). If section is a
                      float, returns a float
    """
    return (3*(1-section)+1)


def map_values(data, x_start1, x_end1, x_start2, x_end2):
    """
    Maps data with boundaries x_start1 and x_end1 to x_start2 and x_start2
    """
    return x_start2 + (data-x_start1)*(x_end2-x_start2)/(x_end1-x_start1)

def tip_root_correction(rotor, a, TSR):
    """
    Applies Prandtl's tip and hub correction to the forces
    """
    mu_r = rotor.ends[0]/rotor.ends[-1]
    cst = np.sqrt(1+TSR**2*rotor.mu**2/((1-a)**2))
    exp_tip = -rotor.num_blades/2*(1-rotor.mu)/rotor.mu*cst
    exp_root = -rotor.num_blades/2*(rotor.mu-mu_r)/rotor.mu*cst
    f_tip = 2/np.pi*np.arccos(np.exp(exp_tip))
    f_root = 2/np.pi*np.arccos(np.exp(exp_root))
    f_corr = f_tip*f_root
    return f_corr

def heavy_loading_induction(a):
    """
    Applies Prandtl's correction for heavily loaded rotors based on induction factors
    """
    CT1 = 1.816

    
    if a>1-np.sqrt(CT1)/2:
        CT = CT1-4*(np.sqrt(CT1)-1)(1-a)
    else:
        CT = 4*a(1-a)
    return CT

def heavy_loading_thrust(CT):
    """
    Applies Prandtl's correction for heavily loaded rotors based on thrust coefficient
    """
    CT1 = 1.816
    CT2 = 2*np.sqrt(CT1) - CT1
    
    if CT>CT2:
        a = 1+(CT-CT1)/(4*np.sqrt(CT)-4)
    else:
        a = 1/2.-np.sqrt(1-CT)/2
    
    return a

def polarvalues(alpha):
    pol = xlrd.open_workbook("polar_DU95W180.xlsx")
    pol = pol.sheet_by_index(0)
    for i in range(2, 62):
        if pol.cell_value(i,0)<= alpha <=pol.cell_value(i+1,0):
            data = []
            for n in range(3):
                b = pol.cell_value(i, n)
                data.append(b)
    return data


