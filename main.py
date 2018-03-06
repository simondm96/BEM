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
        self.twist = twist(self.elements, r_in, r_out)
        self.chord = chord(self.elements, r_in, r_out)
        self.num_blades = num_blades
        
        
def middle_vals(data):
    return np.delete((np.roll(data,1)-data)/2+data,0)


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
    return x_start2 + (data-x_start1)*(x_end2-x_start2)/(x_end1-x_start1)


def polarvalues(alpha):
    pol = xlrd.open_workbook("polar_DU95W180.xlsx")
    pol = pol.sheet_by_index(0)
    for i in range(2, 62):
        if pol.cell_value(i,0)<= alpha <=pol.cell_value(i+1,0):
            xp = []
            data = []
            xp.append(pol.cell_value(i,0))
            xp.append(pol.cell_value(i+1,0))
            for n in range(3):
                fp = []
                fp.append(pol.cell_value(i,n+1))
                fp.append(pol.cell_value(i+1,n+1))
                b = np.interp(alpha, xp, fp)
                data.append(b)
    return data

def liftdragcalc(twist, u_inf, a, aprime, omega, r, chord, rho):
    
    v_ax = u_inf*(1-a)
    v_tan = omega*r*(1-aprime)
    
    ratio = np.v_ax / v_tan
    
    phi = np.arctan(ratio)
    
    alpha = abs(phi - twist)
    
    vp = np.sqrt(v_ax**2 + v_tan**2)
    
    lift = 0.5*chord*rho*(vp**2)*polarvalues(alpha)[1]
    drag = 0.5*chord*rho*(vp**2)*polarvalues(alpha)[2]
    
    f_azim = lift*(v_ax/vp) -drag*(v_tan/vp)
    f_axial = lift*(v_tan/vp) +drag*(v_ax/vp)
    
    return alpha, lift, drag, f_azim, f_axial


def inductioncalc(f_azim, f_axial, nblades, rho, u_inf, r, deltar, lamda, R):
    #axial induction
    A_a = 2*np.pi()*r*deltar
    
    CT = (f_axial *nblades *deltar)/(0.5*rho*(u_inf**2)*A_a)
    a = 0.5*(1-np.sqrt(1-CT))
    #azimuthal induction
    aprime = (f_azim*nblades)/(2*rho*(2*np.pi()*r)*(u_inf**2)*(1-a)*lamda*(r/R))
    
    return CT, a, aprime
    