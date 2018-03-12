#!/usr/bin/env py -3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 17:30:34 2018

@author: Simon, Zyanya, Wessel

About: Main file for a simple BEM code for AE4135
"""




import numpy as np
import matplotlib.pyplot as plt
import xlrd
import xlwt

pol = xlrd.open_workbook("polar_DU95W180.xlsx")
pol = pol.sheet_by_index(0)

class rotor:
    """
    A class which contains properties and a discretisation for a rotor for the 
    BEM model.
    
    Attributes:
            elements    = ndarray, contains the r-coordinates of the middle 
                          sections of the blade
            ends        = ndarray, contains the ends of each section, has length 
                          len(blades)+1
            mu          = ndarray, contains the elements normalised to the radius of the rotor
            twist       = ndarray, contains the twist of each section in radians
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
    
    def tip_root_correction(self, a, TSR):
        """
        Applies Prandtl's tip and hub correction to the forces
        """
        mu_r = self.ends[0]/self.ends[-1]
        cst = np.sqrt(1+TSR**2*self.mu**2/((1-a)**2))
        exp_tip = -self.num_blades/2*(1-self.mu)/self.mu*cst
        exp_root = -self.num_blades/2*(self.mu-mu_r)/self.mu*cst
        f_tip = 2/np.pi*np.arccos(np.exp(exp_tip))
        f_root = 2/np.pi*np.arccos(np.exp(exp_root))
        return f_tip*f_root
    
    def loadpolar(self, filename):
        pol = xlrd.open_workbook("polar_DU95W180.xlsx")
        self.polar = pol.sheet_by_index(0)
    
    def polarvalues(self, alpha):
        for i in range(2, 62):
            if self.polar.cell_value(i,0)<= alpha <=self.polar.cell_value(i+1,0):
                xp = []
                data = []
                xp.append(self.polar.cell_value(i,0))
                xp.append(self.polar.cell_value(i+1,0))
                for n in range(3):
                    fp = []
                    fp.append(self.polar.cell_value(i,n+1))
                    fp.append(self.polar.cell_value(i+1,n+1))
                    b = np.interp(alpha, xp, fp)
                    data.append(b)
        return data
    
    @staticmethod
    def heavy_loading_induction(a):
        """
        Applies Prandtl's correction for heavily loaded rotors based on induction factors
        """
        CT1 = 1.816
    
        CT = np.where(a>=1-np.sqrt(CT1)/2, CT1-4*(np.sqrt(CT1)-1)*(1-a), 4*a*(1-a))
        return CT

    def heavy_loading_thrust(CT):
        """
        Applies Prandtl's correction for heavily loaded rotors based on thrust coefficient
        """
        CT1 = 1.816
        CT2 = 2*np.sqrt(CT1) - CT1
        
        a = np.where(CT>=CT2, 1+(CT-CT1)/(4*np.sqrt(CT1)-4), 1/2.-np.sqrt(1-CT)/2)
    
        return a
        

def middle_vals(data):
    return np.diff(data)+np.delete(data, -1)


def twist(section):
    """
    Generates a twist distrubution
    
    Input:
        section     = ndarray, the section(s) for which the twist should be
                      calculated, normalised to the radius of the rotor
                      
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
                      calculated, normalised to the radius of the rotor
                      
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



def liftdragcalc(twist, u_inf, a, aprime, omega, r, chord, rho):
    
    v_ax = u_inf*(1-a)
    v_tan = omega*r*(1-aprime)
    
    ratio = v_ax / v_tan
    
    phi = np.arctan(ratio)
    
    alpha = abs(phi - twist)
    
    vp = np.sqrt(v_ax**2 + v_tan**2)
    polar = polarvalues(alpha)
    lift = 0.5*chord*rho*(vp**2)*polar[1]
    drag = 0.5*chord*rho*(vp**2)*polar[2]
    
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
    CP = 4*a*((1-a)**2)
    
    return CT, CP, a, aprime

