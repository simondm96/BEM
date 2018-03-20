#!/usr/bin/env py -3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 17:30:34 2018

@author: Simon, Zyanya, Wessel

About: Main file for a simple BEM code for AE4135
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Bookman'], 'size':10})
rc('text', usetex=True)


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
            pitch       = float, pitch of the rotor in radians
            
    """
    def __init__(self, num_blades, r_in, r_out, twist, chord, pitch, spacing = 'norm', N=100):
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
        if spacing != 'norm':
            self.ends = 0.5*(1-np.cos(np.linspace(r_in, r_out, num=N)))
        else:
            self.ends = np.linspace(r_in, r_out, num=N)
        self.elements = middle_vals(self.ends)
        self.mu = map_values(self.elements, r_in, r_out, r_in/r_out, 1)
        self.twist = twist(self.mu)
        self.chord = chord(self.mu)
        self.num_blades = num_blades
        self.pitch = pitch
    
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
        """
        Loads a lift and drag polar from a file with 'filename'
        """
        self.polar = np.genfromtxt(filename, delimiter = ",", skip_header = 2)
    
    def polarvalues(self, alpha):
        """
        Linearly interpolates the values from the polar for a given angle of attack
        """
        cl = np.interp(alpha, self.polar[:,0], self.polar[:,1])
        cd = np.interp(alpha, self.polar[:,0], self.polar[:,2])
        return cl, cd
    
    def liftdragcalc(self, u_inf, a, aprime, TSR, rho):
        """
        Calculates the lift and drag from polar data fr given conditions
        """
        omega = TSR * u_inf /self.ends[-1]
        v_ax = u_inf*(1-a)
        v_tan = omega*self.elements*(1+aprime)
        
        ratio = v_ax / v_tan
        
        phi = np.arctan(ratio)
        
        alpha = (phi - self.twist + self.pitch)*180/np.pi #Angle of attack in degrees
        
        vp = np.sqrt(v_ax**2 + v_tan**2)
        polar = self.polarvalues(alpha)
        lift = 0.5*self.chord*rho*(vp**2)*polar[0]
        drag = 0.5*self.chord*rho*(vp**2)*polar[1]
        
        f_azim = lift*(v_ax/vp) -drag*(v_tan/vp)
        f_axial = lift*(v_tan/vp) +drag*(v_ax/vp)
        return alpha, lift, drag, f_azim, f_axial, phi
    
    def inductioncalc(self, rho, u_inf, a, aprime, TSR):
        """
        Calculates the axial and azimuthal induction for given conditions
        """
        #axial induction
        out = self.liftdragcalc(u_inf, a, aprime, TSR, rho)
        f_azim, f_axial = out[3], out[4]
        
        deltar = np.diff(self.ends)
        A_a = 2*np.pi*self.elements*deltar
        
        CT = (f_axial *self.num_blades*deltar)/(0.5*rho*(u_inf**2)*A_a)
        a = heavy_loading_thrust(CT)
        a *= 1/self.tip_root_correction(a, TSR)
        
        #azimuthal induction
        aprime = (f_azim*self.num_blades)/(2*rho*(2*np.pi*self.elements)*(u_inf**2)*(1-a)*TSR*(self.mu))
        aprime *= 1/self.tip_root_correction(a, TSR)
        
        CP = 4*a*((1-a)**2)
        
        return CT, CP, a, aprime, out
    
    
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
    return np.diff(data)/2+np.delete(data, -1)


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




def run(TSR):
    #Input parameters
    u_inf = 10.
    N_blades = 3
    hubrR = 0.2
    R = 50.
    pitch = 2*np.pi/180
    a = 0.3 #starting value
    aprime = 0.0 #starting value
    rho = 1.225
    #loop parameters
    n_max = 50
    n=0
    diff_a = 1
    diff_aprime = 1
    #convergence history
    a_list = []
    aprime_list = []
    # initialising the rotor class
    rotor_BEM = rotor(N_blades, hubrR*R, R, twist, chord, pitch)
    rotor_BEM.loadpolar("polar_DU95W180.csv")
    
    while (diff_a>0.0001 and diff_aprime>0.0001) and n<n_max:
        #save old induction factors
        a_old = a
        aprime_old = aprime
        #Calculate the solution for the given induction factors
        CT, CP, a, aprime, out = rotor_BEM.inductioncalc(rho, u_inf, a, aprime, TSR)
        #Correct the new induction factors
        a = 0.25*a + 0.75*a_old
        aprime = 0.25*aprime + 0.75*aprime_old
        #Find the difference with the old induction factors
        diffa = np.abs(a_old-a)
        diffaprime = np.abs(aprime_old-aprime)
        #Use the maximum difference as a measure of the convergence
        diff_a = np.amax(diffa)
        diff_aprime = np.amax(diffaprime)
        a_list.append(diff_a)
        aprime_list.append(diff_aprime)
        #Next iteration!        
        n+=1
    
    
    print("Iterations:",n) #uncomment on python 3.x
    conv = np.vstack((np.array(a_list), np.array(aprime_list)))
    return CT, CP, a, aprime, out, conv
    
    
    
def plotdata_single(xdata, ydata, xname, yname, fname, extension=".eps"):
    """
    Plots data with given names and saves it to a filename with extension
    """
    legend = ["TSR = 10", "TSR = 8", "TSR = 6"]
    fig = plt.figure(figsize=(6,4.5))
    ax = fig.add_subplot(111)
    for i in range(3):
        ax.plot(xdata, ydata[:,i], label = legend[i])
    ax.set_xlim([0, xdata[-1]*1.05])
    ax.set_ylim([0, np.amax(ydata)*1.05])
    ax.set_xlabel(xname)
    ax.set_ylabel(yname)
    ax.legend()
    ax.grid(b=True,alpha=0.5,linestyle='--')
    plt.savefig(fname+extension)
    plt.show()

def plotdata_double(xdata, ydata, xname, yname, fname, extension=".eps"):
    """
    Plots data with given names and saves it to a filename with extension, use two datasets with 3 TSRs
    """
    legend = ["TSR10", "TSR8", "TSR6"]
    legend2 = ["a", "a'"]
    fname_l = []
    for leg in legend:
        fname_l.append(fname+leg)

        
    fig1 = plt.figure(figsize=(6,4.5))
    ax1 = fig1.add_subplot(111)
    for i in range(2):
        ax1.plot(xdata, ydata[:,i], label = legend2[i])
    ax1.set_xlim([0, xdata[-1]*1.05])
    ax1.set_ylim([0, np.amax(ydata[:,0:2])*1.05])
    ax1.set_xlabel(xname)
    ax1.set_ylabel(yname)
    ax1.legend()
    ax1.grid(b=True,alpha=0.5,linestyle='--')
    plt.savefig(fname_l[0]+extension)
    
    fig2 = plt.figure(figsize=(6,4.5))
    ax2 = fig2.add_subplot(111)
    for i in range(2):
        ax2.plot(xdata, ydata[:,i+2], label = legend2[i])
    
    ax2.set_xlim([0, xdata[-1]*1.05])
    ax2.set_ylim([0, np.amax(ydata[:,2:4])*1.05])
    ax2.set_xlabel(xname)
    ax2.set_ylabel(yname)
    ax2.legend()
    ax2.grid(b=True,alpha=0.5,linestyle='--')
    plt.savefig(fname_l[1]+extension)
    
    
    fig3 = plt.figure(figsize=(6,4.5))
    
    ax3 = fig3.add_subplot(111)
    for i in range(2):
        ax3.plot(xdata, ydata[:,i+4], label = legend2[i])
    
    ax3.set_xlim([0, xdata[-1]*1.05])
    ax3.set_ylim([0, np.amax(ydata[:,4:6])*1.05])
    ax3.set_xlabel(xname)
    ax3.set_ylabel(yname)
    ax3.legend()
    ax3.grid(b=True,alpha=0.5,linestyle='--')
    plt.savefig(fname[2]+extension)
    plt.show()
    
