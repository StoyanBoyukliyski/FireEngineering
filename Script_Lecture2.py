# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:59:58 2020

@author: srb119
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

class ParametricTimeCurveCalc():
    
    '''
    __init__ function triggering global class parameters upon creating an instance.
    Assign all parameters related to the compartment and the temperature time curve.
    
    '''
    def __init__(self, h,w,l,ao,d1,sh1,tc1,d2,sh2,tc2,d3,sh3,tc3,d4,sh4,tc4,d5,sh5,tc5,d6,sh6,tc6,wavg):
        #Define geometry of the compartment including openings
        self.heigth = h
        self.width = w
        self.length = l
        self.AreaOpening = ao
            
        #Define the material properties of the walls and roof and floor of the compartment
        self.density1 = d1
        self.specheat1 = sh1
        self.thermalcon1 = tc1
        self.density2 = d2
        self.specheat2 = sh2
        self.thermalcon2 = tc2
        self.density3 = d3
        self.specheat3 = sh3
        self.thermalcon3 = tc3
        self.density4 = d4
        self.specheat4 = sh4
        self.thermalcon4 = tc4
        self.density5 = d5
        self.specheat5 = sh5
        self.thermalcon5 = tc5
        self.density6 = d6
        self.specheat6 = sh6
        self.thermalcon6 = tc6
        
        
        #the average height of all of the openings in the compartment
        self.hwinavg = wavg
        
        #Define the path where files are stored
        self.path = r"C:\Users\StoyanBoyukliyski"
        
        #Define the time steps in subsequent analysis
        self.step = 0.001
        
        
        #Assign the areas of the compartment
        self.AreaW1 = self.heigth*self.width
        self.AreaW2 = self.heigth*self.width
        self.AreaW3 = self.heigth*self.length
        self.AreaW4 = self.heigth*self.length
        self.AreaFloor = self.length*self.width
        self.AreaCelling = self.length*self.width
        self.AreaTotal = self.AreaW1 + self.AreaW2 + self.AreaW3 + self.AreaW4 + self.AreaFloor + self.AreaCelling
        
        #Assign the b values used in subsequent analysis
        self.bw1 = np.sqrt(self.density1*self.specheat1*self.thermalcon1)
        self.bw2 = np.sqrt(self.density2*self.specheat2*self.thermalcon2)
        self.bw3 = np.sqrt(self.density3*self.specheat3*self.thermalcon3)
        self.bw4 = np.sqrt(self.density4*self.specheat4*self.thermalcon4)
        self.bwfl = np.sqrt(self.density5*self.specheat5*self.thermalcon5)
        self.bwce = np.sqrt(self.density6*self.specheat6*self.thermalcon6)
        self.Bavg = (self.AreaW1*self.bw1 + self.AreaW2*self.bw2 + 
                     self.AreaW3*self.bw3 + self.AreaW4*self.bw4 - self.AreaOpening*self.bw4 +
                     self.AreaFloor*self.bwfl + 
                     self.AreaCelling*self.bwce)/(self.AreaW1 + 
                                              self.AreaW2 + self.AreaW3 + self.AreaW4 - self.AreaOpening +
                                              self.AreaFloor + self.AreaCelling)
        
        #Assign O and G in the analysis
        self.O = self.AreaOpening*np.sqrt(self.hwinavg)/self.AreaTotal
        if self.O< 0.02:
            self.O = 0.02
        elif self.O > 0.2:
            self.O = 0.2
        self.G = 8.41*10**8*(self.O/self.Bavg)**2
        
    '''
    Classify the fire in the compartment based on most likely loading scenarios taken
    from EN 1991-1-2 and all required factors taken into account related to fire
    Establish limiting parameters based on fuel controlled fire.
    '''
        
    def ClassifyFire(self, qfk, m, delq1, delq2, delqn,tlim, gamolam):
        #qfk is the fire loading multiplied by the factor specified in the exercise   
        self.qfk = qfk*gamolam
        #m factor is defaulted to 0.8 unless otherwise specified
        self.m = m
        #delq1, delq2, delqn are factors related to the risk related to the structure
        self.delq1 = delq1
        self.delq2 = delq2 
        self.delqn = delqn
        
        
        #Governing equation for the applied fire loading used in subsequent equations
        self.qtd = self.qfk*self.m*self.delq1*self.delq2*self.delqn*self.AreaFloor/self.AreaTotal
        
        
        #taumax is used for ventilation controlled fires
        self.taumax = 0.2*10**-3*self.qtd/self.O
        
        #tlim is based on building type
        self.tlim =tlim
        
        #Decision about type of fire control
        self.tmax = max(self.taumax, self.tlim)
        if self.tmax == self.taumax:
            self.x = 1
            print("VentialtionControlled")
        elif self.tmax == self.tlim:
            self.x = self.tlim/self.taumax
            print("Fuel controlled")
        else:
            pass
        self.Olim = 0.1*10**-3*self.qtd/self.tlim
        self.Glim = 8.41*10**8*(self.Olim/self.Bavg)**2
        
        
    '''
    Build the actual Temperature-Time curve 
    after defining all of the parameters.
    and store the values in a CSV file.
    '''

    def DeriveCurves(self):
        self.Tempvector = []
        self.TimeVector = []
        #Start by defining the maximum temperature if the fire is fuel controlled
        if self.tmax == self.tlim:
            self.tempmax = 20 + 1325*(1-0.324*np.exp(-0.2*self.tmax*self.Glim)-
                                      0.204*np.exp(-1.7*self.tmax*self.Glim)-
                                      0.427*np.exp(-19*self.tmax*self.Glim))
            
            #Define curve points up to the maximum temperature, time point
            for j in np.arange(0,self.tmax,self.step):
                self.t = 20 + 1325*(1-0.324*np.exp(-0.2*j*self.Glim)-
                                    0.204*np.exp(-1.7*j*self.Glim)-
                                    0.427*np.exp(-19*j*self.Glim))
                self.Tempvector.append(self.t)
                self.TimeVector.append(j*60)
        
        #Start by defining the maximum temperature if the fire is ventillation controlled
        elif self.tmax == self.taumax:
            self.tempmax = 20 + 1325*(1-0.324*np.exp(-0.2*self.tmax*self.G)-
                                      0.204*np.exp(-1.7*self.tmax*self.G)-
                                      0.427*np.exp(-19*self.tmax*self.G))
            
            #Define curve points up to the maximum temperature, time point
            for j in np.arange(0,self.tmax,self.step):
                self.t = 20 + 1325*(1-0.324*np.exp(-0.2*j*self.G)-
                                    0.204*np.exp(-1.7*j*self.G)-
                                    0.427*np.exp(-19*j*self.G))
                self.Tempvector.append(self.t)
                self.TimeVector.append(j*60)
                
        #Define the decending gradient from maximum point downwards
        for j in np.arange(self.tmax, 2, self.step):
            
            #Define condition is taumax*G <=0.5       
            if self.taumax*self.G <= 0.5:
                self.t = self.tempmax - 625*self.G*(j-self.x*self.taumax)
                if self.t < 0: 
                    self.t = 0
                self.Tempvector.append(self.t)
                self.TimeVector.append(j*60)
            #Define condition is taumax*G <=0.5 and taumax*G>2      
            elif self.taumax*self.G > 0.5 and self.taumax*self.G <= 2:
                self.t = self.tempmax - 250*self.G*(3-self.taumax*self.G)*(j-self.x*self.taumax)
                if self.t < 0: 
                    self.t = 0
                self.Tempvector.append(self.t)
                self.TimeVector.append(j*60)
            #Define condition is taumax*G >=2
            elif self.taumax*self.G >= 2:
                self.t = self.tempmax - 250*self.G*(j - self.x*self.taumax)
                if self.t < 0: 
                    self.t = 0
                self.Tempvector.append(self.t)
                self.TimeVector.append(j*60)
                
        self.PanTemp = pd.DataFrame({"Temperature": self.Tempvector,
                                     "Time" : self.TimeVector})
    
        self.PanTemp.to_csv(self.path + "\ParametricTimeCurve.csv",index =False)
    
    '''
    Build the Fire Curve not related
    to any of the previous parameters
    Use only a simple function
    '''
    
    def StandardFire(self):
        self.Temperature =[]
        self.time = []
        #Define the temperature-time equation 
        for j in np.arange(0, 120, self.step*60):
            temp= 20+ 345*np.log10(8*j+1)
            self.Temperature.append(temp-20)
            self.time.append(j)
            
        self.FirePD = pd.DataFrame({"Time" : self.time,
                             "Temperature" : self.Temperature})
        self.FirePD.to_csv(self.path + "\StandardFire.csv",index =False)
        
    '''
    Create plots to be inserted in the final report
    '''
    
    def Plotters(self):
        #Joint figure showing the Fire Curve and the Temperature Parametric Time curve together
        fig1, (ax1,ax2) = plt.subplots(2,1)
        fig1.set_figheight(10)
        fig1.set_figwidth(10)
        fig1.tight_layout(pad=3.0)
        ax1.grid(linestyle = "--")
        ax1.plot(self.TimeVector, self.Tempvector, "r-", linewidth = 4)
        ax1.set_title("Parametric Time Curve", {'fontsize': 12,
                                                'fontweight' : 12,
                                                'verticalalignment': 'baseline'})
        ax1.set_ylabel("Temperature (C)""\N{DEGREE SIGN}""",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad=20)
        ax1.set_xlabel("Time (min)",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad = 10)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(which='minor', width=1, length = 3)
        ax1.tick_params(which='major', width=1, length = 6)
        
        ax2.plot(self.time, self.Temperature, linewidth = 4)
        ax2.grid(linestyle = "--")
        ax2.set_title("Standard Fire Curve",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'})
        ax2.set_ylabel("Temperature (C)""\N{DEGREE SIGN}""",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad = 10)
        ax2.set_xlabel("Time (min)",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad = 15)
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(which='minor', width=1, length = 3)
        ax2.tick_params(which='major', width=1, length = 6)
        
        
        
        
        #Figure showing the Fire Curve
        fig2, ax1 = plt.subplots(1,1, linewidth = 2)
        fig2.set_figheight(10)
        fig2.set_figwidth(10)
        ax1.plot(self.TimeVector, self.Tempvector, "r-", linewidth = 4)
        ax1.grid(linestyle = "--")
        ax1.set_title("Parametric Time Curve",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'})
        ax1.set_ylabel("Temperature (C)""\N{DEGREE SIGN}""",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad = 10)
        ax1.set_xlabel("Time (min)",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad = 10)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(which='minor', width=1, length = 3)
        ax1.tick_params(which='major', width=1, length = 6)
        
        
        
        #Figure showing the Temperature Parametric Time
        fig3, ax1 = plt.subplots(1,1, linewidth = 2)
        fig3.set_figheight(10)
        fig3.set_figwidth(10)
        ax1.plot(self.time, self.Temperature,  "r-", linewidth = 4)
        ax1.grid(linestyle = "--")
        ax1.set_title("Parametric Time Curve",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'})
        ax1.set_ylabel("Temperature (C)""\N{DEGREE SIGN}""",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad =10)
        ax1.set_xlabel("Time (min)",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad = 10)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(which='minor', width=1, length = 3)
        ax1.tick_params(which='major', width=1, length = 6)
        
            
        #One plot joint figure showing the Fire Curve and the Temperature Parametric Time curve together
        fig4, ax1 = plt.subplots(1,1)
        fig4.set_figheight(10)
        fig4.set_figwidth(10)
        ax1.plot(self.time, self.Temperature, linewidth = 4)
        ax1.grid(linestyle = "--")
        ax1.plot(self.TimeVector, self.Tempvector, "r-", linewidth = 4)
        ax1.set_title("Parametric Time Curve",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'})
        ax1.set_ylabel("Temperature (C)""\N{DEGREE SIGN}""",{'fontsize': 12,
                                               'fontweight' :12,
                                               'verticalalignment': 'baseline'},labelpad = 10)
        ax1.set_xlabel("Time (min)",{'fontsize': 12,
                                               'fontweight' : 12,
                                               'verticalalignment': 'baseline'},labelpad = 10)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(which='minor', width=1, length = 3)
        ax1.tick_params(which='major', width=1, length = 6)
        
        plt.show()