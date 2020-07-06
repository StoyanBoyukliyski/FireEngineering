# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:59:58 2020

@author: srb119
"""

import numpy as np
import matplotlib.pyplot as plt


class ParametricTimeCurveCalc():
    def __init__(self, h,w,l,ao,d1,sh1,tc1,d2,sh2,tc2,d3,sh3,tc3,d4,sh4,tc4,d5,sh5,tc5,d6,sh6,tc6,wavg):
        self.heigth = h
        self.width = w
        self.length = l
        self.AreaOpening = ao
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
        self.hwinavg = wavg
        

    def CalculateAreasB(self):
        self.AreaW1 = self.heigth*self.width
        self.AreaW2 = self.heigth*self.width
        self.AreaW3 = self.heigth*self.length
        self.AreaW4 = self.heigth*self.length
        self.AreaFloor = self.length*self.width
        self.AreaCelling = self.length*self.width
        self.AreaTotal = self.AreaW1 + self.AreaW2 + self.AreaW3 + self.AreaW4 + self.AreaFloor + self.AreaCelling
        self.bw1 = np.sqrt(self.density1*self.specheat1*self.thermalcon1)
        self.bw2 = np.sqrt(self.density2*self.specheat2*self.thermalcon2)
        self.bw3 = np.sqrt(self.density3*self.specheat3*self.thermalcon3)
        self.bw4 = np.sqrt(self.density4*self.specheat4*self.thermalcon4)
        self.bwfl = np.sqrt(self.density5*self.specheat5*self.thermalcon5)
        self.bwce = np.sqrt(self.density6*self.specheat6*self.thermalcon6)
        self.Bavg = (self.AreaW1*self.bw1 + self.AreaW2*self.bw2 + 
                     self.AreaW3*self.bw3 + self.AreaW4*self.bw4 + 
                     self.AreaFloor*self.bwfl + 
                     self.AreaCelling*self.bwce)/(self.AreaW1 + 
                                              self.AreaW2 + self.AreaW3 + self.AreaW4 + 
                                              self.AreaFloor + self.AreaCelling)
        self.O = self.AreaOpening*np.sqrt(self.hwinavg)/self.AreaTotal
        self.G = 8.41*10**8*(self.O/self.Bavg)**2
        
    def ClassifyFire(self, qfk, m, delq1, delq2, delqn,tlim, gamolam):
        self.qfk = qfk*gamolam
        self.m = m
        self.delq1 = delq1
        self.delq2 = delq2 
        self.delqn = delqn
        self.qtd = self.qfk*self.m*self.delq1*self.delq2*self.delqn*self.AreaFloor/self.AreaTotal
        self.taumax = 0.2*10**-3*self.qtd/self.O
        self.tlim =tlim
        self.tmaxh = max(self.taumax, self.tlim)
        if self.tmaxh == self.taumax:
            self.x = 1
        elif self.tmaxh == self.tlim:
            self.x = self.tlim/self.taumax
        else:
            pass
        
    def DeriveCurves(self):
        self.Tempvector = []
        self.TimeVector = []
        self.tmax = self.tmaxh*60
        self.tempmax = 20 + 1325*(1-0.324*np.exp(-0.2*self.tmax*self.G)-
                                  0.204*np.exp(-1.7*self.tmax*self.G)-
                                  0.427*np.exp(-19*self.tmax*self.G))
        for j in np.arange(0,self.tmax,10):
            self.t = 20 + 1325*(1-0.324*np.exp(-0.2*j*self.G)-
                                0.204*np.exp(-1.7*j*self.G)-0.427*np.exp(-19*j*self.G))
            self.Tempvector.append(self.t)
            self.TimeVector.append(j)
            
        for j in np.arange(self.tmax, 60, 10):
            if self.taumax*self.G <= 0.5:
                self.t = self.tempmax - 625*self.G*(j-self.x*self.taumax)
                self.Tempvector.append(self.t)
                self.TimeVector.append(j)
            elif self.taumax*self.G > 0.5 and self.taumax*self.G <= 2:
                self.t = self.tempmax - 250*self.G*(3-self.taumax*self.G)*(j-self.x*self.taumax)
                self.Tempvector.append(self.t)
                self.TimeVector.append(j)
            elif self.taumax*self.G >= 0.5:
                self.t = self.tempmax - 250*self.G*(j - self.x*self.taumax)
                self.Tempvector.append(self.t)
                self.TimeVector.append(j)
                
        plt.plot(self.TimeVector, self.Tempvector)
        plt.show()
        
        
PC = ParametricTimeCurveCalc(4.191, 6.096, 9.144, 2,900,1000, 0.25, 900, 1000, 0.25, 900, 1000, 0.25, 2560, 840,0.96, 2300,1000,1.6, 900, 1000, 0.25, 2)
PC.CalculateAreasB()
PC.ClassifyFire(511,0.8,1.1,1,1, 0.33, 0.43)
PC.DeriveCurves()
        
        