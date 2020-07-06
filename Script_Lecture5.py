
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 19:56:32 2020

@author: srb119
"""

import Script_Lecture2 as SL2
import Script_Lecture3 as SL3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

''' 

Calling the ParametricTimeCurveCalcs() class with the following parameters:
ParametricTimeCurveCalc(h,w,l,ao,d1,sh1,tc1,d2,sh2,tc2,d3,sh3,tc3,d4,sh4,tc4,d5,sh5,tc5,d6,sh6,tc6,wavg)

Where:
    h = height of the compartment
    w = width of the compartment
    l = length of the compartment
    ao = area of the openings (calculated on hand depending on existing windows and doors in compartment)
    d1 = density of one of the walls
    sh1 = specific heat of one of the walls
    tc1 = thermal capacity of one of the walls
    d2 = density of one of the walls
    sh2 = specific heat of one of the walls
    tc2 = thermal capacity of one of the walls
    d3 = density of one of the walls
    sh3 = specific heat of one of the walls
    tc3 = thermal capacity of one of the walls
    d4 = density of one of the walls
    sh4 = specific heat of one of the walls
    tc4 = thermal capacity of one of the walls
    d5 = density of one of the walls
    sh5 = specific heat of one of the walls
    tc5 = thermal capacity of one of the walls
    d6 = density of one of the walls
    sh6 = specific heat of one of the walls
    tc6 = thermal capacity of one of the walls
    wavg = average hight of the openings
    
'''

'''
In EXERCISE:
    
    PC1:
        h = height of the compartmnt                    = 4.191 m 
        w = width of the compartment                    = 9.144 m
        l = length of the compartment                   = 6.096 m
        ao = area of the openings                       = 1 m x 2 m = 2 m^2
        d1 = density of one of the walls                = 900 kg/m^3
        sh1 = specific heat of one of the walls         = 0.25 W/mK
        tc1 = thermal conductivity of one of the walls  = 1000 J/kgK
        d2 = density of one of the walls                = 900 kg/m^3 
        sh2 = specific heat of one of the walls         = 0.25 J/kgK
        tc2 = thermal conductivity of one of the walls  = 1000 W/mK
        d3 = density of one of the walls                = 900 kg/m^3 
        sh3 = specific heat of one of the walls         = 0.25 J/kgK
        tc3 = thermal conductivity of one of the walls  = 1000 W/mK
        d4 = density of one of the walls                = 2560 kg/m^3 
        sh4 = specific heat of one of the walls         = 0.96 J/kgK
        tc4 = thermal conductivity of one of the walls  = 1000 W/mK
        d5 = density of one of the walls                = 2300 kg/m^3
        sh5 = specific heat of one of the walls         = 1000 W/mK
        tc5 = thermal conductivity of one of the walls  = 1.6 J/kgK
        d6 = density of one of the walls                = 2300 kg/m^3
        sh6 = specific heat of one of the walls         = 1000 W/mK
        tc6 = thermal conductivity of one of the walls  = 1.6 J/kgK
        wavg = average hight of the openings            = 2 m

'''
PC1 = SL2.ParametricTimeCurveCalc(4.191, 9.144, 6.096, 2, 900, 0.25, 1000, 900, 0.25, 1000, 900, 0.25, 1000, 2560, 0.96, 1000, 2300, 1000, 1.6, 2300, 1000, 1.6, 2)

'''
Aditional Parameters to quantify the fire:
    ClassifyFire(qfk, m, delq1, delq2, delqn,tlim, gamolam)
    qfk = tabulated value that needs to be looked up from Annex E 1991-1-2
    m = 0.8 for cellulosic fuels
    delq1 = tabulated value that needs to be looked up from Annex E 1991-1-2
    delq2 = tabulated value that needs to be looked up from Annex E 1991-1-2
    delqn = tabulated value that needs to be looked up from Annex E 1991-1-2
    tlim = tabulated value that needs to be looked up from Annex E 1991-1-2
    gamolam = tabulated value that needs to be looked up from Annex E 1991-1-2
'''


'''
IN EXERCISE:
    
    PC1:
        qfk = 511 for the 80% fractile for an office building
        m = 0.8 for cellulosic fuels
        delq1 = 1.5 for a building in between 25m^2 and 250m^2
        delq2 = 1.0 for an office building
        delqn = 1.0 taken as one as I am not engaging 
                in making detailed assumptions about fire services
        tlim = 
        gamolam = 0.8 from the exercise description
        
'''

PC1.ClassifyFire(511, 0.8, 1.5, 1.0, 1.0, 0.33, 0.8)
PC1.DeriveCurves()
PC1.StandardFire()
PC1.Plotters()


MTC = SL3.MaterialCarbonSteel()

'''
PC2 = SL2.ParametricTimeCurveCalc(2.6,6.4, 3.2, 4.31, 900, 1000, 0.25, 900, 1000, 0.25, 900, 1000, 0.25, 900, 1000, 0.25, 2300, 1000, 1.6, 2300, 1000, 1.6, 1.55)
PC2.CalculateAreasB()
PC2.ClassifyFire(377,0.8,1,1,1, 0.33, 0.43)
PC2.DeriveCurves()
PC2.StandardFire()
'''

class FireDesign():
    def __init__(self, dens = 7890, Fi = 1, epsf = 1.0, epsm = 0.8):
        self.Fi = Fi
        self.epsf= epsf
        self.epsm = epsm
        self.dens = dens
        self.boltz = 5.669*10**-8
        
        self.tw = 0.005842
        self.tf = 0.008509
        self.h = 0.348
        self.w = 0.127
        self.sigmayield = 355*10**6
        self.sigmault = 470*10**6
        self.E = 205*10**9
        self.Ixx = 8.28*10**-5
        
        self.Am = (self.h-2*self.tf)*self.tw+2*self.w*self.tf
        self.P = 2*self.h+self.w + 2*(self.w-self.tw)
        self.Pb = 2*self.h+self.w
        
        self.ksh = 0.9*(self.Pb/self.P)
        
        self.AppliedMoment = 161.7*10**3
        
        
    def UnprotectedSection(self, Time, Temperature):
        self.YieldStr = [1.0]
        self.Elasticity = [1.0]
        self.UltimateStr = [1.0]
        self.RealTemp = Temperature
        self.RealTime = Time
        self.Temperature = [20]
        self.SpecHeatCust = [MTC.CustomSHC(self.Temperature[-1])[1]]
        self.Conductivity = [MTC.CustomSHC(self.Temperature[-1])[0]]
        self.dTemp = [0]
        self.TimeDif = [(self.RealTime[1]-self.RealTime[0])*60]
        for j in range(len(self.RealTime)-1):
            self.SpecHeatCust.append(MTC.CustomSHC(self.Temperature[j])[1])
            self.Conductivity.append(MTC.CustomSHC(self.Temperature[j])[0])
            self.YieldStr.append(self.YStrengthLoss(self.Temperature[j]))
            self.UltimateStr.append(self.UlStrengthLoss(self.Temperature[j]))
            self.Elasticity.append(self.ElasticLoss(self.Temperature[j]))
            self.TimeDif.append((-self.RealTime[j] + self.RealTime[j+1])*60)
            self.dTemp.append(self.UnprotectedCurve(self.RealTemp[j], self.Temperature[-1],
                                                    (-self.RealTime[j]*60 + self.RealTime[j+1]*60),
                                                    MTC.CustomSHC(self.RealTemp[j])[0]))
            
            
            
            self.Temperature.append(self.Temperature[-1]+self.dTemp[-1])
            
        self.GlobalSim = pd.DataFrame({"Time": self.RealTime,
                                       "Time Difference": self.TimeDif,
                                       "Temperature Amb": self.RealTemp,
                                       "Temperature Comp": self.Temperature,
                                       "Temperature Dif": self.dTemp,
                                       "Specific Heat": self.SpecHeatCust,
                                       "Conductivity": self.Conductivity,
                                       "Yield Strength Loss": self.YieldStr,
                                       "Ultimate Strength Loss": self.UltimateStr,
                                       "Elasticity Loss": self.Elasticity})
        self.GlobalSim.to_csv("Exercise22.csv")
        
        self.fig, (ax1,ax2,ax3) = plt.subplots(1,3)
        self.fig.set_figheight(8)
        self.fig.set_figwidth(12)
        Criteria1 = [self.AppliedMoment,self.AppliedMoment]
        Criteria2= [90,90]
        time = [0,120]
        stress = [min(min(Criteria1),min([self.MaximumStressYield()*j for j in self.YieldStr]))*0.8, 
                  max(max(Criteria1),max([self.MaximumStressUlt()*j for j in self.UltimateStr]))*1.2]
        ax3.plot(time, Criteria1, "r-", linewidth =3, label = "Applied Moment")
        ax3.plot(Criteria2, stress, "g-.", linewidth = 1, label = "90 minute-criteria")
        
        ax3.plot(self.RealTime, [self.MaximumStressUlt()*j for j in self.UltimateStr], linewidth = 3, label = "Ultimate Resistance")
        ax3.plot(self.RealTime, [self.MaximumStressYield()*j for j in self.YieldStr], linewidth = 3, label = "Elastic Resistance")
        
        
        ax2.plot(self.RealTime,self.YieldStr, linewidth = 3, label = "Yield Strength")
        ax2.plot(self.RealTime,self.UltimateStr, linewidth = 3, label = "Ultimate Strength ")
        ax2.plot(self.RealTime,self.Elasticity, linewidth = 3, label = "Young's Modulus")
        
        
        ax1.plot(self.RealTime, self.RealTemp, linewidth =3, label = "Temperature-Time Curve")
        ax1.plot(self.RealTime, self.Temperature, linewidth = 3, label = "Component Temperature")
        
        
        ax1.grid(linestyle = "--")
        ax1.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax1.set_ylabel("Temperature (C)""\N{DEGREE SIGN}""",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax1.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(which='minor', width=1, length = 3)
        ax1.tick_params(which='major', width=1, length = 6)
        ax1.legend(loc = "upper right", prop = {"size": 9})
        
        
        ax2.grid(linestyle = "--")
        ax2.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax2.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax2.set_ylabel("Capacity as a ratio from ambient (/)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(which='minor', width=1, length = 3)
        ax2.tick_params(which='major', width=1, length = 6)
        ax2.legend(loc = "upper right", prop = {"size": 9})
        
        ax3.grid(linestyle = "--")
        ax3.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax3.set_ylabel("Moment (kN/m)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax3.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax3.xaxis.set_minor_locator(AutoMinorLocator())
        ax3.tick_params(which='minor', width=1, length = 3)
        ax3.tick_params(which='major', width=1, length = 6)
        ax3.legend(loc = "upper right", prop = {"size": 9})
        
        ax3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        self.fig.savefig("Temperature Time Curve Non-Insulation.png")
        plt.show()
    
    def ProtectedSection(self, Time, Temperature):
        self.RealTempPr = Temperature
        self.RealTimePr = np.array(Time)*60
        
        
        self.YieldStrPr = [1.0]
        self.ElasticityPr = [1.0]
        self.UltimateStrPr = [1.0]
        self.TemperaturePr = [20]
        
        
        self.SpecHeatCustPr = [MTC.CustomSHC(self.TemperaturePr[-1])[0]]
        self.ConductivityPr = [MTC.CustomSHC(self.TemperaturePr[-1])[1]]
        self.dTempPr = [0]
        
        self.TimeDifPr = [(self.RealTimePr[1]-self.RealTimePr[0])]
        for j in range(len(self.RealTimePr)-1):
            self.YieldStrPr.append(self.YStrengthLoss(self.TemperaturePr[j]))
            self.UltimateStrPr.append(self.UlStrengthLoss(self.TemperaturePr[j]))
            self.ElasticityPr.append(self.ElasticLoss(self.TemperaturePr[j]))
            
            self.SpecHeatCustPr.append(MTC.CustomSHC(self.TemperaturePr[j])[0])
            self.ConductivityPr.append(MTC.CustomSHC(self.TemperaturePr[j])[1])
            
            self.TimeDifPr.append(self.RealTimePr[j+1]-self.RealTimePr[j])
            self.dTempPr.append(self.ProtectedCurve(self.RealTempPr[j], self.TemperaturePr[j], 
                                                    self.TimeDifPr[j], MTC.CustomSHC(self.RealTempPr[j])[0]))
            
            self.TemperaturePr.append(self.TemperaturePr[-1]+self.dTempPr[-1])
            
        self.GlobalSim = pd.DataFrame({"Time": self.RealTimePr,
                                       "Time Difference": self.TimeDifPr,
                                       "Temperature Amb": self.RealTempPr,
                                       "Temperature Comp": self.TemperaturePr,
                                       "Temperature Dif": self.dTempPr,
                                       "Specific Heat": self.SpecHeatCustPr,
                                       "Conductivity": self.ConductivityPr,
                                       "Yield Strength Loss": self.YieldStrPr,
                                       "Ultimate Strength Loss": self.UltimateStrPr,
                                       "Elasticity Loss": self.ElasticityPr})
        self.GlobalSim.to_csv("Exercise21" + ".csv")
        
        self.fig2, (ax1,ax2,ax3) = plt.subplots(1, 3)
        self.fig2.set_figheight(8)
        self.fig2.set_figwidth(12)
        Criteria1 = [self.AppliedMoment,self.AppliedMoment]
        Criteria2= [90,90]
        time = [0,120]
        stress = [min(min(Criteria1),min([self.MaximumStressYield()*j for j in self.YieldStrPr]))*0.8, 
                  max(max(Criteria1),max([self.MaximumStressUlt()*j for j in self.UltimateStrPr]))*1.2]
        ax3.plot(time, Criteria1, "r-.", linewidth = 1)
        ax3.plot(Criteria2, stress, "r-.", linewidth = 1)
        
        self.RealTimePr = np.array(self.RealTimePr)/60
        Criteria1 = [self.AppliedMoment,self.AppliedMoment]
        Criteria2= [90,90]
        time = [0,120]
        stress = [min(min(Criteria1),min([self.MaximumStressYield()*j for j in self.YieldStrPr]))*0.8, 
                  max(max(Criteria1),max([self.MaximumStressUlt()*j for j in self.UltimateStrPr]))*1.2]
        ax3.plot(time, Criteria1, "r-", linewidth =3, label = "Applied Moment")
        ax3.plot(Criteria2, stress, "g-.", linewidth = 1, label = "90 minute-criteria")
        
        ax3.plot(self.RealTimePr, [self.MaximumStressUlt()*j for j in self.UltimateStrPr], linewidth = 3, label = "Ultimate Resistance")
        ax3.plot(self.RealTimePr, [self.MaximumStressYield()*j for j in self.YieldStrPr], linewidth = 3, label = "Elastic Resistance")
        
        
        ax2.plot(self.RealTimePr,self.YieldStrPr, linewidth = 3, label = "Yield Strength")
        ax2.plot(self.RealTimePr,self.UltimateStrPr, linewidth = 3, label = "Ultimate Strength ")
        ax2.plot(self.RealTimePr,self.ElasticityPr, linewidth = 3, label = "Young's Modulus")
        
        
        ax1.plot(self.RealTimePr, self.RealTempPr, linewidth =3, label = "Temperature-Time Curve")
        ax1.plot(self.RealTimePr, self.TemperaturePr, linewidth = 3, label = "Component Temperature")
        
        
        ax1.grid(linestyle = "--")
        ax1.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax1.set_ylabel("Temperature (C)""\N{DEGREE SIGN}""",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax1.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(which='minor', width=1, length = 3)
        ax1.tick_params(which='major', width=1, length = 6)
        ax1.legend(loc = "upper right", prop = {"size": 7})
        
        
        ax2.grid(linestyle = "--")
        ax2.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax2.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax2.set_ylabel("Capacity as a ratio from ambient (/)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(which='minor', width=1, length = 3)
        ax2.tick_params(which='major', width=1, length = 6)
        ax2.legend(loc = "upper right", prop = {"size": 7})
        
        ax3.grid(linestyle = "--")
        ax3.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax3.set_ylabel("Moment (kN/m)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax3.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax3.xaxis.set_minor_locator(AutoMinorLocator())
        ax3.tick_params(which='minor', width=1, length = 3)
        ax3.tick_params(which='major', width=1, length = 6)
        ax3.legend(loc = "upper right", prop = {"size": 7})
        
        ax3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        self.fig2.savefig("Temperature Time Curve Paint Insulation.png")
        plt.show()
        
    def ProtectedSectionNNeg(self, Time, Temperature):
        self.YieldStrPrn = [1.0]
        self.ElasticityPrn = [1.0]
        self.UltimateStrPrn = [1.0]
        self.RealTempPrn = Temperature
        self.RealTimePrn = Time
        self.TemperaturePrn = [20]
        self.dTfire = [self.RealTempPrn[j+1]-self.RealTempPrn[j] for j in range(len(self.RealTempPrn)-1)]
        self.dTfire.append(0)
        self.SpecHeatCustPrn = [MTC.CustomSHC(self.TemperaturePrn[-1])[1]]
        self.ConductivityPrn = [MTC.CustomSHC(self.TemperaturePrn[-1])[0]]
        self.dTempPrn = [0]
        self.TimeDifPrn = [(self.RealTimePrn[1]-self.RealTimePrn[0])*60]
        for j in range(len(self.RealTimePrn)-1):
            n = 20
            self.SpecHeatCustPrn.append(MTC.CustomSHC(self.TemperaturePrn[j])[1])
            self.ConductivityPrn.append(MTC.CustomSHC(self.TemperaturePrn[j])[0])
            self.YieldStrPrn.append(self.YStrengthLoss(self.TemperaturePrn[j]))
            self.UltimateStrPrn.append(self.UlStrengthLoss(self.TemperaturePrn[j]))
            self.ElasticityPrn.append(self.ElasticLoss(self.TemperaturePrn[j]))
            self.TimeDifPrn.append((-self.RealTimePrn[j] + self.RealTimePrn[j+1])*60)
            self.dTempPrn.append(self.ProtectedCurveNNeg(self.RealTempPrn[j], self.TemperaturePrn[-1],
                                                    (-self.RealTimePrn[j]*60 + self.RealTimePrn[j+1]*60),
                                                    MTC.CustomSHC(self.RealTempPrn[j])[0],self.dTfire[j]))
            
            
            
            for i in self.dTempPrn:
                n = n + i    
            self.TemperaturePrn.append(n)
            
        self.GlobalSim = pd.DataFrame({"Time": self.RealTimePrn,
                                       "Time Difference": self.TimeDifPrn,
                                       "Temperature Amb": self.RealTempPrn,
                                       "Temperature Comp": self.TemperaturePrn,
                                       "Temperature Dif": self.dTempPrn,
                                       "Specific Heat": self.SpecHeatCustPrn,
                                       "Conductivity": self.ConductivityPrn,
                                       "Yield Strength Loss": self.YieldStrPrn,
                                       "Ultimate Strength Loss": self.UltimateStrPrn,
                                       "Elasticity Loss": self.ElasticityPrn})
        self.GlobalSim.to_csv("Exercise23" + ".csv")
        
        self.fig3, (ax1,ax2,ax3) = plt.subplots(1, 3)
        self.fig3.set_figheight(8)
        self.fig3.set_figwidth(12)
        self.fig3.tight_layout(pad=3.0)
        Criteria1 = [self.AppliedMoment,self.AppliedMoment]
        Criteria2= [90,90]
        time = [0,120]
        stress = [min(min(Criteria1),min([self.MaximumStressYield()*j for j in self.YieldStrPrn]))*0.8, 
                  max(max(Criteria1),max([self.MaximumStressUlt()*j for j in self.UltimateStrPrn]))*1.2]
        ax3.plot(time, Criteria1, "r-", linewidth =3, label = "Applied Moment")
        ax3.plot(Criteria2, stress, "g-.", linewidth = 1, label = "90 minute-criteria")
        
        ax3.plot(self.RealTimePrn, [self.MaximumStressUlt()*j for j in self.UltimateStrPrn], linewidth = 3, label = "Ultimate Resistance")
        ax3.plot(self.RealTimePrn, [self.MaximumStressYield()*j for j in self.YieldStrPrn], linewidth = 3, label = "Elastic Resistance")
        
        
        ax2.plot(self.RealTimePrn,self.YieldStrPrn, linewidth = 3, label = "Yield Strength")
        ax2.plot(self.RealTimePrn,self.UltimateStrPrn, linewidth = 3, label = "Ultimate Strength ")
        ax2.plot(self.RealTimePrn,self.ElasticityPrn, linewidth = 3, label = "Young's Modulus")
        
        
        ax1.plot(self.RealTimePrn, self.RealTempPrn, linewidth =3, label = "Temperature-Time Curve")
        ax1.plot(self.RealTimePrn, self.TemperaturePrn, linewidth = 3, label = "Component Temperature")
        
        
        ax1.grid(linestyle = "--")
        ax1.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax1.set_ylabel("Temperature (C)""\N{DEGREE SIGN}""",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax1.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(which='minor', width=1, length = 3)
        ax1.tick_params(which='major', width=1, length = 6)
        ax1.legend(loc = "upper right", prop = {"size": 7})
        
        
        ax2.grid(linestyle = "--")
        ax2.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax2.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax2.set_ylabel("Capacity as a ratio from ambient (/)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(which='minor', width=1, length = 3)
        ax2.tick_params(which='major', width=1, length = 6)
        ax2.legend(loc = "upper right", prop = {"size": 7})
        
        ax3.grid(linestyle = "--")
        ax3.set_title("Parametric Time Curve", {'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'})
        ax3.set_ylabel("Moment (kN/m)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad=20)
        ax3.set_xlabel("Time (min)",{'fontsize': 12, 'fontweight' : 12, 'verticalalignment': 'baseline'},labelpad = 10)
        ax3.xaxis.set_minor_locator(AutoMinorLocator())
        ax3.tick_params(which='minor', width=1, length = 3)
        ax3.tick_params(which='major', width=1, length = 6)
        ax3.legend(loc = "upper right", prop = {"size": 7})
        
        ax3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tight_layout()
        self.fig3.savefig("Temperature Time Curve Insulation with HC.png")
        plt.show()
        
    def UnprotectedCurve(self, Fitag, Fitam, dt, ca):
        self.alphac = 50
        dFitam= self.ksh*(self.P/self.Am)*(dt/(ca*self.dens))*(self.alphac*(Fitag-Fitam)+ 
                          self.Fi*self.epsf*self.epsm*self.boltz*((Fitag+273)**4 - (Fitam+273)**4))
        
        return dFitam

    def ProtectedCurve(self, Fitag, Fitam, dt, ca):
        self.alphacins = 0.15
        self.dp = 0.02
        self.Ap = self.Am
        dFitam = ((self.P/self.Ap)*(self.alphacins*(Fitag - Fitam)*dt))/(self.dp*ca*self.dens)
        
        return dFitam
    
    
    def ProtectedCurveNNeg(self, Fitag, Fitam, dt, ca, dFitag):
        self.alphap = 0.15
        self.cp = 1200
        self.densp = 800
        self.dpn = 0.02
        self.Ap = self.w*self.h
        self.f = self.cp*self.densp*self.dpn*self.Pb/(ca*self.dens*self.Ap)
        dFitam = (self.Pb/self.Ap)*self.alphap*(Fitag-Fitam)*dt/(self.dpn*ca*self.dens*(1+
                 self.f/3))-(np.exp(self.f/10)-1)*dFitag
        if dFitam <0 and dFitag > 0:
            dFitam = 0
        else:
            pass
        return dFitam
    
    
    def UlStrengthLoss(self,T):
        self.dataU = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.780, 0.470, 0.230, 0.110, 0.060, 0.040, 0.020, 0.000, 0.000]
        self.Temperatures = [-20, 20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 2000]
        for m in range(len(self.dataU)-1):
            if T > self.Temperatures[m] and T <= self.Temperatures[m+1]:
                self.kyt1 = self.dataU[m]*(1 - (T-self.Temperatures[m])/(self.Temperatures[m+1]-self.Temperatures[m]))
                self.kyt2 = self.dataU[m+1]*((T - self.Temperatures[m])/(self.Temperatures[m+1]-self.Temperatures[m]))
                self.kyt = self.kyt1 + self.kyt2
            else:
                pass
        return self.kyt
        
    def ElasticLoss(self,T):
        self.dataE = [1.0, 1.0, 1.0, 0.900, 0.800, 0.700, 0.600, 0.310, 0.130, 0.090, 0.0675, 0.0450, 0.0225, 0.000, 0.0]
        self.Temperatures = [-20, 20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 2000]
        for m in range(len(self.dataE)-1):
            if T > self.Temperatures[m] and T <= self.Temperatures[m+1]:
                self.kyt1 = self.dataE[m]*(1 - (T-self.Temperatures[m])/(self.Temperatures[m+1]-self.Temperatures[m]))
                self.kyt2 = self.dataE[m+1]*((T - self.Temperatures[m])/(self.Temperatures[m+1]-self.Temperatures[m]))
                self.kyt = self.kyt1 + self.kyt2
            else:
                pass
        return self.kyt
    
    
    def YStrengthLoss(self, T):
        self.datayi = [1.0, 1.0, 1.0, 0.807, 0.613, 0.420, 0.360, 0.180, 0.075, 0.050, 0.0375, 0.0250, 0.0125, 0.000, 0.0]
        self.Temperatures = [-20, 20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 2000]
        for m in range(len(self.datayi)-1):
            if T > self.Temperatures[m] and T <= self.Temperatures[m+1]:
                self.kyt1 = self.datayi[m]*(1 - (T-self.Temperatures[m])/(self.Temperatures[m+1]-self.Temperatures[m]))
                self.kyt2 = self.datayi[m+1]*((T - self.Temperatures[m])/(self.Temperatures[m+1]-self.Temperatures[m]))
                self.kyt = self.kyt1 + self.kyt2
            else:
                pass
        return self.kyt
    
    
    def MaximumStressUlt(self):
        self.d = ((self.tw*self.w)*(self.h/2 - self.tw/2) + ((self.h/2 - self.tw)*1/2)*self.tf*(self.h/2-self.tw))/(self.tw*self.w+self.tf*(self.h/2-self.tw))
        self.F = self.sigmault*((self.tw*self.w)+self.tf*(self.h/2-self.tf))
        self.Mrdult = 2*self.F*self.d
        return self.Mrdult
    
    
    def MaximumStressYield(self):
        self.Mrdel = self.sigmayield*self.Ixx/(self.h/2)
        return self.Mrdel

        
FD1 = FireDesign()
FD1.UnprotectedSection(PC1.TimeVector, PC1.Tempvector)
FD1.ProtectedSection(PC1.TimeVector, PC1.Tempvector)
FD1.ProtectedSectionNNeg(PC1.TimeVector, PC1.Tempvector)


FD2 = FireDesign()
FD2.UnprotectedSection(PC1.time, PC1.Temperature)
FD2.ProtectedSection(PC1.time, PC1.Temperature)
FD2.ProtectedSectionNNeg(PC1.time, PC1.Temperature)