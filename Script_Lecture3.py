# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as m


class MaterialCarbonSteel():
    def __init__(self):
        self.fp0 = 350*10**6
        self.fy0 = 490*10**6
        self.Ea0 = 210*10**9
        self.epsp0 = self.fp0/self.Ea0
        self.epsy0 = 0.02
        self.epst0 = 0.15
        self.epsu0 = 0.2
        self.path = r"C:\Users\StoyanBoyukliyski"
        

    def ReductionFactor(self):        
        self.c = (self.fy0-self.fp0)**2/((self.epsy0 -self.epsp0)*self.Ea0 - 2*(self.fy0-self.fp0))
        self.a = np.sqrt((self.epsy0-self.epsp0)*(self.epsy0-self.epsp0+self.c/self.Ea0))
        self.b = np.sqrt(self.c*(self.epsy0-self.epsp0)*self.Ea0 + self.c**2)
        
        self.Eq = []
        
        for eps in np.arange(0,self.epsp0,0.0001):
            self.eqn1 = eps*self.Ea0
            self.Eq.append(self.eqn1)
        
        for eps in np.arange(self.epsp0 + 0.0001 ,self.epsy0,0.0001):
            self.eqn2 = self.fp0 - self.c + (self.b/self.a)*(self.a**2 - (self.epsy0 - eps)**2)**0.5
            self.Eq.append(self.eqn2)
            
        for eps in np.arange(self.epsy0 + 0.0001, self.epst0, 0.0001):
            self.eqn3 = self.fy0
            self.Eq.append(self.eqn3)
            
        for eps in np.arange(self.epst0,self.epsu0, 0.0001):
            self.eqn4 = self.fy0*(1-(eps-self.epst0)/(self.epsu0-self.epst0))
            self.Eq.append(self.eqn4)
            
        plt.plot(np.arange(0, self.epsu0, 0.0001), self.Eq)
        plt.plot()
        
        self.data = {"ky0 = fy0/fy": [1.0, 1.0, 1.0, 1.0, 1.0, 0.780, 0.470, 0.230, 0.110, 0.060, 0.040, 0.020, 0.000],
                "kp0 = fp0/fy": [1.0, 1.0, 0.807, 0.613, 0.420, 0.360, 0.180, 0.075, 0.050, 0.0375, 0.0250, 0.0125, 0.000],
                "ke0 = Ea0/Ea": [1.0, 1.0, 0.900, 0.800, 0.700, 0.600, 0.310, 0.130, 0.090, 0.0675, 0.0450, 0.0225, 0.000]}
        
        self.Temperatures = [20, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
        
        
        self.ReductionFactors = pd.DataFrame(self.data, index = self.Temperatures)
        
        self.ReductionFactors.to_csv(self.path + "\ReductionFactorsCarbonSteel.csv")
            
        self.Eq1 = []
        self.Equation = []
        self.ElasticityReduction = []
        
        
        self.eps = np.linspace(0, self.epsu0, 300)
        self.t = self.ReductionFactors.index
        self.T, self.EPS = np.meshgrid(self.t, self.eps)
        for j in self.t:
            self.fy = self.fy0*self.ReductionFactors.loc[j]["ky0 = fy0/fy"]
            self.fp = self.fp0*self.ReductionFactors.loc[j]["kp0 = fp0/fy"]
            self.Ea = self.Ea0*self.ReductionFactors.loc[j]["ke0 = Ea0/Ea"]
            self.ElasticityReduction.append(self.Ea)
            
            self.c = (self.fy-self.fp)**2/((self.epsy0 -self.epsp0)*self.Ea - 2*(self.fy-self.fp))
            self.a = np.sqrt((self.epsy0-self.epsp0)*(self.epsy0-self.epsp0+self.c/self.Ea))
            self.b = np.sqrt(self.c*(self.epsy0-self.epsp0)*self.Ea + self.c**2)
            
            self.Eq1 = []
            
            for i in self.eps:
                if 0< i <= self.epsp0:
                    self.Eq1.append(i*self.Ea)
                elif self.epsp0 < i <= self.epsy0:
                    self.Eq1.append(self.fp - self.c + (self.b/self.a)*(self.a**2 - (self.epsy0 - i)**2)**0.5)
                elif self.epsy0 < i <= self.epst0:
                    self.Eq1.append(self.fy)
                elif self.epst0 < i <= self.epsu0:
                    self.Eq1.append(self.fy*(1-(i-self.epst0)/(self.epsu0-self.epst0)))
                elif i < self.epsu0:
                    self.Eq1.append(0)
                else:
                    pass
            self.Equation.append(self.Eq1)
    
        figure1, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize =(10,10))
        ax1 = plt.axes(projection = "3d")
        ax1.plot_surface(self.EPS, self.T, np.transpose(self.Equation), rstride=1, cstride=1, cmap='RdGy', edgecolor='none')
        "ax1.contour3D(EPS, T, np.transpose(Equation), 100, cmap='RdGy')"
        
        ax1.set_xlabel('Strain ()')
        ax1.set_ylabel('Temperature (C)')
        ax1.set_zlabel('Stress')
        
        plt.show()
        
    def ThermalElongation(self):
        self.Elong = []
        self.Temperature = []
        for j in range(20, 750, 1):    
            dlol = 1.2*10**-5*j + 0.4*10**-8*j**2 - 2.416*10**-4
            self.Elong.append(dlol)
            self.Temperature.append(j)
            
        for j in range(750, 860, 1):
            dlol = 1.1*10**-2
            self.Elong.append(dlol)
            self.Temperature.append(j)
        
        for j in range(860, 1200, 1):
            dlol = 2*10**-5*j - 6.2*10**-3
            self.Elong.append(dlol)
            self.Temperature.append(j)
        
        ThermEl = pd.DataFrame({"Thermal Elongation": self.Elong,
                                       "Temperature ": self.Temperature})
            
        
        ThermEl.to_csv(self.path + "\ThermalElongationCarbonSteel.csv",index =False)
        plt.plot(self.Temperature, self.Elong)
        plt.show()
        
    def CustomSHC(self,temp):
        ca = 0
        la =0
        if temp< 600:
            ca = 425 + 7.73*10**(-1)*temp-1.69*10**(-3)*temp**(2) + 2.22*10**(-6)*temp**(3)
        elif temp>=600 and temp < 735:
            ca = 666 + 13002/(738-temp)
        elif temp>=735 and temp < 900:
            ca = 545 + 17820/(temp-731)
        elif temp >=900 and temp <1200:
            ca = 650
        else:
            pass
        if temp <800:
            la = 54 -3.33*10**-2*temp
        elif temp >=800 and temp < 1200:
            la =27.3
        else:
            pass
        return [ca, la]
        
    def SpecificHeat(self):
        self.SpecHeat = []
        self.Temperature = []
        for j in range(20,1200,5):
            if j< 600:
                ca = 425 + 7.73*10**(-1)*j-1.69*10**(-3)*j**(2) + 2.22*10**(-6)*j**(3)
            elif j>=600 and j < 735:
                ca = 666 + 13002/(738-j)
            elif j>=735 and j < 900:
                ca = 545 + 17820/(j-731)
            elif j>=900 and j <1200:
                ca = 650
            else:
                pass
            self.SpecHeat.append(ca)
            self.Temperature.append(j-20)
        
        SpecificHeatDF = pd.DataFrame({"Specific Heat": self.SpecHeat,
                                       "Temperature ": self.Temperature})
        
        SpecificHeatDF.to_csv(self.path + "\SpecificHeatCarbonSteel.csv",index =False)
        plt.plot(self.Temperature, self.SpecHeat)
        plt.show()
    
    def ThermalConductivity(self):
        self.ThCon = []
        self.Temperature = []
        for j in range(20,1200, 10):
            if j <800:
                la = 54 -3.33*10**-2*j
            elif j >=800 and j < 1200:
                la =27.3
            else:
                pass
            
            self.ThCon.append(la)
            self.Temperature.append(j-20)
        ThermalConductivityDF = pd.DataFrame({"Thermal Conductivity" : self.ThCon,
                                              "Temperature": self.Temperature})
        ThermalConductivityDF.to_csv(self.path + "\ThermalConductivityCarbonSteel.csv",index =False)
        plt.plot(self.Temperature, self.ThCon)
        plt.show()
        
        
class MaterialStainlessSteel():
    def __init__(self):
        self.density= 7850
        self.Ea = 205*10**9
        self.fy = 350*10**6
        self.fu = 490*10**6
        self.epscfita = 0.02
        self.epsufita = 0.2
    
    def ReductionFactors(self):
        pass
    
    def ThermalElongation(self):
        pass
    
    def SpecificHeat(self):
        pass
    
    def ThermalConductivity(self):
        pass
    
        
class MaterialConcrete():
    def __init__(self):
        pass
    
    def ReductionFactors(self):
        pass
    
    def ThermalElongation(self):
        pass
    
    def SpecificHeat(self,i1):
        MC = i1
        if MC == 0:
            self.SpecH = []
            self.Temperature = []
            for j in range(20,1200, 10):
                if j < 100:
                    cp= 900
                elif j>=100 and j < 120:
                    cp = 900+(j-100) + 0
                elif j>= 120 and j< 200:
                    cp = 900 +(j-100) + 0 - 0*(j-120)/80
                elif j >= 200 and j < 400:
                    cp = 1000 + (j-200)/2
                elif j>=400 and j < 1200:
                    cp = 1100
                else:
                    pass
                self.SpecH.append(cp)
                self.Temperature.append(j-20)
            self.SpecificHeatDF = pd.DataFrame({"Specific Conductivity" : self.SpecH,
                                           "Temperature": self.Temperature})
            self.SpecificHeatDF.to_csv(self.path + "\SpecificHeatConcrete00.csv",index =False)
            
            plt.plot(self.Temperature, self.SpecH)
            plt.show()
            
        elif MC == 1.5:
            self.SpecH = []
            self.Temperature = []
            for j in range(20,1200, 10):
                if j < 100:
                    cp= 900
                elif j>=100 and j < 120:
                    cp = 900+(j-100) + 570
                elif j>= 120 and j< 200:
                    cp = 900 +(j-100) + 570 - 570*(j-120)/80
                elif j >= 200 and j < 400:
                    cp = 1000 + (j-200)/2
                elif j>=400 and j < 1200:
                    cp = 1100
                else:
                    pass
                self.SpecH.append(cp)
                self.Temperature.append(j-20)
            self.SpecificHeatDF = pd.DataFrame({"Specific Conductivity" : self.SpecH,
                                           "Temperature": self.Temperature})
            self.SpecificHeatDF.to_csv(self.path + "\SpecificHeatConcrete15.csv",index =False)
            
            plt.plot(self.Temperature, self.SpecH)
            plt.show()
        
        elif MC == 3.0:
            self.SpecH = []
            self.Temperature = []
            for j in range(20,1200, 10):
                if j < 100:
                    cp= 900
                elif j>=100 and j < 120:
                    cp = 900+(j-100) + 1020
                elif j>= 120 and j< 200:
                    cp = 900 +(j-100) + 1020 - 1020*(j-120)/80
                elif j >= 200 and j < 400:
                    cp = 1000 + (j-200)/2
                elif j>=400 and j < 1200:
                    cp = 1100
                else:
                    pass
                self.SpecH.append(cp)
                self.Temperature.append(j-20)
            self.SpecificHeatDF = pd.DataFrame({"Specific Conductivity" : self.SpecH,
                                           "Temperature": self.Temperature})
            self.SpecificHeatDF.to_csv(self.path + "\SpecificHeatConcrete30.csv",index =False)
            
            plt.plot(self.Temperature, self.SpecH)
            plt.show()
            
    
    def ThermalConductivity(self):
        self.UPPER = []
        self.LOWER = []
        self.Temperature = []
        for j in range(20,1200,10):
            upper = 2 - 0.2451*(j/100) +0.0107*(j/100)**2
            lower = 1.36-0.136*(j/100) + 0.0057*(j/100)**2
            self.UPPER.append(upper)
            self.LOWER.append(lower)
            self.Temperature.append(j-20)
        self.ThermalConductivityDF = pd.DataFrame({"Thermal Conductivity UPPER Siliceous Concrete" : self.UPPER,
                                              "Thermal Conductivity LOWER Calcareous Concrete" : self.LOWER,
                                              "Temperature": self.Temperature})
        self.ThermalConductivityDF.to_csv(self.path + "\ThermalConductivityConcrete.csv",index =False)
        plt.plot(self.Temperature, self.UPPER)
        plt.plot(self.Temperature, self.LOWER)
        plt.show()
            
    def Density(self):
        self.DensityPar = []
        self.Temperature = []
        densityreal = 2300
        for j in range(20,1200,10):
            if j < 115:
                den = densityreal
            elif j >= 115 and j < 200: 
                den = densityreal*(1-0.02*(j-115)/85)
            elif j >= 200 and j <400:
                den = densityreal *(0.98 -0.03*(j-200)/200)
            elif j >= 400 and j <1200:
                den = densityreal*(0.95-0.07*(j-400)/800)
            else: 
                pass
            self.DensityPar.append(den)
            self.Temperature.append(j-20)
        self.ThermalConductivityDF = pd.DataFrame({"Density Concrete" : self.DensityPar,
                                              "Temperature": self.Temperature})
        self.ThermalConductivityDF.to_csv(self.path + "\DensityConcrete.csv",index =False)
        plt.plot(self.Temperature, self.DensityPar)
        plt.show()
        
        
class ReinforcementSteel():
    def __init__(self):
        pass
    
    def ReductionFactor(self):
        pass
    
    def ThermalElongation(self):
        pass