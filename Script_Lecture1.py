# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 20:13:49 2020

@author: srb119
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
from mpl_toolkits.mplot3d import Axes3D

class Conduction():
    def __init__(self, k, c, rho, th, Temp1, Temp0,h):
        self.k = k
        self.c = c
        self.rho = rho
        self.distance = th
        self.T0 = Temp0
        self.T1 = Temp1
        self.nmax = 100
        self.endtime= 100
        self.h = h
        self.Tfire = 500
        self.Tamb = 20
    
    def CondSteadyState(self):
        self.heatflux = self.k*(self.T1 - self.T0)/self.distance
        print(self.heatflux)
        return self.heatflux
    
    def Convection(self, Tav, Twall):
        self.Tav = Tav
        self.Twall = Twall
        self.ConvTemp= self.h * (self.Tav - self.Twall)
        print(self.ConvTemp)
    
    def CondConvSteadyState(self):
        self.aheat = np.array([self.h*self.Tfire, self.h*self.Tamb, 0])
        self.MatCoef = np.array([[1, self.h, 0],
                                 [1, 0, -self.h],
                                 [-1, self.k/self.distance, -self.k/self.distance]])
        self.unkn = np.matmul(np.linalg.inv(self.MatCoef),self.aheat)
    
        
    def CondConvInsul(self,k, Lh):
        self.kins = k
        self.Lh = Lh
        self.aheatins = np.array([self.h*self.Tfire, 0, 0, self.h*self.Tamb])
        self.MatCoefins = np.array([[1, self.h, 0,0],
                                    [-1, self.kins/self.Lh, -self.kins/self.Lh, 0],
                                    [-1, 0, self.k/self.distance, -self.k/self.distance],
                                    [-1, 0, 0, self.h]])
        self.unknins = np.matmul(np.linalg.inv(self.MatCoefins), self.aheatins)
    
    def CondTransient(self):
        self.alpha = self.k/(self.rho*self.c)
        x = np.linspace(0.01, self.distance, 50)
        t = np.linspace(1, self.endtime, 50)
        X, T= np.meshgrid(x, t)
        self.fita = (self.T1)*(1-X/self.distance)
        self.heatfluxvec = -(self.T1)/self.distance
        
        for n in range(1,self.nmax):
            self.lam = np.array(m.sqrt(self.alpha)*n*m.pi/self.distance)
            self.fita = self.fita - (2*(self.T1)/m.pi)*(1/n)*(np.exp(-T*self.lam**2)*np.sin(n*m.pi*X/self.distance))
            self.heatfluxvec = self.heatfluxvec - (2*self.T1/self.distance)*(np.exp(-T*self.lam**2)*np.cos(n*m.pi*X/self.distance))
            
            
        figure1, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize =(10,10))
        ax1 = plt.axes(projection = "3d")
        ax1.plot_surface(X, T, self.fita, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        #ax1.contour3D(X, T, fita, 100, cmap='RdGy')
        
        ax1.set_xlabel('Distance (m)')
        ax1.set_ylabel('Time (min)')
        ax1.set_zlabel('Temperature (C)')
        
        
        figure2, ax2 = plt.subplots(nrows = 1, ncols = 1, figsize =(10,10))
        ax2 = plt.axes(projection = "3d")
        ax2.plot_surface(X, T, self.heatfluxvec, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        #ax2.contour3D(X, T, heatflux, 600, cmap='RdGy')
        ax2.set_xlabel('Distance (m)')
        ax2.set_ylabel('Time (min)')
        ax2.set_zlabel('Heat Flux (W/mK)')
        print("HEATFLUX VECTOR= ", self.heatfluxvec)
        print("TIME = ", T)
        print("Temperature = ", self.fita)
        print("Distance = ", X)
        plt.show()


class Radiation():
    def __init__(self):
        self.c = 2.998*10**8
        self.h = 6.626*10**-34
        self.k = 1.381*10**-23
        self.StBolt = 5.660*10**-8

    def BlackBodyRadiation(self,eps):
        self.WVL = []
        self.Energy = []
        self.eps =eps
        for Temp in range(273, 5800, 500):
            self.WVL = []
            self.Energy = []
            for lam in np.arange(10**-7, 10**-4, 10**-9):     
                E = (self.eps*2*m.pi*(self.c**2)*(self.h)*(lam**(-5)))/(m.exp((self.c*self.h)/(lam*self.k*Temp))-1)
                self.WVL.append(lam)
                self.Energy.append(E)
                
            plt.plot(self.WVL, self.Energy, label = str(Temp) + "K" )
            plt.yscale("log")
            plt.xscale("log")
            plt.legend()
        
        plt.show()
        
    def TotalEnergy(self,eps):
        self.eps = eps
        self.Etotal = []
        self.Temperature = []
        self.eps = eps
        for Temp in range(273,5800, 500): 
            ent = self.StBolt*self.eps*Temp**4
            self.Etotal.append(ent)
            self.Temperature.append(Temp)
        plt.plot(self.Temperature, self.Etotal)
        plt.ylabel("Heat flux")
        plt.xlabel("Teperature (K)")
        
        plt.show()
    
    def HeatTransfer(self, AC, BD, AD, BC, CD, eps1,eps2, Te, Tr):
        self.Fi = np.sqrt((AC + BD - (AD+BC))**2)/(2*CD)
        self.eps1 = eps1
        self.eps2 = eps2
        self.Te = Te
        self.Tr = Tr
        self.he = self.Fi*self.eps1*self.eps2*self.StBolt*(self.Te**4-self.Tr**4)
        
    
        
        
Inst1 = Conduction(45.8,1400,20,0.05, 40, 20,8)
Inst1.CondConvSteadyState()
print(Inst1.unkn)
Inst1.CondConvInsul(0.15, 0.01)
print(Inst1.unknins)