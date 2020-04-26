# -*- coding: utf-8 -*-
"""
2019-9-7 杨福桂：从SPECTRA读取二维功率密度分布，画图

"""
import numpy as np
import matplotlib.pyplot as plt
import os, sys; sys.path.append(os.path.join('..', '..', '..'))  
from mpl_toolkits.mplot3d import Axes3D


#SPECTRA软件计算结果分析：一定口径大小
filename = 'Power_10keV_normal'
Import_data = np.loadtxt(filename + '.dta', skiprows=2)

thetax = np.unique(Import_data[:,0])
thetay = np.unique(Import_data[:,1])

X, Y = np.meshgrid(thetax, thetay)
I1 = np.reshape(Import_data[:,2],(len(thetay), len(thetax)))

fig = plt.figure()
ax1 = plt.axes(projection='3d')
ax1.plot_surface(X*1e3, Y*1e3, I1, alpha=0.7,cmap='rainbow') 
ax1.set_xlabel('X ($\mu$rad)', fontsize=14)
ax1.set_ylabel('Y ($\mu$rad)', fontsize=14)
ax1.set_zlabel('P.Density (kW/mrad^2)', fontsize=14)
#plt.xticks(fontsize=13)
#plt.yticks(fontsize=13)

plt.savefig(filename + '_Power density.png', dpi = 200)

filename = 'Power_10keV_taper'
Import_data = np.loadtxt(filename + '.dta', skiprows=2)

thetax = np.unique(Import_data[:,0])
thetay = np.unique(Import_data[:,1])

X, Y = np.meshgrid(thetax, thetay)
I2 = np.reshape(Import_data[:,2],(len(thetay), len(thetax)))

fig = plt.figure()
ax1 = plt.axes(projection='3d')
ax1.plot_surface(X*1e3, Y*1e3, I2, alpha=0.7,cmap='rainbow') 
ax1.set_xlabel('X ($\mu$rad)', fontsize=14)
ax1.set_ylabel('Y ($\mu$rad)', fontsize=14)
ax1.set_zlabel('P.Density (kW/mrad^2)', fontsize=14)
#plt.xticks(fontsize=13)
#plt.yticks(fontsize=13)

plt.savefig(filename + '_Power density.png', dpi = 200)

fig = plt.figure()
ax1 = plt.axes(projection='3d')
ax1.plot_surface(X*1e3, Y*1e3, I2-I1, alpha=0.7,cmap='rainbow') 
ax1.set_xlabel('X ($\mu$rad)', fontsize=14)
ax1.set_ylabel('Y ($\mu$rad)', fontsize=14)
ax1.set_zlabel('P.Density (kW/mrad^2)', fontsize=14)
#plt.xticks(fontsize=13)
#plt.yticks(fontsize=13)

plt.savefig('_Power density.png', dpi = 200)