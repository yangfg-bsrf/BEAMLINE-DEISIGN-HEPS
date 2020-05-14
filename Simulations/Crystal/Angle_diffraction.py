# -*- coding: utf-8 -*-
r"""
Beamline 6, wave propagation caculation: 
V.2018.12.22  
基本功能实现：
（1）改进通量 
（2）

"""
__author__ = "Fugui Yang",
__date__ = "12 Aug 2018"

import os, sys; 
import xrt.backends.raycing.materials as rmats
import numpy as np
import matplotlib.pyplot as plt

#%%%%Material properties
#反射镜材料
Pt = rmats.Material(elements=r"Pt",kind=r"mirror",rho=21.4,name=None)
Rh = rmats.Material(elements=r"Rh",kind=r"mirror",rho=12.41,name=None)
Si = rmats.Material(elements=r"Si",kind=r"mirror",rho=2.33,name=None)
Rh_mirror = rmats.Coated(coating=Rh,cThickness=500,surfaceRoughness=3,substrate=Si,name=None)
Pt_mirror = rmats.Coated(coating=Pt,cThickness=500,surfaceRoughness=3,substrate=Si,name=None)
crystalSi01 = rmats.CrystalSi(name=None,hkl=[3,1,1])#单色器材料
mBeryllium = rmats.Material('Be', rho=1.845, kind='lens')#透镜材料
material = mBeryllium

C = rmats.Material(elements=r"C",rho=2.26,name=None)
C_mirror = rmats.Coated(coating=C,cThickness=50,surfaceRoughness=50,substrate=Si,name=None)
stripeSiO2 = rmats.Material(('Si', 'O'), quantities=(1, 2), rho=2.2)

mB4C = rmats.Material(('B', 'C'), quantities=(4, 1), rho=2.28)
mW = rmats.Material('W', rho=19.3)
mL = rmats.Multilayer( tLayer= mB4C, tThickness= 19.32, bLayer= mW, bThickness=12.88, nPairs=200,substrate=Si)

E0 = np.arange(4000,30000,100)
theta = crystalSi01.get_Bragg_angle(E0)*180/3.14

plt.figure()
plt.semilogx(E0*1e-3,theta)
plt.xlabel('Energy (keV)', fontsize=14)
plt.ylabel('Bragg angle (degree)', fontsize=14)
plt.legend()
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.grid()
plt.savefig('Bragg_angle.png', dpi = 200)
plt.show()

file = open("Si311.txt","w") 
for i,j in zip(E0, theta):
     file.write("{0:.1f} {1:.6f}\n" .format(i,j))
file.close() 





