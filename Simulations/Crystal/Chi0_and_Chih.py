# -*- coding: utf-8 -*-
"""
__author__ = "QJ,Jia"
__date__ = "2020-03-27"

"""
import os, sys; sys.path.append(os.path.join('..', '..', '..'))  # analysis:ignore
import numpy as np
import xrt.backends.raycing.materials as rmats
import xrt.backends.raycing.physconsts as physconsts


Si111_130K = rmats.CrystalSi(tK=130, factDW=0.993, name=None)
Energy = np.array([10000])
thetaB = Si111_130K.get_Bragg_angle(Energy) # - Si111_130K.get_dtheta_symmetric_Bragg(Energy)  
_,_,_,chi0,chih,chihbar= Si111_130K.get_F_chi(E=Energy,sinThetaOverLambda=np.sin(thetaB)/(physconsts.CH/Energy))
print('**SRW 晶体中的chi0, chih 和chih bar **')
print('chi0 = ', chi0)
print('chih = ', chih)
print('chihbar = ', chihbar)
print('**SRW d_space**')
print('Bragg angle = ', thetaB[0], ' rad',', and ',thetaB[0]*180/np.pi,'degree')
print('d_space = ',Si111_130K.d)