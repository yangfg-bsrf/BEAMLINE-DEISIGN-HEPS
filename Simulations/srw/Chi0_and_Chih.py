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
Energy = 10000
thetaB=Si111_130K.get_Bragg_angle(Energy)
_,_,_,chi0,chih,_= Si111_130K.get_F_chi(E=Energy,sinThetaOverLambda=np.sin(thetaB)/(physconsts.CH/Energy))
print('SRW 晶体中的chih Π 和chih σ 选用下面相同的计算值chih')
print('chi0 = ', chi0)
print('chih = ', chih)