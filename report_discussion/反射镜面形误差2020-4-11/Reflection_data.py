# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 11:57:39 2020

@author: yangfg
"""

import numpy as np
import matplotlib.pyplot as plt

SR = np.linspace(0.8, 0.999, 100)
SR = np.array([0.81,0.85,0.9,0.95,0.97,0.98,0.99])
ratio_height = 1/np.sqrt(1-SR)*4*np.pi
fig = plt.figure(1)
ax = plt.semilogy(SR, ratio_height)
plt.xlabel('Strehl Ratio',fontsize = 14)
plt.ylabel('(λ/θ) / $σ_{M}$',fontsize = 14)
plt.grid(which='both',axis='both')
plt.show()
plt.savefig('SR_ratio.png', dpi = 300)



