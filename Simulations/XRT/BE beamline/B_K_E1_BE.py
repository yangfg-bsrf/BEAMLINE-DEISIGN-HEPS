# -*- coding: utf-8 -*-
r"""
插入件调谐计算
"""
import os, sys; 
import matplotlib.pyplot as plt
import numpy as np
import xrt.backends.raycing  as raycing
import xrt.backends.raycing.sources as rs
precisionOpenCL = 'auto'
K2B = 10.710201593926415


#%%光源参数定义
E0 = 8000
dE = E0*0.0001/2
eMin0,eMax0 = E0-dE,E0+dE
kwargs_B8 = dict(
        eE = 6.0, eI = 0.2, eEspread = 0.00111,
        eEpsilonX = 0.02755, eEpsilonZ = 0.002755,  
        betaX = 10.12, betaZ = 9.64,
        eMin=eMin0, eMax=eMax0, eN=51,
        distE='BW',   
        K = 0.7833,period=32.7,n=151)
kwargs_HEPS = kwargs_B8
source = rs.Undulator(name='BE', **kwargs_HEPS, precisionOpenCL = precisionOpenCL)

energy = np.linspace(5000, 15000, 500)
theta = np.linspace(-1, 1, 11) * 10e-6  
psi = np.linspace(-1, 1, 11) * 10e-6  
dtheta, dpsi = theta[1] - theta[0], psi[1] - psi[0]


"""
调谐计算：基波能量，K值，gap，和磁场
"""
Kmax, gmin = 2.44264, 11  # needed to calculate gap(K)
Ks = np.linspace(0.5, Kmax, 100)

gaps = np.log(Kmax/Ks) / np.pi * source.L0 + gmin
B = K2B * Ks / source.L0

E1 = 9.5*source.Ee**2/(1+Ks**2/2)/source.L0*1e3


plt.figure(1)
plt.plot(E1/1000, Ks)
ax = plt.gca()
ax.set_xlabel(u'1st energy (keV)',fontsize=13)
ax.set_ylabel(u'K',fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
xlims, ylims = (E1[0],E1[-1]), (0, 1)
plt.grid(which='both')
plt.show()
plt.savefig('characteristics_tuning_Ks.png', dpi = 300)

plt.figure(2)
plt.plot(E1/1000, gaps)
ax = plt.gca()
ax.set_xlabel(u'1st energy (keV)',fontsize=13)
ax.set_ylabel(u'gaps (mm)',fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
xlims, ylims = (E1[0],E1[-1]), (0, 1)
plt.grid(which='both')

plt.savefig('characteristics_tuning_gaps.png', dpi = 300)

plt.figure(3)
plt.plot(E1/1000, B)
ax = plt.gca()
ax.set_xlabel(u'1st energy (keV)',fontsize=13)
ax.set_ylabel(u'magnet field (T)',fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
xlims, ylims = (E1[0],E1[-1]), (0, 1)
plt.grid(which='both')
plt.show()
plt.savefig('characteristics_tuning_B.png', dpi = 300)



E1 = []
for iK, K in enumerate(Ks):
    source.Ky = K
    source.reset()  
    E1.append(source.E1)
    

#harmonics = [1, 3, 5, 7, 9, 11]
#colors = ['b', 'g', 'r', 'c', 'm', 'y']
#
#

#tunesE, tunesF = [], []
#
#tmpKy = source.Ky
#for iK, K in enumerate(Ks):
#    if raycing._VERBOSITY_ > 10:
#        print("K={0}, {1} of {2}".format(K, iK+1, len(Ks)))
#    source.Ky = K
#    source.reset()  
#    I0 = source.intensities_on_mesh(energy, theta, psi, harmonics)[0]
#    flux = I0.sum(axis=(1, 2)) * dtheta * dpsi
#    argm = np.argmax(flux, axis=0)
#    fluxm = np.max(flux, axis=0)
#    tunesE.append(energy[argm] / 1000.)
#    tunesF.append(fluxm)
#source.Ky = tmpKy
#source.reset()
#
#
## plot:
#for tuneE, tuneF, harmonic, color in zip(tunesE, tunesF, harmonics, colors):
#    plt.loglog(tuneE, tuneF, 'o-', label='{0}'.format(harmonic),
#               color=color)
## format the graph:
#ax = plt.gca()
#ax.set_xlabel(u'energy (eV)')
#ax.set_ylabel(u'flux through (20 µrad)² (ph/s/0.1% bw)')
##ax.set_xlim(xlims)
##ax.set_ylim(ylims)
##ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
#plt.grid(which='both')
#ax.legend(loc='upper right', title='harmonics')
#
## labels of K and gap:
#for tuneE, tuneF, harmonic in zip(tunesE, tunesF, harmonics):
#    if harmonic == 3:
#        ax.text(tuneE[-1], tuneF[-1]*1.05, 'K=', fontsize=8,
#                ha='right', va='bottom')
#        ax.text(tuneE[-1]/1.05, tuneF[-1]/1.2, 'gap=\n(mm)', fontsize=8,
#                color='b', ha='right', va='center')
#    for x, y, K, gap in zip(tuneE, tuneF, Ks, gaps):
#        if xlims[0] < x < xlims[1] and ylims[0] < y < ylims[1]:
#            ax.text(x, y*1.05, '{0:.2f}'.format(K), fontsize=8)
#            if harmonic == 3:
#                ax.text(x/1.05, y/1.2, '{0:.2f}'.format(gap), fontsize=8,
#                        color='b')
#
#plt.savefig('calc_undulator_tune.png')
#plt.show()