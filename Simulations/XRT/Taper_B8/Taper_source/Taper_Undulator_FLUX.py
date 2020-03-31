# -*- coding: utf-8 -*-
r"""
Taper模式下，插入件光源的性能计算，该程序下给出：
    （1）gap-K-B-E1之间的关系；
    （2）不同K值下，一定口径范围theta\psi内，通量谱的分布曲线；

"""
#%%引入模块
import os, sys; 
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np
import pickle, imageio, time

import xrt.backends.raycing  as raycing
import xrt.backends.raycing.sources as rs

strExDataFolderName = 'output'
cwd = os.getcwd()

precisionOpenCL = 'auto'
K2B = 10.710201593926415

#%%光源参数定义
E0 = 10000
dE = E0*0.0001/2
eMin0,eMax0 = E0-dE,E0+dE
kwargs_B8 = dict(
        eE = 6.0, eI = 0.2, eEspread = 0.00111,
        eEpsilonX = 0.02755, eEpsilonZ = 0.002755,  
        betaX = 10.12, betaZ = 9.64,
        eMin=eMin0, eMax=eMax0, eN=51, #targetE = [10000,3],
        distE='BW',   
        K = 1.9648, period=35, n=142,taper=(1.09, 15.240))#

kwargs_HEPS = kwargs_B8
source = rs.Undulator(name='B8', **kwargs_HEPS, precisionOpenCL = precisionOpenCL)
print('K value = ', source.K)

energy = np.linspace(9600, 10400, 5)
theta = np.linspace(-1, 1, 11) * 9e-6  
psi = np.linspace(-1, 1, 11) * 9e-6  
dtheta, dpsi = theta[1] - theta[0], psi[1] - psi[0]

#%%调谐参数计算
"""
调谐计算：基波能量，K值，gap，和磁场
"""
Kmax, gmin = 2.8759, 11  # needed to calculate gap(K)
#给定K值，确定E1
if False:
    Ks = np.linspace(-1, 1, 1)*0e-2 + 1.9648
    gaps = np.log(Kmax/Ks) / np.pi * source.L0 + gmin
    B = K2B * Ks / source.L0
    E1 = 9.5*source.Ee**2/(1+Ks**2/2)/source.L0*1e3
    print(E1)

#给定E1值，确定K值
if True:
    E1 = np.linspace(-1,1,1)*0e1 + 3.333e3
    Ks =np.sqrt((9.5*source.Ee**2/source.L0/E1*1e3 - 1)*2)
    B = K2B * Ks / source.L0
    gaps = np.log(Kmax/Ks) / np.pi * source.L0 + gmin
    print('gaps = ',gaps)


#%%计算不同K下的通量-能量分布
fluxs = []
fluxs.append(energy)
plt.figure(1)
ax = plt.gca()
for iK, K in enumerate(Ks):
    source.Ky = K
    source.reset()
    flux = []
    for ie, energy_i in enumerate(energy):
        I0, l1, l2, l3 = source.intensities_on_mesh(np.array([energy_i]), theta, psi)
        flux.append(I0.sum(axis=(1, 2)) * dtheta * dpsi)
        print(ie/len(energy))
    fluxs.append(flux)   
    plt.plot(energy, flux, label = 'energy = {0:.2f} keV'.format(energy[iK]/1e3))
plt.xlabel(u'energy (ev)', fontsize =14)    
plt.ylabel(u'flux (ph/s/0.1% bw)', fontsize =14)
plt.legend()
plt.grid() 
# plt.show()  
strTrajOutFileName = "Spectroscopy_FLUX.png"
saveName =  os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName)
plt.savefig(saveName, dpi = 300) 

# %%
#保存数据
filename = os.path.join(os.getcwd(), strExDataFolderName, "XRT_data.txt")
fluxs0 = np.array(fluxs)
Ks0 = np.insert(Ks, 0, [0])
zt = np.vstack((Ks0, fluxs0.T))    
np.savetxt(filename, zt, fmt='%.3e', header='the first line is K value, \
the rows are energy and flux @ given K value')


# %%
#选择几个能量点，画出二维分布曲线
source.Ky = Ks[0]
source.reset()
energy_plot = np.linspace(-1,1,21)*3e2 + 10000 #np.array([10000,11100])
theta = np.linspace(-1, 1, 11) * 9e-6  
psi = np.linspace(-1, 1, 11) * 9e-6  
dtheta, dpsi = theta[1] - theta[0], psi[1] - psi[0]

gif_images = [] 
for ie, energy_i in enumerate(energy_plot):
    I0, l1, l2, l3 = source.intensities_on_mesh(np.array([energy_i]), theta, psi)
    intensity = I0[0,:,:]
    plt.figure(ie+100)
    im =plt.imshow(intensity.T, interpolation='bicubic', cmap='jet',
                extent=[min(theta)*1e6,max(theta)*1e6,min(psi)*1e6,max(psi)*1e6], aspect='auto')#,vmin=vmin, vmax=vmax)
    
    plt.xlabel('x(urad)')
    plt.ylabel('z(urad)')
    plt.xlim([min(theta)*1e6,max(theta)*1e6])
    plt.ylim([min(psi)*1e6,max(psi)*1e6])
    plt.colorbar(im, shrink=1)
    plt.title('energy = {0} eV'.format(energy_i)) 
    filename = os.path.join(os.getcwd(), strExDataFolderName, 'distribution_{0}.png'.format(ie))
    plt.savefig(filename, dpi = 300)  
    gif_images.append(imageio.imread(filename))


#保存动态图
filename = os.path.join(os.getcwd(), strExDataFolderName, "Dynamic_image.gif")
imageio.mimsave(filename,gif_images,fps=1.5)

  
# %%
