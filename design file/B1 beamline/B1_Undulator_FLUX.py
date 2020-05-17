# -*- coding: utf-8 -*-
r"""
Taper模式下，插入件光源的性能计算，该程序下给出：
    （1）gap-K-B-E1之间的关系；
    （2）不同K值下，一定口径范围theta\psi内，通量谱的分布曲线；
"""
#%%引入模块
import os, sys; 
import matplotlib
if os.name == 'nt':
    print('the OS is Windows, interative plot is used')
if os.name == 'posix':
    print('the OS is Linux, interative plot is not used')
    matplotlib.use('Agg')
import matplotlib.pyplot as plt    

import numpy as np
import pickle, imageio, time
import numpy as np
from scipy.interpolate import UnivariateSpline
import xrt.backends.raycing  as raycing
import xrt.backends.raycing.sources as rs

strExDataFolderName = 'output'
cwd = os.getcwd()

precisionOpenCL = 'auto'
K2B = 10.710201593926415

#%%光源参数定义
E0 = 100000
dE = E0*0.001/2
eMin0,eMax0 = E0-dE,E0+dE
kwargs_B1 = dict(
        eE = 6.0, eI = 0.2, eEspread = 0.00111,
        eEpsilonX = 0.02755, eEpsilonZ = 0.002755,  
        betaX = 2.83871, betaZ = 1.91667,
        eMin=eMin0, eMax=eMax0, eN=51, targetE = [100000,11],
#        distE='BW',   
        K = 'auto', period=17, n=251)#

kwargs_HEPS = kwargs_B1
source = rs.Undulator(name='B1', **kwargs_HEPS, precisionOpenCL = precisionOpenCL)
print('K value = ', source.K)

energy = np.linspace(-1,1,3)*2e2 + E0
theta = np.linspace(-1, 1, 31) * 20e-6  
psi = np.linspace(-1, 1, 31) * 20e-6  
dtheta, dpsi = theta[1] - theta[0], psi[1] - psi[0]

# %%
#选择几个能量点，画出二维分布曲线
energy_plot = np.linspace(-1,1,3)*1e1 + E0 #np.array([10000,11100])
theta = np.linspace(-1, 1, 31) * 20e-6  
psi = np.linspace(-1, 1, 31) * 20e-6  
dtheta, dpsi = theta[1] - theta[0], psi[1] - psi[0]
#计算数据
I_data = []
for energy_i in energy_plot:
    I0, l1, l2, l3 = source.intensities_on_mesh(np.array([energy_i]), theta, psi) 
#    print(source.get_SIGMAP(energy_plot))
    intensity = I0[0,:,:]
    I_data.append(intensity)
#画图    
I_data = np.array(I_data)
vmin = np.amin(I_data)
vmax = np.amax(I_data)
#%%
gif_images = [] 
for ie, intensity in enumerate(I_data):
    plt.figure(ie + 10)
    im =plt.imshow(intensity.T, interpolation='bicubic', cmap='jet',
                extent=[min(theta)*1e6,max(theta)*1e6,min(psi)*1e6,max(psi)*1e6], aspect='auto')#,vmin=vmin, vmax=vmax)   
    plt.xlabel('x(urad)')
    plt.ylabel('z(urad)')
    plt.xlim([min(theta)*1e6,max(theta)*1e6])
    plt.ylim([min(psi)*1e6,max(psi)*1e6])
    # 计算sigma值和FWHM值
    Ix = np.sum(intensity,axis =1)#intensity[:, len(psi)//2]
    Iy = np.sum(intensity,axis =0)#intensity[len(theta)// 2, :]
    varIx = ((Ix * theta**2).sum() / Ix.sum())**0.5
    varIy = ((Iy * psi**2).sum() / Iy.sum())**0.5
    spline = UnivariateSpline(theta, Ix-Ix.max()/2, s=0)
    r1, r2 = spline.roots()  # find the roots
    FWHMx = r2 - r1
    spline = UnivariateSpline(psi, Iy-Iy.max()/2, s=0)
    r1, r2 = spline.roots()  # find the roots
    FWHMy = r2 - r1    
    
    print('Sigma:', varIx, varIy)
    print('FWHM:',FWHMx, FWHMy)
    plt.colorbar(im, shrink=1)
    plt.title('energy = {0} eV'.format(energy_plot[ie])) 
    filename = os.path.join(os.getcwd(), strExDataFolderName, 'distribution_{0}.png'.format(ie))
    plt.savefig(filename, dpi = 300)  
    gif_images.append(imageio.imread(filename))
#保存数据
pickleName = os.path.join(os.getcwd(), strExDataFolderName, '2D_data.pickle')
with open(pickleName, 'wb') as f:
    pickle.dump([I_data,energy_plot, theta, psi], f, protocol=2)

#保存动态图
filename = os.path.join(os.getcwd(), strExDataFolderName, "Dynamic_image.gif")
imageio.mimsave(filename,gif_images,fps=1.5)

# %%


I_XRT = I_data[1,:,:]
x_xrt = theta*30

I_SRW = np.loadtxt('100keV.dta', skiprows=2, usecols=[0, 1, 2], unpack=True)
I_SRW = np.reshape(data_SPECTRA[2,:],(51,51))
fig1 = plt.figure(1, figsize=(7, 5))
ax = plt.subplot(111, label='1')
ax.set_xlabel(u'mm')
ax.set_ylabel(u'flux density(a.u.) ')

ax.plot(theta*30*1e3, I_XRT[:,len(psi)//2]/max(I_XRT[:,len(psi)//2]), label='XRT x' , lw=2)
ax.plot(psi*30*1e3, I_XRT[len(theta)// 2, :]/max(I_XRT[len(theta)// 2, :]), label='XRT y' , lw=2)
ax.plot(np.linspace(-0.6,0.6,51), I_SRW[:,25]/max(I_SRW[:,25]), label='SPECTRA x' , lw=2)
ax.plot(np.linspace(-0.6,0.6,51), I_SRW[25, :]/max(I_SRW[25, :]), label='SPECTRA y' , lw=2)



ax.set_ylim([0,1.1])
ax.legend(loc='upper right')
plt.grid()
fig1.savefig('compareTaper1.png',dpi = 300)
plt.show()
(theta, Ix-Ix.max()/2)