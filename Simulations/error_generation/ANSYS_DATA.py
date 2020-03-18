# -*- coding: utf-8 -*- #
__author__ = "Fugui Yang"
__date__ = "18 March 2020"
"""
本程序用于从有限元分析计算产生面形误差数据
并保存为固定格式，同时画图显示面形情况；单独运行
"""
import os, sys; 
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

def read_ANSYS_bench(input_fname):
    f = open(input_fname)
    datas0 = np.loadtxt(input_fname,skiprows=1,delimiter='	')
    f.close()
    datas1 = np.unique(datas0,axis=0)  
    thredhold = 1e-10
    datas1[abs(datas1[:,2])<thredhold,2] = datas1[abs(datas1[:,2])<thredhold,2]*0
    datas1[abs(datas1[:,1])<thredhold,1] = datas1[abs(datas1[:,1])<thredhold,1]*0
    
    ind = np.lexsort((datas1[:, 1],datas1[:, 2]))
    datas2 = datas1[ind]
    y_01 = np.unique(datas2[:,1])          
    for i, t in enumerate(y_01): 
        if i%2 == 1:
            datas2 = datas2[datas2[:,1] != t]            
    y0,x0 = np.unique(datas2[:,1]),np.unique(datas2[:,2]) 
    z0 = datas2[:, 4]
    nX,nY = len(x0),len(y0)                 
    z = z0.reshape((nX,nY)) 
    x,y,z = x0,y0,z
    return x,y,z

#保存面形数据至标准格式
input_fname='B5-DCM-7kev.txt'    ########
name_1='DCM'+ '_'                ########

xi,yi,zi=read_ANSYS_bench(input_fname)
zi0=zi.T
np.savetxt(name_1+'X.txt', xi,fmt='%.8e')
np.savetxt(name_1+'Y.txt', yi,fmt='%.8e')
np.savetxt(name_1+'Z.txt', zi0,fmt='%.8e')
yi00 = np.insert(yi, 0, 0)*1e-3
zi00 = np.insert(zi0*1e-3, 0, xi*1e-3,axis=0)
zi01 = np.insert(zi00, 0, yi00,axis=1)
np.savetxt(name_1+'Z1.dat', zi01,fmt='%.7e', delimiter='\t', newline='\n')
   
####显示数据        
warpX,warpY,warpZ=xi,yi,zi    
warpNX, warpNY = len(warpX), len(warpY)
warpA, warpB = np.gradient(warpZ)
dx = warpX[1] - warpX[0]
dy = warpY[1] - warpY[0]
warpA = np.arctan(warpA/dx)
warpB = np.arctan(warpB/dy)
rmsA = ((warpA**2).sum() / (warpNX*warpNY))**0.5
rmsB = ((warpB**2).sum() / (warpNX*warpNY))**0.5    

warpSplineA = ndimage.spline_filter(warpA)
warpSplineB = ndimage.spline_filter(warpB) 
warpSplineZ = ndimage.spline_filter(warpZ) 

nsamples = 400
xmin, xmax =np.min(warpX),np.max(warpX)
ymin, ymax =np.min(warpY),np.max(warpY)   

###设定分析区域尺寸
Lm = ymax*1.5               ########
Ls = 0.6                    ########

x0 = np.linspace(-Ls/2, Ls/2, nsamples)
y0 = np.linspace(-Lm/2, Lm/2, nsamples)


"""
子午方向
"""
coords = np.array([(x0*0-xmin) /(xmax-xmin) * (warpNX-1),\
                   (y0-ymin)/(ymax-ymin) * (warpNY-1)])
a = ndimage.map_coordinates(warpSplineA, coords, prefilter=True)   #该因子仅仅是为了显示方便
b = ndimage.map_coordinates(warpSplineB, coords, prefilter=True)
z0 = ndimage.map_coordinates(warpSplineZ, coords, prefilter=True) 
z1 = ndimage.map_coordinates(zi, coords, prefilter=True) 
nx, ny=int(warpNX/2), int(warpNY/2) 
  
#面形误差整体作图分析--子午方向 
plt.figure(200,[16, 8]) 

plt.subplot(221)      #Height profile
plt.plot(warpY, warpZ[nx,:]*1.e6) 
plt.plot(y0,z1*1.e6)
plt.xlabel(u'Y (mm)',fontsize=14)
plt.ylabel(u'height (nm)',fontsize=14)
PV0=np.max(z1)-np.min(z1)
plt.title("Tangential profile = %.3f μrad rms and %.3f nm rms %.3f nm PV \n  total rms = %.3f μrad  and %.3f nm "
          %(np.std(b)*1e6,np.std(z1)*1e6, PV0*1e6,np.std(warpB[nx,:])*1e6,np.std(z0)*1e6))  

plt.subplot(223)      #Slope Profile
plt.plot(warpY,warpB[nx,:]*1.e6)
plt.plot(y0,b*1.e6)
plt.xlabel(u'Y (mm)',fontsize=14)
plt.ylabel(u'slope (μrad)',fontsize=14)
 
 
#选择区域的分析计算--子午方向    
plt.figure(211) 
plt.subplot(211)      #Height profile 
p = np.polyfit(y0,z1,2)
z11 = np.polyval(p,y0)
z_1h=z1-z11
PV0=np.max(z1)-np.min(z1)
PV1=np.max(z_1h)-np.min(z_1h)
radius = 1/p[0]/2*1e-3      
plt.plot(y0,(z1-np.mean(z1))*1.e6)
plt.plot(y0,z_1h*1.e6,label = 'residual error = %.3f nm rms %.3f nm PV'%(np.std(z_1h)*1e6, PV1*1e6))
plt.ylabel(u'height (nm)',fontsize=14)
plt.legend()
plt.title("Tangential profile fit radius = %.1f m " %(radius) ) 

plt.subplot(212)      #Slope Profile  
p = np.polyfit(y0,b,1)
z11 = np.polyval(p,y0)
z_1s=b-z11
plt.plot(y0,b*1.e6)   
plt.plot(y0,z_1s*1.e6, label = 'residual error = %.3f μrad rms'%(np.std(z_1s)*1e6))   

plt.xlabel(u'Y (mm)',fontsize=14)
plt.ylabel(u'slope (μrad)',fontsize=14)
plt.legend()

plt.show() 
plt.savefig('TANG_slope.png',dpi=300)  



"""
弧矢方向
"""
coords = np.array([(x0-xmin) /(xmax-xmin) * (warpNX-1),\
                   (y0*0-ymin)/(ymax-ymin) * (warpNY-1)])
a = ndimage.map_coordinates(warpSplineA, coords, prefilter=True)   #该因子仅仅是为了显示方便
b = ndimage.map_coordinates(warpSplineB, coords, prefilter=True)
z0 = ndimage.map_coordinates(warpSplineZ, coords, prefilter=True) 
z1 = ndimage.map_coordinates(zi, coords, prefilter=True) 

#单色器面形误差整体作图分析--弧矢方向 
plt.figure(200) 
plt.subplot(222)      #Height profile
plt.plot(warpX, warpZ[:,ny]*1.e6) 
plt.plot(x0,z1*1.e6)
plt.xlabel(u'X (mm)',fontsize=14)
plt.ylabel(u'height (nm)',fontsize=14)
PV0=np.max(z1)-np.min(z1)
plt.title("Sagittal profile = %.3f μrad rms and %.3f nm rms %.3f nm PV \n  total rms = %.3f μrad  and %.3f nm "
          %(np.std(a)*1e6,np.std(z1)*1e6, PV0*1e6,np.std(warpA[:,ny])*1e6,np.std(z0)*1e6))  

plt.subplot(224)     #Slope Profile
plt.plot(warpX,warpA[:,ny]*1.e6)
plt.plot(x0,a*1.e6)
plt.xlabel(u'X (mm)',fontsize=14)
plt.ylabel(u'slope (μrad)',fontsize=14)
plt.show() 
plt.savefig('total selection.png')


#单色器面形误差画图弧矢方向    
plt.figure(212) 

plt.subplot(211)     #Height profile  
p = np.polyfit(x0,z1,2)
z11 = np.polyval(p,x0)
z_1h=z1-z11
PV0=np.max(z1)-np.min(z1)
PV1=np.max(z_1h)-np.min(z_1h)
radius = 1/p[0]/2*1e-3   
 
plt.plot(x0,(z1-np.mean(z1))*1.e6)
plt.plot(x0,z_1h*1.e6, label = 'residual error = %.3f nm rms %.3f nm PV'%(np.std(z_1h)*1e6, PV1*1e6))
plt.ylabel(u'height (nm)',fontsize=14)
plt.legend()
plt.title("Sagittal profile fit radius = %.1f m " %(radius) )     
     
plt.subplot(212)     #Slope profile  
p = np.polyfit(x0,a,1)
z11 = np.polyval(p,x0)
z_1s=a-z11
radius = 1/p[0]*1e-3

plt.plot(x0,a*1.e6) 
plt.plot(x0,z_1s*1.e6, label = 'residual error = %.3f μrad rms'%(np.std(z_1s)*1e6))   
plt.xlabel(u'X (mm)',fontsize=14)
plt.ylabel(u'slope (μrad)',fontsize=14)
plt.legend()
 
plt.show() 
plt.savefig('Sag_slope.png',dpi=300)  
