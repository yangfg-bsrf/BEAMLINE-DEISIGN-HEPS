# -*- coding: utf-8 -*-
__author__ = "Fugui Yang"
__date__ = "11 Sept 2019"
"""
Version 1:
本程序用于连接ESRF-DABAM数据库中的数据，
生成面形误差数据，并保存为固定格式（供XRT和SRW使用），同时画图显示面形情况；单独运行
"""
import os, sys; sys.path.append(os.path.join('C:\\Users\Lenovo\Desktop\XRT_Simulation\surface error\Error_Generation'))
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import dabam as dabamA

"""
方法：读取测量一维面形数据,操作包括读取元数据
"""
def read_DABAM(x, ii, L):      #x是弧矢方向的采样点位置，，ii是DABAM数据库中的面型序列
    dm = dabamA.dabam()    
    dm.set_input_entryNumber(ii)
    #dm.set_input_outputFileRoot("dabam-002") # write files by default
    dm._set_from_command_line()   # get arguments of dabam command line
    
    if dm.get_input_value("runTests"): # if runTests selected
        dm.set_input_outputFileRoot("")      # avoid output files
        dm.test_dabam_names()
        dm.test_dabam_stdev_slopes()
    elif dm.get_input_value("summary"):
        print(dm.dabam_summary())
    else:
        dm.load()        # access data    
        if dm.get_input_value("plot") != None:
            dm.plot()
    """
    plt.figure(100)
    #Height profile
    plt.subplot(211)
    plt.plot(1e3*dm.y,1e6*dm.zHeights*1e3) 
    plt.xlabel("Y [mm]")
    plt.ylabel("Z [nm]")
    plt.title("M%d profile rms=%.2f nm and %.2f urad, \n Shape: %s, Facility: %s\n"%(dm.inputs.get('entryNumber'),1e9*dm.stdev_profile_heights(),1e6*dm.stdev_profile_slopes(),dm.metadata['SURFACE_SHAPE'],dm.metadata['FACILITY']))
    #Slope Profile
    plt.subplot(212)
    plt.plot(1e3*dm.y,1e6*dm.zSlopes)
    plt.xlabel("Y [mm]")
    plt.ylabel("Zp [urad]")
    plt.savefig('DATA_%d.png'%(dm.inputs.get('entryNumber')))
    plt.show()   
    """
    z_t=1e3*dm.zHeights

    y=1e3*dm.y
    x=x-(np.min(x)+np.max(x))/2  
    y=y-(np.min(y)+np.max(y))/2
    #拟合去除球面项    
    y00=[]
    z00=[]
    for i,t in enumerate(y):
        if np.abs(t)<L/2:
            y00.append(t)
            z00.append(z_t[i])
              
    p = np.polyfit(y00,z00,2)
    z1 = np.polyval(p,y)
    z_t=z_t-z1
    """
    plt.figure()
    plt.plot(y,z_t)
    plt.plot(y00,z00)
    plt.plot(y,z1)
    plt.show()
    """
    z_t=z_t[::-1]
    
    z0=[]
    for i in np.arange(0,len(x),1):
        z0.append(z_t)
    z=np.asarray(z0) 
    return x, y, z

def get_distorted_surface(fname1,fname2,fname3):
    x = np.loadtxt(fname1, unpack=True)
    y = np.loadtxt(fname2, unpack=True)
    z = np.loadtxt(fname3, unpack=True)
    return x, y, z
    
MM = 45   #反射镜序号
L = 550  #目标镜长
xi= np.arange(-5,5,1)       
xi,yi,zi=read_DABAM(x=xi, ii=MM, L=L)  #x是弧矢方向的采样点位置，，ii是DABAM数据库中的面型序列
zi=zi

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####STEP-3：保存数据，（5设置文件名）    
####保存数据，x，y，z（mm）
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
name_1='VCM'+ str(MM)+ '_'
zi = zi*0.3/0.335
zi0=zi.T
np.savetxt(name_1+'X.txt', xi,fmt='%.8e')
np.savetxt(name_1+'Y.txt', yi,fmt='%.8e')
np.savetxt(name_1+'Z.txt', zi0,fmt='%.8e')
yi00 = np.insert(yi, 0, 0)*1e-3
zi00 = np.insert(zi0*1e-3, 0, xi*1e-3,axis=0)
zi01 = np.insert(zi00, 0, yi00,axis=1)
np.savetxt(name_1+'Z1.dat', zi01,fmt='%.7e', delimiter='\t', newline='\n')

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
####STEP-4：显示数据  

###设定分析区域尺寸
xmin, xmax =np.min(xi),np.max(xi)
ymin, ymax =np.min(yi),np.max(yi)   

Lm = L             
Ls = 0.6

      
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
plt.title(name_1+ "Tangential profile fit radius = %.1f m " %(radius) ) 

plt.subplot(212)      #Slope Profile  
p = np.polyfit(y0,b,1)
z11 = np.polyval(p,y0)
z_1s=b-z11
plt.plot(y0,b*1.e6)   
plt.plot(y0,z_1s*1.e6, label = 'residual error = %.3f μrad rms'%(np.std(b)*1e6))   

plt.xlabel(u'Y (mm)',fontsize=14)
plt.ylabel(u'slope (μrad)',fontsize=14)
plt.legend()

plt.show() 
plt.savefig(name_1+ 'TANG_slope.png',dpi=300)  



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
plt.savefig(name_1+'total selection.png')


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
plt.title(name_1+ "Sagittal profile fit radius = %.1f m " %(radius) )     
     
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
plt.savefig(name_1+'Sag_slope.png',dpi=300)  

  
 
