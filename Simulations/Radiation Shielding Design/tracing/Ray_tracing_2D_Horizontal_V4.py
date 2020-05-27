#辐射追迹画图
import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.axisartist.axislines import SubplotZero
plt.rcParams['font.sans-serif']=['SimHei'] # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False # 用来正常显示负号

def pointu(pointa, pointb, x): #两点连线，确定y坐标
    y = (x-pointa[0])/(pointb[0]-pointa[0])*(pointb[1]-pointa[1])+pointa[1]
    return y

def line(pointa, pointb, x_axis, aper):  #确定可以通过口径的光线  
    y_axis = pointu(pointa, pointb, x)
    x_max = [max(x)]
    for aperture in aper:
        y = pointu(pointa, pointb, aperture[2])
        if (y < aperture[0])|(y > aperture[1]):    #在口径外的将不画出
            x_max.append(aperture[2])   
    x_end = min(x_max) #确定最大横坐标位置
    
    for i0, aperture in enumerate(apertures):  #定位挡光位置，为最后一个屏时，indexi为aperture总长度加1
        if x_end == max(x):    
            indexi = len(apertures) + 1
        if x_end == aperture[2]:
            indexi = i0
    x_axis_new = x_axis[x_axis <= x_end]
    y_axis_new = y_axis[x_axis <= x_end]  
    return x_axis_new, y_axis_new, indexi

def plot_Radiation_2p(point_rc, point_rbc, x, aper): #画出所有的线
    x_new, y_new, indexi = line(point_rc, point_rbc, x, aper)
    if max(x_new) == max(x):
        plt.plot(x_new, y_new, 'r-', linewidth=linwidth1) 
    else:
        plt.plot(x_new, y_new, 'k-' ,linewidth=linwidth2) 

def plot_Radiation(sr, C, x, aper):
    plot_Radiation_2p(sr.point_rc, C.point_rbc, x, aper)
    plot_Radiation_2p(sr.point_rc, C.point_lbc, x, aper)
    plot_Radiation_2p(sr.point_lc, C.point_rbc, x, aper)
    plot_Radiation_2p(sr.point_lc, C.point_lbc, x, aper)

def plot_Radiation2(sr, Cs, x, aper): #画出所有的线，近画出分布特征线
    index0 = []
    x0 = []
    y0 =[]
    for C in Cs:  #记录每条光线的截至位置，寻找最大y值光线
        x_new, y_new, indexi = line(sr.point_rc, C.point_rbc, x, aper)
        x0.append(x_new)
        y0.append(y_new)  
        index0.append(indexi)       
        x_new, y_new, indexi = line(sr.point_rc, C.point_lbc, x, aper)
        x0.append(x_new)
        y0.append(y_new) 
        index0.append(indexi) 
        x_new, y_new, indexi = line(sr.point_lc, C.point_rbc, x, aper)
        x0.append(x_new)
        y0.append(y_new) 
        index0.append(indexi)        
        x_new, y_new, indexi = line(sr.point_lc, C.point_lbc, x, aper)
        x0.append(x_new)
        y0.append(y_new)
        index0.append(indexi) 
    print(list(set(index0)))             
    #获得每个屏上的光线极限,步骤：1.把在同一平面截至的光线统计到一起，2.寻找边界光线
    for indexi in list(set(index0)): 
        #1.把在同一平面截至的光线统计到一起
        yi = []
        xi = []
        for x_new, y_new, indexii in zip(x0, y0, index0):
            if indexi == indexii:
                yi.append(y_new)
                xi.append(x_new)
        yi = np.array(yi)
        xi = np.array(xi)
        
        #寻找极限边界光线，在狭缝的上下界
        if (indexi <= len(aper)) & (indexi > 0) :
            mr = Cs[indexi].mr
            yi_max = max(yi[:,-1])
    
            xii = xi[yi[:,-1] == yi_max][0]
            yii = yi[yi[:,-1] == yi_max][0]
            plt.plot(xii,yii, 'k-' ,linewidth = linwidth1)
            if yi_max > aper[indexi][1]:
                plt.plot([xii[-1], xii[-1]],[yii[-1], yii[-1] + mr],color = 'red',linewidth = 2)
                Cs[indexi].wailunkuo = (yii[-1] + mr - Cs[indexi].center[1])*2
    
            yi_min = min(yi[:,-1])
            xii = xi[yi[:,-1] == yi_min][0]
            yii = yi[yi[:,-1] == yi_min][0]
            plt.plot(xii,yii, 'k-' ,linewidth = linwidth1)
            if yi_min < aper[indexi][0]:  
                plt.plot([xii[-1], xii[-1]],[yii[-1], yii[-1] - mr],color = 'red',linewidth = 2)
                a0 = (-yii[-1] + mr + Cs[indexi].center[1])*2
                Cs[indexi].wailunkuo = a0
                # if Cs[indexi].wailunkuo < a0:
                #     Cs[indexi].wailunkuo = a0
                    
        if indexi == (len(aper)+1) : #Hutch最后墙
            yi_max = max(yi[:,-1])  
            xii = xi[yi[:,-1] == yi_max][0]
            yii = yi[yi[:,-1] == yi_max][0]
            plt.plot(xii,yii, 'k-' ,linewidth = linwidth1)       #追迹线 
            plt.plot([xii[-1], xii[-1]],[yii[-1], yii[-1] + mr],color = 'red',linewidth = 2) #辐射线
 
            yi_min = min(yi[:,-1])
            xii = xi[yi[:,-1] == yi_min][0]
            yii = yi[yi[:,-1] == yi_min][0]
            plt.plot(xii,yii, 'k-' ,linewidth = linwidth1)       #追迹线 
            plt.plot([xii[-1], xii[-1]],[yii[-1], yii[-1] - mr],color = 'red',linewidth = 2) #辐射线
   
        if indexi == (len(aper)):
            angle = np.rad2deg(np.arctan(yi[:,1]-yi[:,0])/(xi[:,1]-xi[:,0]))
        angle = [0.1]    
    return angle
                       
def line_range(pointa, pointb, x_axis, aper):  #横坐标x_axis，纵坐标y_axis   
    y_axis = pointu(pointa, pointb, x)
    x_min = [0]
    for aperture in aper:
        y = pointu(pointa, pointb, aperture[2])
        if (y < aperture[0])|(y > aperture[1]):    #在口径外的将不画出
            x_min.append(aperture[2])     #记录位置  
    x_axis_new = x_axis[(x_axis>=max(x_min))]#
    y_axis_new = y_axis[(x_axis>=max(x_min))]# 
    
    return x_axis_new, y_axis_new

def plot_Radiation_2p_range(point_a, point_b, x, aper):
    x_new, y_new = line_range(point_a, point_b, x, aper)
    if min(x_new) == min(x):
        plt.plot(x_new[x_new<=point_a[0]], y_new[x_new<=point_a[0]], 'blue', linewidth=linwidth1) 
    else:#光线被截止
        # print('no plot')
        plt.plot(x_new[x_new<=point_a[0]], y_new[x_new<=point_a[0]], 'orangered' ,linewidth=linwidth1) 
    
def plot_Radiation_range(C1, C2, x, aper):
    plot_Radiation_2p_range(C1.point_rfc, C2.point_lbc, x, aper)
    plot_Radiation_2p_range(C1.point_rfc, C2.point_rbc, x, aper)

def plot_Radiation_range_wall(Ce, x, aper):
    plot_Radiation_2p_range(Ce.point_rfc, Ce.point_lbc, x, aper)
    plot_Radiation_2p_range(Ce.point_lfc, Ce.point_rbc, x, aper)

def plot_Radiation_range2(Ce, Cs, x, C_end, aper):
#C_end最远端的准直器件截至光线    
    x0 = []
    y0 =[]
    for i0, Ci in enumerate(Cs):
        if i0 > 0:
            plot_Radiation_2p_range(Ce.point_rfc, Ci.point_lbc, x, aper)
            plot_Radiation_2p_range(Ce.point_lfc, Ci.point_rbc, x, aper) 
            
        x_new1, y_new1= line_range(Ce.point_rfc, Ci.point_lbc, x, aper)
        x_new2, y_new2 = line_range(Ce.point_rfc, Ci.point_rbc, x, aper)
        if min(x_new1) == min(x):
            x0.append(x_new1)
            y0.append(y_new1)
        if min(x_new2) == min(x):
            x0.append(x_new2)
            y0.append(y_new2)   
             
    y0 = np.array(y0)
    x0 = np.array(x0)
    ymin = min(y0[:,0])

    for i, yi in enumerate(y0): #寻找最下边的光线并画图
        if yi[0] == ymin:
             x_new = x0[i,:]
             y_new = yi
    plt.plot(x_new[x_new<C_end.center[0]], y_new[x_new < C_end.center[0]], 'royalblue', linewidth=linwidth1) 
    angle = np.rad2deg(np.arctan(y_new[1]-y_new[0])/(x_new[1]-x_new[0]))
    return angle
    
#*************************************************************************#
#################################器件定义###################################
#*************************************************************************#
#%%光源
class RSource(object):
    #中心位置p，尺寸q  
    def __init__(self, center, aperture, name):
        self.center = center
        self.width = aperture
        self.name = name
        self.info = self.name + ': aperure = {0:.2f}mm'.format(self.width)       
        self.generate_point()
    def generate_point(self):
        #水平点，前后
        self.point_rc = [self.center[0], self.center[1] + self.width/2]
        self.point_lc = [self.center[0], self.center[1] - self.width/2]        
   
#%%准直器
class Collimator(object):
    #中心位置p，尺寸q，y轴是光束传播方向，xz是截面方向    
    def __init__(self, center, aperture, L, thickness, mr, name):
        self.center = center
        self.width = aperture
        self.thickness = thickness
        self.name = name
        self.mr = mr
        self.wailunkuo = L
        self.info = self.name + ': inner = {0:.2f}mm'.format(self.width) + ', out = {0:.2f}mm'.format(self.wailunkuo)
        self.generate_point()
        self.generate_boundary()
    def generate_point(self):
        #水平点，前后左右x，y
        self.point_rbc = [self.center[0], self.center[1] + self.width/2]
        self.point_lbc = [self.center[0], self.center[1] - self.width/2]
        
        self.point_rfc = [self.center[0] + self.thickness, self.center[1] + self.width/2]        
        self.point_lfc = [self.center[0] + self.thickness, self.center[1] - self.width/2]  
        
    def generate_boundary(self):
        #aperture的 y坐标范围，和x位置
        self.apertureb = [self.center[1] - self.width/2, self.center[1] + self.width/2, self.center[0]]      
        self.aperturef = [self.center[1] - self.width/2, self.center[1] + self.width/2, self.center[0]+self.thickness]   
    def update_info(self):
        self.info = self.name + ': inner = {0:.2f}mm'.format(self.width) + ', out = {0:.2f}mm'.format(self.wailunkuo)      

#%%狭缝
class Slit(object):
    #中心位置p，尺寸q，y轴是光束传播方向，xz是截面方向    
    def __init__(self, center, aperture, name):
        self.center = center
        self.width = aperture
        self.name = name
        self.generate_point()
        self.generate_boundary()
    def generate_point(self):
        #水平点，前后左右
        self.point_rbc = [self.center[0], self.center[1] + self.width/2]  
        self.point_lbc = [self.center[0], self.center[1] - self.width/2]
    def generate_boundary(self):
        #aperture的位置
        self.apertureb = [self.center[1] - self.width/2, self.center[1] + self.width/2, self.center[0]]     


#%%反射镜
class Mirror(object):
    #反射镜    
    def __init__(self, center, length, pitch, name):
        self.center = center
        self.pitch = pitch
        self.length = length
        self.width = self.length*np.sin(self.pitch)
        self.name = name
        self.generate_point()
        
    def generate_point(self):
        #水平点，前后左右
        self.point_rbc = [self.center[0] + self.length*np.cos(self.pitch)/2, self.center[1] + self.width/2]  
        self.point_lbc = [self.center[0] - self.length*np.cos(self.pitch)/2, self.center[1] - self.width/2]

#%%系统参数的定义
#原则性参数
residual_r = 30
mr1 = 30
mr2 = 50
#%%辐射防护器件定义
if True:
    sr = RSource([3e3, 0], 22, 'Source for BR')# 真空盒
    C1 = Collimator([20.3e3 + sr.center[0], 0], 20, 300, 300, mr2, 'C1')#位置、内口径、外口经、厚度、3*Moliere_radius，name
    C2 = Collimator([26.35e3 + sr.center[0], 0], 20, 150, 300, mr2, 'C2')
    C3 = Collimator([27.04e3 + sr.center[0], 0], 30, 100, 200, mr2, 'SS & C3')
    C4 = Collimator([27.8e3 + sr.center[0], 0], 20, 1000,700, mr2, 'Ratched Wall & C4')
    C5 = Collimator([35.5e3+300, 0], 15, 95, 300, mr2, 'C5')  
    C6 = Collimator([45.2e3-150+300, 0], 20, 200, 300, mr2, 'C6')

    #准直器件组定义1：用于追迹    
    fs_components = [C1, C2, C4, C5, C6] # 这里C3不考虑挡光效应
    apertures = [] #确定每个口径的边界位置
    for component in fs_components:
        apertures.append(component.apertureb) 
        
    #准直器件组定义2：用于画图 
    fs_components_all = [C1,C2,C3,C4,C5,C6] 
    xticks = [] #确定每个器件的中心位置
    for component in fs_components_all:
        xticks.append(component.center[0] + component.thickness/2)    


#%%光学器件定义
if True:
    centers = []
    #插入件光源
    ur = RSource([0, 0], 40e-3, 'Undulator for BD') #插入件，位置和尺寸口径
    #白光狭缝
    angle_accep = 25e-6
    p_WB = 35e3
    aper_WB = Slit([p_WB, 0], p_WB*angle_accep ,'White Beam Slit')

    #白光准直镜
    theta_WM = 1.7e-3
    p_WBM = 36.5e3+300
    length = angle_accep*p_WBM/theta_WM
    WBM = Mirror([p_WBM, 0], length, theta_WM, 'white beam mirror')

    #单色器
    theta_CCM = np.deg2rad(25)
    # CCM
    p_CCM = 40.5e3+300
    CCM_gap = 8
    fixedOffset = CCM_gap/np.sin(theta_CCM)*np.sin(2*(theta_CCM))

    p_CCM01 = p_CCM
    CCM01y = (p_CCM01-p_WBM)*np.tan(2*theta_WM) 
    CCM02y = CCM01y + CCM_gap/np.sin(theta_CCM)*np.sin(2*(theta_CCM+theta_WM))
    CCM02x = CCM_gap/np.sin(theta_CCM)*np.cos(2*(theta_CCM + theta_WM)) 
    p_CCM02 = p_CCM01 + CCM02x
    CCM01 = Mirror([p_CCM01, 0], 1, 0, 'CCM01')
    CCM02 = Mirror([p_CCM02, 0], 1, 0, 'CCM02')
    
    #谐波镜HRM
    p_HRM = 45.9e3+300
    theta_HRM = 1.7e-3
    # fixedOffset = 6
    # HRM_gap = fixedOffset*np.sin(theta_HRM)/np.sin(2*(theta_HRM))
    HRM_gap = 3
    fixedOffset = HRM_gap/np.sin(theta_HRM)*np.sin(2*(theta_HRM))
    
    p_HRM01 = p_HRM
    HRM01y = CCM02y + (p_HRM01-p_CCM02)*np.tan(2*theta_WM)  
    HRM02y = HRM01y + HRM_gap/np.sin(theta_HRM)*np.sin(2*(theta_HRM+theta_WM))
    HRM02x = HRM_gap/np.sin(theta_HRM)*np.cos(2*(theta_HRM + theta_WM)) 
    p_HRM02 = p_HRM01 + HRM02x
    HRM01 = Mirror([p_HRM01, 0], 300, 0, 'HRM01')
    HRM02 = Mirror([p_HRM02, 0], 300, 0, 'HRM02')    
    
    #VFM镜
    p_VFM = 48.9e3 +300
    theta_VFM = 1.7e-3
    VFMy = HRM02y + (p_VFM-p_HRM02)*np.tan(2*theta_WM) 
    VFM = Mirror([p_VFM, 0], 300, theta_VFM, 'VFM')
    
    #二次光源
    p_secslit= p_VFM + 4.5e3
    secslit = Slit([p_secslit, 0], 1 ,'White Beam Slit')
    
    centers.append([0, 0])
    centers.append(aper_WB.center)    
    centers.append(WBM.center)
    centers.append(CCM01.center)
    centers.append(CCM02.center)
    centers.append(HRM01.center)
    centers.append(HRM02.center)  
    centers.append(VFM.center) 
    centers.append(secslit.center)
    
    xticks.append(aper_WB.center[0])
    xticks.append(WBM.center[0])
    xticks.append(CCM01.center[0])
    xticks.append(HRM01.center[0])
    xticks.append(HRM02.center[0])
    xticks.append(VFM.center[0])  
    
    
    # DCM
    theta_DCM = theta_CCM
    p_DCM = 42.55e3+300
    fixedOffset = 20.5
    DCM_gap = fixedOffset*np.sin(theta_DCM)/np.sin(2*(theta_DCM))
    
    p_DCM01 = p_DCM
    DCM01y = (p_DCM01-p_WBM)*np.tan(2*theta_WM) 
    DCM02y = DCM01y + DCM_gap/np.sin(theta_DCM)*np.sin(2*(theta_DCM+theta_WM))
    DCM02x = DCM_gap/np.sin(theta_DCM)*np.cos(2*(theta_DCM + theta_WM)) 
    p_DCM02 = p_DCM01 + DCM02x
    DCM01 = Mirror([p_DCM01, 0], 30, theta_DCM + 2*theta_WM, 'DCM01')
    DCM02 = Mirror([p_DCM02, 0], 30, theta_DCM + 2*theta_WM, 'DCM02')
    #VFM镜
    p_VFM = 48.9e3+300
    theta_VFM = 1.7e-3
    VFMy = DCM02y + (p_VFM-p_DCM02)*np.tan(2*theta_WM) 
    VFM = Mirror([p_VFM, 0], 300, theta_VFM, 'VFM')    
    
    centers_DCM = []
    centers_DCM.append([0, 0])
    centers_DCM.append(aper_WB.center)    
    centers_DCM.append(WBM.center)
    centers_DCM.append(DCM01.center)
    centers_DCM.append(DCM02.center)
    centers_DCM.append(VFM.center) 
    centers_DCM.append(secslit.center)
    
    xticks.append(DCM01.center[0])

    #画图参数
    centers = np.array(centers)
    centers_DCM = np.array(centers_DCM)

#*************************************************************************#
#################################系统画图###################################
#*************************************************************************#

#%%坐标轴设计
if True:

    fig = plt.figure(1, (16, 8))
    ax = SubplotZero(fig, 1, 1, 1)
    fig.add_subplot(ax)    

    fig.set_tight_layout('tight')
    x = np.linspace(sr.center[0], 50e3, 1000)
    x_limit = [0, max(x)]
    y_limit = [-150, 150]
    plt.xlim(x_limit)
    plt.ylim(y_limit)
    
    #隐藏右边和上边的边框，使之没有颜色,隐藏默认坐标轴（上下左右边框），并新建坐标轴X-Y，同时设置刻度标识方向
    ax.axis[:].set_visible(False)
    ax.axis["bottom","top"].set_visible(True)
    ax.axis["bottom","top"].major_ticklabels.set_rotation('vertical')
    ax.axis["bottom","top"].major_ticklabels.set_pad(10)
    ax.axis["x"] = ax.new_floating_axis(0, 0)
    ax.axis["y"] = ax.new_floating_axis(1, sr.center[0])

    ax.axis["x"].set_axisline_style("->", size = 2.0)
    ax.axis["y"].set_axisline_style("->", size = 2.0)

    new_axisline = ax.get_grid_helper().new_fixed_axis
    ax.axis["new2"] = new_axisline(loc="right", offset=(0, 0), axes=ax)
    
    new_axisline = ax.get_grid_helper().new_fixed_axis
    ax.axis["new3"] = new_axisline(loc="left", offset=(0, 0), axes=ax)     
    #画图参数
    linwidth1 = 1
    linwidth2 = 0.1
    
 

#%%#画束线光学器件
if True:
    #插入件光源
    plt.text(ur.point_lc[0]+600,y_limit[1]*0.5,ur.name,fontdict={'size':'12','color':'b','rotation':'90'})
    plt.gca().add_patch(plt.Rectangle((sr.point_lc[0],sr.point_lc[1]), 100, sr.width,  color = 'red'))   

    #白光狭缝
    L = 20 #显示的外轮廓尺寸
    plt.gca().add_patch(plt.Rectangle((aper_WB.point_rbc[0]-150, aper_WB.point_rbc[1]), 300, L, color='gray'))
    plt.gca().add_patch(plt.Rectangle((aper_WB.point_lbc[0]-150, aper_WB.point_lbc[1] - L), 300, L, color='gray'))   
    plt.text(aper_WB.point_lbc[0]-600,y_limit[1]*0.3,"White Slit",fontdict={'size':'12','color':'g','rotation':'90'}) 
    #白光镜
    plt.plot([WBM.point_lbc[0],WBM.point_rbc[0]], [WBM.point_lbc[1],WBM.point_rbc[1]],linewidth =2,color = 'gray')
    plt.text(WBM.point_lbc[0]-00,y_limit[1]*0.2,"WBM",fontdict={'size':'12','color':'g','rotation':'90'}) 
    #单色器CCM
    plt.plot([CCM01.point_lbc[0],CCM01.point_rbc[0]], [CCM01.point_lbc[1],CCM01.point_rbc[1]],linewidth =2,color = 'g')
    plt.plot([CCM02.point_lbc[0],CCM02.point_rbc[0]], [CCM02.point_lbc[1],CCM02.point_rbc[1]],linewidth =2,color = 'g')
    plt.text(CCM02.point_lbc[0]-600,y_limit[1]*0.6,"CCM",fontdict={'size':'12','color':'g','rotation':'90'}) 
    
    #单色器DCM
    plt.plot([DCM01.point_lbc[0],DCM01.point_rbc[0]], [DCM01.point_lbc[1],DCM01.point_rbc[1]],linewidth =2,color = 'g')
    plt.plot([DCM02.point_lbc[0],DCM02.point_rbc[0]], [DCM02.point_lbc[1],DCM02.point_rbc[1]],linewidth =2,color = 'g')
    plt.text(DCM02.point_lbc[0]-600,y_limit[1]*0.6,"DCM",fontdict={'size':'12','color':'g','rotation':'90'}) 
    #谐波镜
    plt.plot([HRM01.point_lbc[0],HRM01.point_rbc[0]], [HRM01.point_lbc[1],HRM01.point_rbc[1]],linewidth =2,color = 'g')
    plt.plot([HRM02.point_lbc[0],HRM02.point_rbc[0]], [HRM02.point_lbc[1],HRM02.point_rbc[1]],linewidth =2,color = 'g')
    plt.text(HRM02.point_lbc[0]-600,y_limit[1]*0.8,"HRMS",fontdict={'size':'12','color':'g','rotation':'90'})   
    #VFM
    plt.plot([VFM.point_lbc[0],VFM.point_rbc[0]], [VFM.point_lbc[1],VFM.point_rbc[1]],linewidth =2,color = 'g')
    plt.text(VFM.point_lbc[0]-600,y_limit[1]*0.8,"VFM",fontdict={'size':'12','color':'g','rotation':'90'}) 

#%%准直器画线 
if True:   
    #辐射追迹option2：画出部分光线
    angle1 = plot_Radiation2(sr, fs_components, x, apertures)

    #从墙器件C4出发，投射光线  
    angle2 = plot_Radiation_range2(C4, [C1, C2], x, C5, apertures[0:2]) 

    #墙C4内外侧发出光线在准直器SS & C3上的分布
    plot_Radiation_range_wall(C4, x, [C3.apertureb])

    #视场角
    field_angle = max([max(angle1), angle2]) - min(angle1)    
    
#%%画光轴,主要辅助光线
if True:
    #光轴
    plt.plot(centers[:,0], centers[:,1], '-.', color = 'lawngreen',linewidth = 1)
    plt.plot(centers_DCM[:,0], centers_DCM[:,1], '-.',color = 'darkorchid', linewidth = 1)

    #横坐标ticks位置
    labelx = []
    xticks.append(max(x))
    for item in xticks:
        labelx.append('{0:.1f}m'.format(item*1e-3)) 
   
    plt.xticks(xticks,labelx,rotation='vertical')
    plt.grid(axis='x',linewidth = 0.2)
    ax.axis["x"].major_ticklabels.set_visible(False)

#%%补充辅助信息——TEXT文档
#%%补充辅助信息——TEXT文档
if True:
    text_positionx, text_positiony= 7e3, y_limit[1]*0.3
    Note_text = 'Vertical red line: 3*Moliere_radius' + '\n' + 'Field_angle = {0:.2f} Degree'.format(field_angle)
    text_indicator = 'Radiation Shielding Design Info.:'+ '\n' + 'Collimators: C1, C2, C5, C6' + '\n' + '\n' + sr.info 
    
    for C in [C1, C2, C3, C4, C5, C6]:
        C.update_info()
        if C.name == C1.name:
            C.info =   C.name + ': inner = {0:.2f}mm'.format(C.width) + ', out = --'
        # if C.name == C2.name:
        #     C.info =   C.name + ': inner = {0:.2f}mm'.format(C.width) + ', out = 150mm'            
        if C.name == C3.name:
            C.info =   C.name + ': inner = {0:.2f}mm'.format(C.width) + ', out = --'  
        if C.name == C4.name:
            C.info =   C.name + ': inner = {0:.2f}mm'.format(C.width) + ', out = ∞'          
        # if C.name == C5.name:
        #     C.info =   C.name + ': inner = {0:.2f}mm'.format(C.width) + ', out = 135mm'
        # if C.name == C6.name:
        #     C.info =   C.name + ': inner = {0:.2f}mm'.format(C.width) + ', out = 230mm'   
        text_indicator = text_indicator + '\n' + C.info

    text_indicator = text_indicator+ '\n' + '\n' + Note_text
    plt.text(text_positionx,text_positiony, text_indicator,color = 'black',  fontsize = 12,bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})

#%%#最后画准直器件
if True:
    #光源
    plt.text(sr.point_lc[0]+600,y_limit[1]*0.5,sr.name,fontdict={'size':'12','color':'b','rotation':'90'})
    plt.gca().add_patch(plt.Rectangle((sr.point_lc[0],sr.point_lc[1]), 100, sr.width,  color = 'darkblue'))
    
    #准直器类(前端固定)
    for C in [C1, C2, C3, C4, C5, C6]:
        if C.name == 'C4 & Ratched Wall':
            color_block = 'black'
            C.wailunkuo = 1000
        elif C.name == 'C3 & SS':
            color_block = 'yellow'
        else:
            color_block = 'darkblue'          
        plt.gca().add_patch(plt.Rectangle((C.point_rbc[0], C.point_rbc[1]), C.thickness, (C.wailunkuo-C.width)/2, color = color_block))
        plt.gca().add_patch(plt.Rectangle((C.point_lbc[0], C.point_lbc[1] - (C.wailunkuo-C.width)/2), C.thickness, (C.wailunkuo-C.width)/2, color = color_block))
        plt.text(C.point_rbc[0]-400,y_limit[1]*0.6,"{}".format(C.name),fontdict={'size':'12','color':'b','rotation':'90'})   

plt.show()
plt.savefig('Radiation_protection_Horizontal.png', dpi =300)





