#辐射追迹画图
import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits.axisartist.axislines import SubplotZero
plt.rcParams['font.sans-serif']=['SimHei'] # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False # 用来正常显示负号
#原则性参数
mr = 50
residual_r = 30


def pointu(pointa, pointb, x): #两点连线，确定y坐标
    y = (x-pointa[0])/(pointb[0]-pointa[0])*(pointb[1]-pointa[1])+pointa[1]
    return y

def line(pointa, pointb, x_axis, apertures):  #确定可以通过口径的光线  
    y_axis = pointu(pointa, pointb, x)
    x_max = [max(x)]
    for aperture in apertures:
        y = pointu(pointa, pointb, aperture[2])
        if (y < aperture[0])|(y > aperture[1]):    #在口径外的将不画出
            x_max.append(aperture[2])
    
    x_end = min(x_max) #最大横坐标位置
    for i0, aperture in enumerate(apertures):  #定位挡光位置
        if x_end == max(x):
            indexi = len(apertures)+1
        if x_end == aperture[2]:
            indexi = i0
    x_axis_new = x_axis[x_axis <= x_end]
    y_axis_new = y_axis[x_axis <= x_end]  
    return x_axis_new, y_axis_new, indexi

def plot_Radiation_2p(point_rc, point_rbc, x, apertures): #画出所有的线
    x_new, y_new, indexi = line(point_rc, point_rbc, x, apertures)
    if max(x_new) == max(x):
        plt.plot(x_new, y_new, 'r-', linewidth=linwidth1) 
    else:
        plt.plot(x_new, y_new, 'k-' ,linewidth=linwidth2) 

def plot_Radiation(sr, C, x, apertures):
    plot_Radiation_2p(sr.point_rc, C.point_rbc, x, apertures)
    plot_Radiation_2p(sr.point_rc, C.point_lbc, x, apertures)
    plot_Radiation_2p(sr.point_lc, C.point_rbc, x, apertures)
    plot_Radiation_2p(sr.point_lc, C.point_lbc, x, apertures)


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
        if (indexi < len(aper)) & (indexi > 0) :
            yi_max = max(yi[:,-1])
            xii = xi[yi[:,-1] == yi_max][0]
            yii = yi[yi[:,-1] == yi_max][0]
            plt.plot(xii,yii, 'k-' ,linewidth = linwidth1)
            if yi_max > aper[indexi][1]:
                plt.plot([xii[-1], xii[-1]],[yii[-1], yii[-1] + mr],color = 'red',linewidth = 2)
                Cs[indexi].wailunkuo = yii[-1] + mr - Cs[indexi].center[1]
    
            yi_min = min(yi[:,-1])
            xii = xi[yi[:,-1] == yi_min][0]
            yii = yi[yi[:,-1] == yi_min][0]
            plt.plot(xii,yii, 'k-' ,linewidth = linwidth1)
            if yi_min < aper[indexi][0]:  
                plt.plot([xii[-1], xii[-1]],[yii[-1], yii[-1] - mr],color = 'red',linewidth = 2)
                Cs[indexi].wailunkuo = -yii[-1] + mr + Cs[indexi].center[1]
                

def line_range(pointa, pointb, x_axis, apertures):  #横坐标x_axis，纵坐标y_axis   
    y_axis = pointu(pointa, pointb, x)
    x_min = [0]
    for aperture in apertures:
        y = pointu(pointa, pointb, aperture[2])
        if (y < aperture[0])|(y > aperture[1]):    #在口径外的将不画出
            x_min.append(aperture[2])     #记录位置  
    x_axis_new = x_axis[(x_axis>=max(x_min))& (x_axis<=pointa[0])]#
    y_axis_new = y_axis[(x_axis>=max(x_min)) & (x_axis<=pointa[0])]# 
    
    return x_axis_new, y_axis_new

def plot_Radiation_2p_range(point_a, point_b, x, apertures):
    x_new, y_new = line_range(point_a, point_b, x, apertures)
    if min(x_new) == min(x):
        plt.plot(x_new, y_new, 'blue', linewidth=linwidth1) 
    else:#光线被截止
        print('no plot')
        # plt.plot(x_new, y_new, 'b' ,linewidth=linwidth2) 

        
def plot_Radiation_range(C1, C2, x, apertures):
    plot_Radiation_2p_range(C1.point_rfc, C2.point_lbc, x, apertures)
    plot_Radiation_2p_range(C1.point_rfc, C2.point_rbc, x, apertures)
    #plot_Radiation_2p_range(C1.point_lfc, C2.point_lbc, x, apertures)
    #plot_Radiation_2p_range(C1.point_lfc, C2.point_rbc, x, apertures)

def plot_Radiation_range2(Ce, Cs, x, apertures):
    x0 = []
    y0 =[]
    for Ci in Cs:
        x_new1, y_new1= line_range(Ce.point_rfc, Ci.point_lbc, x, apertures)
        x_new2, y_new2 = line_range(Ce.point_rfc, Ci.point_rbc, x, apertures)
        if min(x_new1) == min(x):
            x0.append(x_new1)
            y0.append(y_new1)
        if min(x_new2) == min(x):
            x0.append(x_new2)
            y0.append(y_new2)   
             
    y0 = np.array(y0)
    x0 = np.array(x0)
    ymin = min(y0[:,0])

    for i, yi in enumerate( y0):
        if yi[0] == ymin:
             x_new = x0[i,:]
             y_new = yi
    plt.plot(x_new, y_new, 'blue', linewidth=linwidth1)  
        


#%%光源
class RSource(object):
    #中心位置p，尺寸q  
    def __init__(self, center, aperture, name):
        self.center = center
        self.width = aperture
        self.name = name
        self.generate_point()
    def generate_point(self):
        #水平点，前后
        self.point_rc = [self.center[0], self.center[1] + self.width/2]
        self.point_lc = [self.center[0], self.center[1] - self.width/2]        
   
#%%准直器
class Collimator(object):
    #中心位置p，尺寸q，y轴是光束传播方向，xz是截面方向    
    def __init__(self, center, aperture, thickness, name):
        self.center = center
        self.width = aperture
        self.thickness = thickness
        self.name = name
        self.wailunkuo = self.width + 100
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


#%%辐射防护器件定义
if True:
    fs_components = []         
    ur = RSource([0, 0],40e-3, 'Undulator for BD') #插入件，位置和尺寸口径
    sr = RSource([3e3, 0], 15, 'Source for BR')# 真空盒
    C1 = Collimator([20.3e3 + sr.center[0], 0], 20, 300, 'C1')#位置、口径、厚度
    C2 = Collimator([26.35e3 + sr.center[0], 0], 20, 300, 'C2')#位置、口径、厚度
    C3 = Collimator([27.04e3 + sr.center[0], 0], 30, 200, 'C3')#位置、口径、厚度
    C4 = Collimator([27.8e3 + sr.center[0], 0], 20, 700, 'Ratched Wall& C4')#位置、口径、厚度
    C5 = Collimator([30.84e3 + sr.center[0], 0], 30, 300, 'C5')#位置、口径、厚度  
    C6 = Collimator([43e3, 46], 5, 300, 'C6')#位置、口径、厚度
    fs_components = [C1, C2, C4, C6] #准直器件组
    fs_components_all = [C1,C2,C3,C4,C5,C6]
    
#确定每个y轴位置的边界
if True:  
    xticks = [] 
    apertures = [] 
    for component in fs_components:
        apertures.append(component.apertureb)
    for component in fs_components_all:
        xticks.append(component.center[0] + component.thickness/2)

#%%光学器件定义
if True:
    centers = []
    #白光狭缝
    angle_accep = 25e-6
    p_WB = 35e3
    aper_WB = Slit([p_WB, 0], p_WB*angle_accep ,'White Beam Slit')

    #白光准直镜
    theta_WM = 2e-3
    p_WBM = 36.5e3
    length = angle_accep*p_WBM/theta_WM
    WBM = Mirror([p_WBM, 0], length, theta_WM, 'white beam mirror')

    #单色器
    p_DCM = 42.55e3
    
    theta_DCM = np.deg2rad(45)
    fixedOffset = 20
    DCM_gap = fixedOffset*np.sin(theta_DCM)/np.sin(2*(theta_DCM))
    # DCM_gap = 5
    # fixedOffset = DCM_gap/np.sin(theta_DCM)*np.sin(2*(theta_DCM))

    p_DCM01 = p_DCM
    DCM01y = (p_DCM01-p_WBM)*np.tan(2*theta_WM) 
    DCM02y = DCM01y + DCM_gap/np.sin(theta_DCM)*np.sin(2*(theta_DCM+theta_WM))
    DCM02x = DCM_gap/np.sin(theta_DCM)*np.cos(2*(theta_DCM + theta_WM)) 
    p_DCM02 = p_DCM01 + DCM02x
    DCM01 = Mirror([p_DCM01, DCM01y], 30, theta_DCM + 2*theta_WM, 'DCM01')
    DCM02 = Mirror([p_DCM02, DCM02y], 30, theta_DCM + 2*theta_WM, 'DCM02')
    
    #谐波镜HRM
    p_HRM = 45.9e3
    theta_HRM = 2e-3
    fixedOffset = 1
    HRM_gap = fixedOffset*np.sin(theta_HRM)/np.sin(2*(theta_HRM))

    p_HRM01 = p_HRM
    HRM01y = DCM02y + (p_HRM01-p_DCM02)*np.tan(2*theta_WM)  
    HRM02y = HRM01y + HRM_gap/np.sin(theta_HRM)*np.sin(2*(theta_HRM+theta_WM))
    HRM02x = HRM_gap/np.sin(theta_HRM)*np.cos(2*(theta_HRM + theta_WM)) 
    p_HRM02 = p_HRM01 + HRM02x
    HRM01 = Mirror([p_HRM01, HRM01y], 300, theta_HRM + 2*theta_WM, 'HRM01')
    HRM02 = Mirror([p_HRM02, HRM02y], 300, theta_HRM + 2*theta_WM, 'HRM02')    
    
    #VFM镜
    p_VFM = 48.9e3
    theta_VFM = 1.7e-3
    VFMy = HRM02y + (p_VFM-p_HRM02)*np.tan(2*theta_WM) 
    VFM = Mirror([p_VFM, VFMy], 300, theta_HRM, 'HRM02')
    
    #二次光源
    p_secslit= p_VFM + 4.5e3
    secslit = Slit([p_secslit, VFMy], 1 ,'White Beam Slit')
    
    centers.append([0, 0])
    centers.append(aper_WB.center)    
    centers.append(WBM.center)
    centers.append(DCM01.center)
    centers.append(DCM02.center)
    centers.append(HRM01.center)
    centers.append(HRM02.center)  
    centers.append(VFM.center) 
    centers.append(secslit.center)
    
    #画图参数
    xticks.append(aper_WB.center[0])
    xticks.append(WBM.center[0])
    xticks.append(DCM01.center[0])
    xticks.append(HRM01.center[0])
    xticks.append(VFM.center[0])    
    centers = np.array(centers)

#%%坐标轴设计
if True:
#    fig, ax = plt.subplots(figsize = (16,8))
    fig = plt.figure(1, (16, 8))
    ax = SubplotZero(fig, 1, 1, 1)
    fig.add_subplot(ax)    

    fig.set_tight_layout('tight')
    x = np.linspace(sr.center[0], 50e3, 1000)
    x_limit = [0,max(x)]
    y_limit = [-200,200]
    plt.xlim(x_limit)
    plt.ylim(y_limit)
    
    #隐藏右边和上边的边框，使之没有颜色,隐藏默认坐标轴（上下左右边框），并新建坐标轴X-Y，同时设置刻度标识方向
    ax.axis[:].set_visible(False)
    ax.axis["bottom","top"].set_visible(True)
    ax.axis["x"] = ax.new_floating_axis(0, 0)
    ax.axis["y"] = ax.new_floating_axis(1, sr.center[0])
    #new_floating_axis(self, nth_coord, value, axis_direction='bottom')
    #
    ax.axis["x"].set_axis_direction('top')
    ax.axis["y"].set_axis_direction('left')
    
    ax.axis["x"].set_axisline_style("->", size = 2.0)
    ax.axis["y"].set_axisline_style("->", size = 2.0)

    new_axisline = ax.get_grid_helper().new_fixed_axis
    ax.axis["new2"] = new_axisline(loc="right", offset=(0, 0), axes=ax)
    
    new_axisline = ax.get_grid_helper().new_fixed_axis
    ax.axis["new3"] = new_axisline(loc="left", offset=(0, 0), axes=ax)     
    #画图参数
    linwidth1 = 1
    linwidth2 = 0.1
    
#%%#画器件
if True:
    #%%准直器件
    plt.text(ur.point_lc[0]+600,y_limit[1]*0.5,ur.name,fontdict={'size':'12','color':'b','rotation':'90'})
    plt.gca().add_patch(plt.Rectangle((sr.point_lc[0],sr.point_lc[1]), 100, sr.width,  color = 'red'))    
    #光源
    plt.text(sr.point_lc[0]+600,y_limit[1]*0.5,sr.name,fontdict={'size':'12','color':'b','rotation':'90'})
    plt.gca().add_patch(plt.Rectangle((sr.point_lc[0],sr.point_lc[1]), 100, sr.width,  color = 'red'))
    
    #准直器类
    L = 100
    for C in [C1, C2, C3, C4, C5, C6]:
        if C.name == 'Ratched Wall& C4':
            color_block = 'black'
        elif C.name == 'C3':
            color_block = 'yellow'
        elif C.name == 'C5':
            color_block = 'yellow'
        else:
            color_block = 'blue'          
        plt.gca().add_patch(plt.Rectangle((C.point_rbc[0], C.point_rbc[1]), C.thickness, L, color = color_block))
        plt.gca().add_patch(plt.Rectangle((C.point_lbc[0], C.point_lbc[1] - L), C.thickness, L, color = color_block))
        plt.text(C.point_rbc[0]-400,y_limit[1]*0.7,"{}".format(C.name),fontdict={'size':'12','color':'b','rotation':'90'})            

    #%%光学器件
    #白光狭缝
    plt.gca().add_patch(plt.Rectangle((aper_WB.point_rbc[0], aper_WB.point_rbc[1]), 300, L, color='gray'))
    plt.gca().add_patch(plt.Rectangle((aper_WB.point_lbc[0], aper_WB.point_lbc[1] - L), 300, L, color='gray'))   
    plt.text(aper_WB.point_lbc[0]-600,y_limit[1]*0.3,"White Slit",fontdict={'size':'12','color':'g','rotation':'90'}) 
    #白光镜
    plt.plot([WBM.point_lbc[0],WBM.point_rbc[0]], [WBM.point_lbc[1],WBM.point_rbc[1]],linewidth =2,color = 'gray')
    plt.text(WBM.point_lbc[0]-00,y_limit[1]*0.2,"WBM",fontdict={'size':'12','color':'g','rotation':'90'}) 
    #单色器
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
    #
    #辐射追迹option1：画出所有线
#    for C in fs_components:
#        plot_Radiation(sr, C, x, apertures[0:4])
    #辐射追迹option2：画出部分光线
    plot_Radiation2(sr, fs_components, x, apertures)
    
    #从墙器件C4出发，投射光线  
    plot_Radiation_range2(C4, [C1, C2], x,apertures[0:2])      
        
    # for C in [C1,C2,C3]:
    #     plot_Radiation_range(C4, C, x, apertures[0:3])

#%%画光轴,主要辅助光线
if True:
    plt.plot(centers[:,0], centers[:,1], 'k-', linewidth = 2)

    labelx = []
    for item in xticks:
        labelx.append('{0:.1f}m'.format(item*1e-3))         
    plt.xticks(xticks,labelx,rotation='vertical')
    plt.grid(axis='x')
    ax.axis["x"].major_ticklabels.set_rotation('vertical')
    plt.annotate(s='', xy=(VFM.center[0],VFM.center[1]), xytext=(VFM.center[0]+4e3,VFM.center[1]), arrowprops=dict(arrowstyle='<-',linewidth =2))
    #plt.arrow(VFM.center[0], VFM.center[1], 4e3, 0,length_includes_head=True, head_width=0.2, lw=2) #画箭头方式2

    plt.axhline(y = VFM.center[1]-25, color = 'k', linestyle = '--', linewidth = 2)
    plt.axvline(x = ur.center[0], color = 'k', linestyle = '--', linewidth = 2)

#%%补充辅助信息
if True:
    text_positionx = 5e3
    text_positiony = y_limit[1]*0.6
    text_indicator = 'Collimator Info:'
    for C in fs_components_all:
        if C.name == 'C1':
            C.wailunkuo = 0
        if C.name == 'C3':
            C.wailunkuo = 0    
        if C.name == 'C5':
            C.wailunkuo = 0 
        C.update_info()
        text_indicator = text_indicator+ '\n' + C.info
    plt.text(text_positionx,text_positiony, text_indicator,color = 'black',  fontsize = 12,bbox={'facecolor':'green', 'alpha':0.5, 'pad':10})

plt.show()
plt.savefig('Radiation_protection.png', dpi =300)

"""
#    ax.spines['right'].set_color('black')
#    ax.yaxis.set_label_position("right")
#    ax.spines['top'].set_color('none')
    #默认ax里面的x轴和y轴
#    ax.xaxis.set_ticks_position('bottom')
#    ax.yaxis.set_ticks_position('left')
    #移动x轴y轴的位置
#    ax.spines['bottom'].set_position(('data',0))#0,就是移到0的位置
#    ax.spines['bottom'].set_linewidth(1.5)
#    ax.spines['left'].set_position(('data',0))#0,就是移到0的位置
#    ax.spines['left'].set_linewidth(1.5)

#    #墙
#    C = C4
#    plt.gca().add_patch(plt.Rectangle((C.point_rbc[0], C.point_rbc[1]), C.thickness, L, color = 'black'))
#    plt.gca().add_patch(plt.Rectangle((C.point_lbc[0], C.point_lbc[1] - L), C.thickness, L, color = 'black'))
#    plt.text(C.point_rbc[0]-600,y_limit[1]*0.5,"{}".format(C.name),fontdict={'size':'12','color':'b','rotation':'90'})   


"""