# -*- coding: utf-8 -*-
"""
2020.2.28：
    几何追迹
2020.3.5：
    波动传播
2020.3.27：
    Taper模式；
    增添自动识别系统功能；
    指定文件路径输出
    修正KB镜参数设置
2020.4.23：
    修改为BE束线，修正
    
"""
__author__ = "Fugui Yang"
__date__ = "28 February 2020"
import_lib = 1
if import_lib == 1:
    import os, sys
    import numpy as np
    import pickle, imageio, time
    #交互式画图选项
    import matplotlib as mpl
    if os.name == 'posix':
        print('the OS is Linux, interative plot is not used')
        mpl.use('Agg')
    import matplotlib.pyplot as plt
    # xrt函数引用
    import xrt.backends.raycing as raycing
    import xrt.runner as xrtr
    import xrt.plotter as xrtp
    import xrt.backends.raycing.apertures as ra
    import xrt.backends.raycing.oes as roes
    import xrt.backends.raycing.screens as rsc
    import xrt.backends.raycing.waves as rw
    import Distorted_Element as Distorted

raycing.precisionOpenCL = 'float64'
raycing.targetOpenCL = "auto"
#指定输出路径
strExDataFolderName = 'output'
cwd = os.getcwd()
Error_path = './surface error_DATA/'
#从保存的图片中读取图片，并保存动态图片
def dynamic_images(Plots): 
    gif_images = []    
    for plot in Plots:
        filename = plot.saveName
        gif_images.append(imageio.imread(filename)) 
    strTrajOutFileName = "dynamic_images.gif"
    saveName_dynamic =  os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName)
    imageio.mimsave(saveName_dynamic,gif_images,fps=1.5) 
    print("dynamic images : Done!")

#离焦曲线画图
def defocus_plot(Plots, dqs): 
    qCurve = []
    xCurve = []
    zCurve = []
    czCurve = []
    cxCurve = []
    for dq, plot in zip(dqs, Plots):
        qCurve.append(dq)            
        xCurve.append(plot.dx)
        zCurve.append(plot.dy)
        cxCurve.append(plot.cx)
        czCurve.append(plot.cy)

    plt.figure(100)
    ax1 = plt.subplot(111)
    ax1.set_title(r'FWHM size of beam cross-section near focal position')
    ax1.plot(qCurve, xCurve, 'o', label='Horizontal direction')  
    ax1.plot(qCurve, zCurve, '+', label='Vertical direction')
    ax1.set_xlabel(u'Shift (mm)', fontsize=14)
    ax1.set_ylabel(u'FWHM size ($\mu$m)', fontsize=14)
    # ax1.title(u'focal spot size', fontsize=14)
    plt.legend()          
    plt.show()
    strTrajOutFileName = "Depth_Of_Focus.png"
    path =  os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName)
    plt.savefig(path, dpi = 200)

    plt.figure(300)
    ax1 = plt.subplot(111)
    ax1.plot(qCurve, cxCurve, 'o', label='Horizontal direction')  
    ax1.plot(qCurve, czCurve, '+', label='Vertical direction')
    ax1.set_xlabel(u'Shift (mm)', fontsize=14)
    ax1.set_ylabel(u'Position ($\mu$m)', fontsize=14)
    # ax1.title(u'Position shift', fontsize=14)
    plt.legend()          
    plt.show()
    strTrajOutFileName = "position.png"
    path =  os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName)
    plt.savefig(path, dpi = 200)

    print('plot_generator Done')

    #保存截面图数据
    dump = []
    I_x = []
    I_z = []
    for plot in Plots:
        x = plot.xaxis.binCenters
        y = plot.yaxis.binCenters
        Ixy = plot.total2D
        ny = len(y)
        nx = len(x)
        I_x.append(Ixy[ny//2,:])
        I_z.append(Ixy[:,nx//2])
    dump.append([x, y, dqs, I_x,I_z])
    strTrajOutFileName = "image_section.pickle"
    pickleName =  os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName)
    with open(pickleName, 'wb') as f:
        pickle.dump(dump, f, protocol=2)
    print(plot.title," Defocus Data save :Done")

#画图数据保存,强度分布及坐标轴信息
def save_data(plot):
    dump = []
    x = plot.xaxis.binCenters
    y = plot.yaxis.binCenters
    Ixy = plot.total2D
    dump.append([x, y, Ixy])
    fileName = '{0}'.format(plot.title)
    strTrajOutFileName = '{0}.pickle'.format(fileName)
    pickleName =  os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName)
    with open(pickleName, 'wb') as f:
        pickle.dump(dump, f, protocol=2)
    print(plot.title," Data save :Done")

#**********************************************************#
""" ********************全局角度******************** """
mode = 0
if mode == 0:
    accept_ang_x, accept_ang_z = 40e-6, 40e-6 #亚微米模式
else:
    accept_ang_x, accept_ang_z = 10e-6, 1e-6 #微米聚焦模式
E0, n0 = 8000, 1  # 插入间工作能量点，及谐波级次
nrays, repeats =5e4, 1
n1 = 1 # 待分析谐波级次
# what = 'rays'
what = 'hybrid'

#器件误差设置
oe_error_WM = 0
oe_error_DCM = 0
oe_error_TM = 0

""" ********************材料特性******************** """
if True:
    #单色器材料
    crystalSi01 = raycing.materials.CrystalSi(name=None,hkl= (1, 1, 1))
    crystalSi0n = raycing.materials.CrystalSi(name=None,hkl= (n1/n0, n1/n0, n1/n0))
    #反射镜材料
    Si = raycing.materials.Material(elements=r"Si",kind=r"mirror",rho=2.33,name=None)
    Pt = raycing.materials.Material('Pt', rho=21.45, kind='mirror')
    Rh = raycing.materials.Material(elements=r"Rh",kind=r"mirror",rho=12.41,name=None)

""" 光源 """
if True:
    E1 = E0/n0*n1  # 待分析的能量点
    if what == 'rays':
        uniformRayDensity = False
        filamentBeam = False
        isMono = False
        prefix = 'BE-rays-'
    else:
        uniformRayDensity = True
        filamentBeam = True
        isMono = True
        prefix = 'BE-hybr-'
   
    dE = E1*0.001/2
    eMin0, eMax0 = E1-dE, E1+dE
    kwargs_SR = dict(
        nrays = nrays, center = [0, 0, 0],
        eE = 6.0, eI = 0.2, 
        # eEspread = 0.00111,
        # eEpsilonX = 0.02728, eEpsilonZ = 0.00276, 
        eEspread = 0,
        eEpsilonX = 0.00, eEpsilonZ = 0.00, 
        betaX = 10.12, betaZ = 9.64,
        xPrimeMax = accept_ang_x, zPrimeMax = accept_ang_z,
        xPrimeMaxAutoReduce = False, zPrimeMaxAutoReduce = False,
        targetE = [E0+123, n0], eMin = eMin0, eMax = eMax0,
        uniformRayDensity = uniformRayDensity,
        filamentBeam = filamentBeam, period=32.7,n=152,
        #K =  0.7576907940983183, taper(1.5,22.82838746))
        K = r"auto")
   

    beamLine = raycing.BeamLine(alignE=E0)
    beamLine.source = raycing.sources.Undulator(beamLine, **kwargs_SR)
    print('K = ',beamLine.source.K)

""" 第一光学元件：白光狭缝@35m """
if True:
    p_slitWB = 33000
    slitDx, slitDz = accept_ang_x*p_slitWB, accept_ang_z*p_slitWB
    opening = [-slitDx/2, slitDx/2, -slitDz/2, slitDz/2]
    beamLine.slit_WB = ra.RectangularAperture(
        beamLine, 'Slit_WB',  [0, p_slitWB, 0], 
        ('left', 'right', 'bottom', 'top'), opening = opening )

""" 第二光学元件：白光镜@36m """  
if True:
    p_WM = 36000
    theta_WM = 3e-3
    Surface_name='WM_'

    kwargs_WMirr = dict(
            name='WhiteMirror',center=[0, p_WM, 0],
            pitch = theta_WM,
            material = Rh,
            limPhysX=[-1, 1],limPhysY=[-300, 300]) 
    if oe_error_WM==0:
        #不带面形误差 
        beamLine.WMirr = roes.OE(bl=beamLine,**kwargs_WMirr)       
    else:       
        #带面形误差   
        kwargs_WMirr['get_distorted_surface'] ='error'
        kwargs_WMirr['fname1'] =Surface_name+'X.txt'
        kwargs_WMirr['fname2'] =Surface_name+'Y.txt'
        kwargs_WMirr['fname3'] =Surface_name+'Z.txt'
        beamLine.WMirr = Distorted.PlaneMirrorDistorted(bl=beamLine, **kwargs_WMirr) 

""" 第三光学元件：单色器@45m """
if True:
    p_DCM = 45000  
    Surface_name = 'DCM_'
    theta_DCM = crystalSi01.get_Bragg_angle(E1)-crystalSi01.get_dtheta_symmetric_Bragg(E1) 
    fixedOffset = 20
    DCM_gap = fixedOffset*np.sin(theta_DCM)/np.sin(2*(theta_DCM))
    # DCM_gap = 5
    # fixedOffset = DCM_gap/np.sin(theta_DCM)*np.sin(2*(theta_DCM))

    p_DCM01 = p_DCM
    DCM01x = 0
    DCM01z = (p_DCM01-p_WM)*np.tan(2*theta_WM)
    
    DCM02z = DCM01z + DCM_gap/np.sin(theta_DCM)*np.sin(2*(theta_DCM+theta_WM))
    DCM02y = DCM_gap/np.sin(theta_DCM)*np.cos(2*(theta_DCM + theta_WM)) 
    p_DCM02 = p_DCM01 + DCM02y
    DCM02x = 0 

    DCMDx,DCMDz = accept_ang_x*p_DCM/2*1.1, accept_ang_z*p_DCM/2/np.sin(theta_DCM)*1.1

    #使用双平晶代表，在tracing部分的代码也是不同的
    Surface_name= Error_path + 'DCM_'   
    kwargs_DCM01 = dict(
            name='DCM01',center=[DCM01x,p_DCM01, DCM01z],
            pitch = theta_DCM+theta_WM*2,
            material=crystalSi0n,
            limPhysX=[-DCMDx, DCMDx],limPhysY=[-DCMDz, DCMDz]) 
    kwargs_DCM02 = dict(
            name='DCM02',center=[DCM02x, p_DCM02, DCM02z],
            pitch = -theta_DCM-theta_WM*2, 
            positionRoll=-np.pi,
            material=crystalSi0n,
            limPhysX=[-DCMDx, DCMDx],limPhysY=[-DCMDz, DCMDz])    
    if oe_error_DCM == 0:
        #不带面形误差 
        beamLine.DCM01 = roes.OE(bl=beamLine,**kwargs_DCM01)      
    else:       
        #带面形误差   
        kwargs_DCM01['get_distorted_surface'] ='error'
        kwargs_DCM01['fname1'] =Surface_name+'X.txt'
        kwargs_DCM01['fname2'] =Surface_name+'Y.txt'
        kwargs_DCM01['fname3'] =Surface_name+'Z.txt'
        beamLine.DCM01 = Distorted.PlaneMirrorDistorted(bl=beamLine, **kwargs_DCM01)       
    beamLine.DCM02 = roes.OE(bl=beamLine,**kwargs_DCM02) 

""" （忽略）第四光学元件：单色器狭缝@46m,狭缝尺寸41.5*8.3 = 344um；37*5.5 = 228um"""
if True:
    p_slitDCM = p_DCM+1
    dx_slitDCM, dz_slitDCM = p_slitDCM*40e-6, p_slitDCM*40e-6
    print(dx_slitDCM, dz_slitDCM)
    beamLine.slit_DCM = raycing.apertures.RectangularAperture(
        beamLine, 'Slit_DCM',   [0, p_slitDCM,
                                DCM02z + (p_slitDCM-p_DCM02)*np.tan(theta_WM*2)], 
                                ('left', 'right', 'bottom', 'top'),
        [-dx_slitDCM/2, dx_slitDCM/2, -dz_slitDCM/2, dz_slitDCM/2])

""" 第五光学元件：超环面聚焦镜@48m, 反射镜长1m """
if True: 
    p_TM = 48000
    p_samp = p_TM + 8000
    theta_TM = 3e-3
    TM_z = DCM02z + (p_TM- p_DCM02)*np.tan(2*theta_WM)
   
    R_TM = 1/(1/p_TM+1/(p_samp -p_TM))*2/np.sin(theta_TM)
    r_TM = 1/(1/p_TM+1/(p_samp -p_TM))*2*np.sin(theta_TM)
    Surface_name= Error_path + 'TM_'

    kwargs_TM = dict(
               name = 'TM',center = [0, p_TM, TM_z],
               pitch = -theta_TM, positionRoll = -np.pi,
               material = Rh,
               limPhysX = [-1, 1], limPhysY = [-500, 500],
               R = R_TM, r = r_TM,
               targetOpenCL=r"auto")   
    if oe_error_TM == 0:        
        beamLine.TM = roes.ToroidMirror(bl=beamLine,**kwargs_TM)
    else: 
        kwargs_TM['get_distorted_surface'] ='error'
        kwargs_TM['fname1'] =Surface_name+'X.txt'
        kwargs_TM['fname2'] =Surface_name+'Y.txt'
        kwargs_TM['fname3'] =Surface_name+'Z.txt'
        beamLine.TM = Distorted.ToroidMirrorDistorted(bl=beamLine, **kwargs_TM)

""" 第六样品处狭缝 @56m""" 
if True:
    p_imag = p_samp
    xs_slit,  zs_slit = 0.05, 0.05  
    beamLine.slitimag = ra.RectangularAperture(
                    beamLine, 'Slit_image', [0, p_imag, TM_z],
                    ('left', 'right', 'bottom', 'top'),
                    [-xs_slit/2, xs_slit/2, -zs_slit/2, zs_slit/2])     

""" 观察屏组件 """
if True:
    zbins, zppb = 128, 2 
    xbins, xppb = 128, 2

    beamLine.fsm0 = rsc.Screen(beamLine, 'FSM0')
    beamLine.fsm1 = rsc.Screen(beamLine, 'FSM1')
    beamLine.fsm2 = rsc.Screen(beamLine, 'FSM2')
    beamLine.fsmsamp = rsc.Screen(beamLine, 'samp')
    beamLine.fsmSecSlit = rsc.Screen(beamLine, 'SecSlit')

    x_lim, z_lim = 20, 20
    edges = np.linspace(-x_lim, x_lim, xbins+1)*1e-3
    fsmExpX = (edges[:-1] + edges[1:]) * 0.5
    edges = np.linspace(-z_lim, z_lim, zbins+1)*1e-3
    fsmExpZ = (edges[:-1] + edges[1:]) * 0.5  

    dqs = np.linspace(-0, 0, 1)
    beamLine.fsmn = rsc.Screen(beamLine, name = 'FSMn') 
    
""""光线追迹模块"""
def run_process_rays(beamLine):
    print('***************************')
    print('run_process_rays  user')
    if isMono:
        fixedEnergy = E0
    else:
        fixedEnergy = False
    waveOnSlit = beamLine.slit_WB.prepare_wave(beamLine.source, nrays)
    beamSource = beamLine.source.shine(wave=waveOnSlit, fixedEnergy=fixedEnergy)  # 直接得到白光狭缝slit后的光场分布    
    beamLine.fsm0.center = beamLine.slit_WB.center
    beamFSM0 = beamLine.fsm0.expose(beam=beamSource)

    """ White collimating mirror """  
    WMGlobal, WMLocal = beamLine.WMirr.reflect(beam=beamSource)

    """ DCM """ 
    DCM_in = WMGlobal
    DCM01Global, DCM01Local = beamLine.DCM01.reflect(beam=DCM_in)
    DCM02Global, DCM02Local = beamLine.DCM02.reflect(beam=DCM01Global)
    # beamLine.slit_DCM.propagate(beam=DCM02Global)
    beamLine.fsm1.center = beamLine.slit_DCM.center
    beamFSM1 = beamLine.fsm1.expose(beam=DCM02Global)
    DCM_OUT = DCM02Global

    """  Toroidal mirror  """       
    TMGlobal, TMLocal = beamLine.TM.reflect(DCM_OUT)

    """  Sample position  """  
    #1. 样品附近多观察点屏
    beam_imags = []
    for i, dq in enumerate(dqs):
        lsamp = [0, p_samp+dq, TM_z]     
        beamLine.fsmn.center = lsamp
        beam_imag = beamLine.fsmn.expose(TMGlobal)
        beam_imags.append(beam_imag)

    #2. 样品位置处
    lsamp = [0, p_samp, TM_z]   
    beamLine.slitimag.center = lsamp
    beamLine.slitimag.propagate(TMGlobal)
    beamLine.fsmsamp.center = lsamp
    beamsamp = beamLine.fsmsamp.expose(TMGlobal)

    outDict = {
        'beamSource':beamSource,
        'beamFSM0':beamFSM0,
        'WMGlobal': WMGlobal,'WMLocal': WMLocal,
        'DCM':beamFSM1,
        'DCM01Global': DCM01Global,'DCM01Local': DCM01Local, 'DCM02Global': DCM02Global,'DCM02Local': DCM02Local,
        'TMGlobal': TMGlobal, 'TMLocal': TMLocal,
        }  
    outDict['Samp'] = beamsamp
    for i, dq in enumerate(dqs):
        outDict['beamFSMn_{0:02d}'.format(i)] = beam_imags[i]  
    return outDict  

""""波动传播模块"""
def run_process_hybr(beamLine):
    print('***************************')
    print('run_process_hybr  user')
    if isMono:
        fixedEnergy = E0
    else:
        fixedEnergy = False
    waveOnSlit = beamLine.slit_WB.prepare_wave(beamLine.source, nrays)
    beamSource = beamLine.source.shine(wave=waveOnSlit, fixedEnergy=fixedEnergy)  # 直接得到白光狭缝slit后的光场分布
    beamFSM0 = waveOnSlit

    # 单色光狭缝到 -> 白光准直镜：几何追迹
    WMGlobal, WMLocal = beamLine.WMirr.reflect(beam=beamSource)

    #白光准直镜 -> 单色器一晶：波动传播
    if False:
        waveOnDCM01 = beamLine.DCM01.prepare_wave(beamLine.WMirr, nrays)
        beamToDCM01 = rw.diffract(WMLocal, waveOnDCM01)
        DCM01Global, DCM01Local = beamLine.DCM01.reflect(
            beamToDCM01, noIntersectionSearch=True)    
    else:
        DCM01Global, DCM01Local = beamLine.DCM01.reflect(WMGlobal)
        
    #单色器一晶 -> 单色器二晶 -> 单色狭缝：几何追迹
    DCM02Global, DCM02Local = beamLine.DCM02.reflect(beam=DCM01Global)
    beamLine.fsm1.center = beamLine.slit_DCM.center
    beamFSM1 = beamLine.fsm1.expose(beam=DCM02Global)    
    slit_DCM_Local1b = beamLine.slit_DCM.propagate(beam=DCM02Global)
    DCM_OUT = DCM02Global

    #单色器二晶 -> TM：波动传播或者几何追迹
    """  Toroidal mirror  """      
    TMGlobal, TMLocal = beamLine.TM.reflect(DCM_OUT)
 
    """  Sample position  """  
    #1. 样品附近多观察点屏 
    
    #HKB -> 样品处：波动传播
    #1. Single wave - 样品处-焦点
    lsamp = [0, p_samp, TM_z]   
    beamLine.fsmsamp.center = lsamp
    waveOnSample = beamLine.fsmsamp.prepare_wave(beamLine.TM, fsmExpX, fsmExpZ)
    rw.diffract(TMLocal, waveOnSample)

    #2. Waves array @ some porition
    waveOnImags = []
    for dq in dqs:  # prepare the wave
        lsamp = [0, p_samp+dq, TM_z]     
        beamLine.fsmn.center = lsamp
        waveOnImag = beamLine.fsmn.prepare_wave(beamLine.TM, fsmExpX, fsmExpZ)
        rw.diffract(TMLocal, waveOnImag)
        waveOnImags.append(waveOnImag) 


    outDict = {
        'beamSource':beamSource,
        'beamFSM0':beamFSM0,
        'WMGlobal': WMGlobal,'WMLocal': WMLocal,
        'DCM':beamFSM1,
        'DCM01Global': DCM01Global,'DCM01Local': DCM01Local, 'DCM02Global': DCM02Global,'DCM02Local': DCM02Local,
        'TMGlobal': TMGlobal, 'TMLocal': TMLocal,
        'Samp': waveOnSample}  
    for ic, waveOnImag in enumerate(waveOnImags):
        outDict['beamFSMn_{0:02d}'.format(ic)] = waveOnImag     

    return outDict  

if what == 'rays':
    raycing.run.run_process = run_process_rays
elif what == 'hybrid':
    raycing.run.run_process = run_process_hybr

""""画图"""
plots = []
plots_FSMns = []
Indicator_plot = 1
if Indicator_plot == 1:
    #Source
    plot = xrtp.XYCPlot(
        'beamSource', aspect='auto',
        xaxis=xrtp.XYCAxis(r'$x$', unit=r"$\mu$m"),
        yaxis=xrtp.XYCAxis(r'$z$', unit=r"$\mu$m"),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0,
                            limits=[eMin0-1, eMax0+1]),
        title='Source size'
        )
    plot.xaxis.limits = [-30, 30]
    plot.yaxis.limits = [-30, 30]  
    plot.baseName = plot.title
    plots.append(plot)

    plot = xrtp.XYCPlot(
        'beamSource', aspect='auto',
        xaxis=xrtp.XYCAxis(r"x'", unit=r"$\mu$rad"),
        yaxis=xrtp.XYCAxis(r"z'", unit=r"$\mu$rad"),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0,
                            limits=[eMin0-1, eMax0+1]),
        title='source angle', fluxFormatStr=r"%g")
    plot.baseName = plot.title
    plots.append(plot)
    #WM
    plot = xrtp.XYCPlot(
        'WMLocal', aspect='auto',
        xaxis=xrtp.XYCAxis(r'$x$', unit=r"mm", bins=xbins, ppb=xppb),
        yaxis=xrtp.XYCAxis(r'$y$', unit=r"mm", bins=zbins, ppb=zppb),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0, limits=[eMin0-1, eMax0+1]),
        title='VCMLocal')    
    plot.baseName = plot.title
    plots.append(plot)
    #DCM
    plot = xrtp.XYCPlot(
        'DCM', aspect='auto',
        xaxis=xrtp.XYCAxis(r'$x$', unit=r"$\mu$m", bins=xbins, ppb=xppb),
        yaxis=xrtp.XYCAxis(r'$z$', unit=r"$\mu$m", bins=zbins, ppb=zppb),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0, limits=[eMin0-1, eMax0+1]),
        title='DCM02Global')    
    plot.baseName = plot.title
    plots.append(plot)
    #TM
    plot = xrtp.XYCPlot(
        'TMLocal', aspect='auto',
        xaxis=xrtp.XYCAxis(r'$x$', unit=r"$\mu$m", bins=xbins, ppb=xppb),
        yaxis=xrtp.XYCAxis(r'$y$', unit=r"$\mu$m", bins=zbins, ppb=zppb),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0, limits=[eMin0-1, eMax0+1]),
        title='TMLocal')    
    plot.baseName = plot.title
    plots.append(plot)

    #Sample
    plot = xrtp.XYCPlot(
        'Samp', aspect='auto',
        xaxis=xrtp.XYCAxis(r'$x$', unit=r"$\mu$m", bins=xbins, ppb=xppb),
        yaxis=xrtp.XYCAxis(r'$z$', unit=r"$\mu$m", bins=zbins, ppb=zppb),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0, limits=[eMin0-1, eMax0+1]),
        title='Sample')    
    plot.baseName = plot.title
    plot.xaxis.limits = [-x_lim, x_lim]
    plot.yaxis.limits = [-z_lim, z_lim]      
    plots.append(plot)
   
    #Imags @ different position
    for i, dq in enumerate(dqs):
        plot = xrtp.XYCPlot(
            'beamFSMn_{0:02d}'.format(i), aspect='auto',
            xaxis=xrtp.XYCAxis(
                r'$x$', r"$\mu$m", bins=xbins, ppb=xppb),#µm
            yaxis=xrtp.XYCAxis(
                r'$z$', r"$\mu$m", bins=zbins, ppb=zppb),#
            title = beamLine.fsmn.name+'_{0:02d}'.format(i),
            caxis=xrtp.XYCAxis('energy', 'eV',offset=E0,limits=[E0-2,E0+2])
            )
        plot.xaxis.fwhmFormatStr = '%.2f'
        plot.yaxis.fwhmFormatStr = '%.2f'
        plot.textPanel = plot.fig.text(
            0.2, 0.75, '', transform=plot.fig.transFigure, size=14, color='r',
            ha='left')
        plot.textPanelTemplate = 'd$q=${0:.3f} mm'.format(dq)
        plot.xaxis.limits = [-x_lim, x_lim]
        plot.yaxis.limits = [-z_lim, z_lim] 
        plots.append(plot)
        plots_FSMns.append(plot)

    for plot in plots:
        plot.saveName = [prefix + plot.title + '.png', ]
        if plot.caxis.label.startswith('energy'):
            plot.caxis.limits = eMin0, eMax0
            plot.caxis.offset = E0
        if plot.fluxKind.startswith('power'):
            plot.fluxFormatStr = '%.0f'
        else:
            plot.fluxFormatStr = '%.1p'

""" Generator """
def plot_generator(plots, plots_FSMns, beamLine):
    for plot in plots:
        plot.xaxis.fwhmFormatStr = '%.2f'
        plot.yaxis.fwhmFormatStr = '%.2f'         
        plot.fluxFormatStr = '%.2p'        
        fileName = '{0}'.format(plot.title)
        strTrajOutFileName = fileName + '.png'
        plot.saveName =  os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName)
        try:
            plot.textPanel.set_text(
                plot.textPanelTemplate)
        except AttributeError:
            pass

    yield
    #生成动态图像
    dynamic_images(plots_FSMns)
    defocus_plot(plots_FSMns, dqs)

""" afterScript """
def afterScript(plots):
    for plot in plots:
        save_data(plot)
    print('AfterScript done')

""""run program"""
def main():
    # xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine)

    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine, 
            generator=plot_generator, generatorArgs=[plots, plots_FSMns, beamLine],
            afterScript=afterScript, afterScriptArgs=[plots])

if __name__ == '__main__':
    main()     