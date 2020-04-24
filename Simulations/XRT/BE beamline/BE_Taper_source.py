# -*- coding: utf-8 -*-
"""
2020.1.30：
    修改1：增添几何光源
    修改2：增添KB镜聚焦
2020.2.1：   
    修改1：增加波动计算
2020.3.17：
    使用几何光源，研究离焦误差的影响
2020.4.24：
    BE 束线，检验光源的参数性能
"""
__author__ = "Fugui Yang"
__date__ = "17 April 2020"
import_lib = 1
if import_lib == 1:
    import os, sys
    import pickle, imageio, time
    import matplotlib.pyplot as plt
    import numpy as np
    # xrt函数引用
    import xrt.backends.raycing as raycing
    import xrt.runner as xrtr
    import xrt.plotter as xrtp
    import Distorted_Element as Distorted
    import xrt.backends.raycing.apertures as ra
    import xrt.backends.raycing.screens as rsc
    import xrt.backends.raycing.waves as rw

strExDataFolderName = 'output' #example data sub-folder name    
#从保存的图片中读取图片，并保存动态图片
def dynamic_images(Plots):
    
    gif_images = []    
    for plot in Plots:
        filename = plot.saveName
        gif_images.append(imageio.imread(filename)) 

    path =  os.path.join(os.getcwd(), strExDataFolderName, "test.gif")
    imageio.mimsave(path, gif_images,fps=1.5) 
    print("dynamic images : Done!")   

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
    plt.legend()          
    plt.show()
    path =  os.path.join(os.getcwd(), strExDataFolderName, "Depth_Of_Focus.png")
    plt.savefig(path, dpi = 200)

    plt.figure(300)
    ax1 = plt.subplot(111)
    ax1.plot(qCurve, cxCurve, 'o', label='Horizontal direction')  
    ax1.plot(qCurve, czCurve, '+', label='Vertical direction')
    ax1.set_xlabel(u'Shift (mm)', fontsize=14)
    ax1.set_ylabel(u'Position ($\mu$m)', fontsize=14)
    plt.legend()          
    plt.show()
    path =  os.path.join(os.getcwd(), strExDataFolderName, "position.png")
    plt.savefig(path, dpi = 200)
    print(plot.title," Defocus Data save :Done")
    
#画图数据保存
def save_data(plot):
    dump = []
    x = plot.xaxis.binCenters
    y = plot.yaxis.binCenters
    Ixy = plot.total2D
    dump.append([x, y, Ixy])
    pickleName = '{0}.pickle'.format(plot.title)
    path =  os.path.join(os.getcwd(), strExDataFolderName, pickleName)
    with open(path, 'wb') as f:
        pickle.dump(dump, f, protocol=2)
    print(plot.title," Data save :Done")
    
#画图数据保存
#**********************************************************#
raycing.precisionOpenCL = 'float64'
raycing.targetOpenCL = "auto"

""" ********************材料特性******************** """
Si = raycing.materials.Material(elements=r"Si",kind=r"mirror",rho=2.33,name=None)
crystalSi01 = raycing.materials.CrystalSi(name=None)

Error_path = '../surface error_DATA/'

""" ********************全局角度******************** """
accept_ang_x, accept_ang_z = 40e-6, 40e-6

""" 光源 """
Indicator_source = 1
if Indicator_source == 1:
    E0, n0 = 8000, 1  # 插入间工作能量点，及谐波级次
    n1 = 1  # 待分析谐波的能量点
    E1 = E0/n0*n1  # 待分析的能量点
    nrays, repeats = 1e5, 20
    what = 'rays'
#    what = 'hybrid'
    if what == 'rays':
        uniformRayDensity = False
        filamentBeam = False
        isMono = False
        prefix = 'Taper-rays-'
    else:
        uniformRayDensity = True
        filamentBeam = True
        isMono = False
        prefix = 'Taper-hybr-'

    dE = E1*0.0005/2
    eMin0, eMax0 = E1-dE, E1+dE
    kwargs_SR = dict(
        eE = 6.0, eI = 0.2,
        eEspread = 0.00111,
        eEpsilonX = 0.02755, eEpsilonZ = 0.002755,  
#        eEspread = 0,
#        eEpsilonX = 0, eEpsilonZ = 0,         
        betaX = 10.12, betaZ = 9.64,
        eMin=eMin0, eMax=eMax0, eN=51, #targetE = [10000,3],
        distE='BW',  
        xPrimeMax = accept_ang_x, zPrimeMax = accept_ang_z,
        xPrimeMaxAutoReduce = False, zPrimeMaxAutoReduce = False, 
        uniformRayDensity = uniformRayDensity,
        filamentBeam = filamentBeam,
        K = 0.75769, period=32.7,n=151,taper=(4, 23.17))

    kwargs_Geo = dict(
        name = 'geometric source', nrays = nrays, center = [0, 0, 0],
        distx='normal', dx=0.02, distz='normal', dz=0.02, 
        distxprime='flat', dxprime=accept_ang_x, distzprime='flat', dzprime=accept_ang_z, 
        distE='flat', energies=(eMin0,  eMax0),
        # distE='lines', energies=(E0,), 
        polarization='horizontal', 
        filamentBeam=filamentBeam, uniformRayDensity=uniformRayDensity)

    case = 'point source'
    if case == 'point source':  # point source
        kwargs_Geo.update(dict(
            dx=0.00, dz=0.00, distxprime='flat', 
            distzprime='flat'))
    elif case == 'collimated source':  # collimated source
        kwargs_Geo.update(dict(
            dx=1, dz=1, distx='flat', distz='flat',
            distxprime=None, distzprime=None)) 
 
    beamLine = raycing.BeamLine(alignE=E0)
    """  选择光源类型  """
    source_und = raycing.sources.Undulator(beamLine, **kwargs_SR)
    # source_geo = raycing.sources.GeometricSource(beamLine, **kwargs_Geo)
    beamLine.source = source_und

""" 第一光学元件：白光狭缝 """
Indicator_slit_WB = 1
if Indicator_slit_WB ==1:
    p_slitWB = 38000
    slitDx, slitDz = accept_ang_x*p_slitWB, accept_ang_z*p_slitWB
    opening = [-slitDx/2, slitDx/2, -slitDz/2, slitDz/2]
    beamLine.slit_WB = ra.RectangularAperture(
        beamLine, 'Slit_WB',  [0, p_slitWB, 0], 
        ('left', 'right', 'bottom', 'top'), opening = opening )

""" 第二光学元件：单色器 """
Indicator_DCM = 1
if Indicator_DCM ==1:
    p_DCM = 39000
    DCM_MODEL = 0  # 选择自己建模的DCM：0; 选择系统自带的DCM：1；
    DCM_error = 0
    Surface_name = 'DCM_OE2_'
    fixedOffset = 0.5
    DCM01x = 0
    p_DCM01 = p_DCM
    theta_DCM = crystalSi01.get_Bragg_angle(E0)-crystalSi01.get_dtheta_symmetric_Bragg(E0)
    DCM02x = DCM01x+fixedOffset/np.sin(2*theta_DCM)*np.sin(2*(theta_DCM))
    DCM02z = DCM01x+fixedOffset/np.sin(2*theta_DCM)*np.cos(2*(theta_DCM))
    p_DCM02 = p_DCM01+DCM02z

    # 使用传统单色器模型
    kwargs = dict(
        name='DCM',
        center=[0, p_DCM, 0],
        bragg=r"auto",  positionRoll=np.pi/2,
        material=crystalSi01, material2=crystalSi01,
        fixedOffset=fixedOffset)
    beamLine.dcM1 = raycing.oes.DCM(bl=beamLine, **kwargs)

    # 使用双平晶代表，在tracing部分的代码也是不同的
    kwargs_DCM01 = dict(
        name='DCM01', center=[DCM01x, p_DCM01, 0],
        pitch=theta_DCM,#[E0],
        positionRoll=np.pi/2,
        material=crystalSi01,
        limPhysX=[-25, 25], limPhysY=[-25, 25])
    kwargs_DCM02 = dict(
        name='DCM02', center=[DCM02x, p_DCM02, 0],
        pitch=-theta_DCM,#[E0],
        positionRoll=-np.pi/2,
        material=crystalSi01,
        limPhysX=[-25, 25], limPhysY=[-25, 25])
    if DCM_error == 0:
        # 不带面形误差
        beamLine.DCM01 = raycing.oes.OE(bl=beamLine, **kwargs_DCM01)
    else:
        # 带面形误差
        kwargs_DCM01['get_distorted_surface'] = 'error'
        kwargs_DCM01['fname1'] = Surface_name+'X.txt'
        kwargs_DCM01['fname2'] = Surface_name+'Y.txt'
        kwargs_DCM01['fname3'] = Surface_name+'Z.txt'
        beamLine.DCM01 = Distorted.PlaneMirrorDistorted(
            bl=beamLine, **kwargs_DCM01)

    beamLine.DCM02 = raycing.oes.OE(bl=beamLine, **kwargs_DCM02)

""" 第三光学元件：单色器狭缝 """
Indicator_slit_DCM = 1
if Indicator_slit_DCM == 1:
    p_slitDCM = p_DCM + 500
    beamLine.slit_DCM = raycing.apertures.RectangularAperture(
        beamLine, 'Slit_DCM',  [fixedOffset, p_slitDCM,
                                0], ('left', 'right', 'bottom', 'top'),
        [-slitDx/2, slitDx/2, -slitDz/2, slitDz/2])

""" 聚焦光学元件：KB聚焦镜  """ 
Indicator_KB = 1
if Indicator_KB == 1:
    theta_KB = 4e-3
    y_VKB = p_DCM + 1000
    p_VKB = y_VKB    
    q_VKB= 40000
    # y_samp = 40000
    
    oe_error_V = 0  
    Surface_name= Error_path+'VKB' 
    globalRoll = 0
    inclination = 0
    theta_VKB = theta_KB

    sourceCenter = [fixedOffset, 0, 0]
    kwargs_OE = dict(
                name ='VKB', 
                center=[fixedOffset, y_VKB, 0],
#                material = Si,
                limPhysX=[-1, 1], limPhysY=[-150, 150],
                rotationSequence='RyRzRx',
                pitch= theta_VKB+ inclination*np.cos(globalRoll), 
                positionRoll = globalRoll,
                yaw = inclination*np.sin(globalRoll),
                q = q_VKB, p = p_VKB,
                isCylindrical=True)    
    if oe_error_V==0:        
        beamLine.VKB = raycing.oes.EllipticalMirrorParam(bl=beamLine,**kwargs_OE)
    else:            
        kwargs_OE['get_distorted_surface'] ='error'
        kwargs_OE['fname1'] =Surface_name+'X.txt'
        kwargs_OE['fname2'] =Surface_name+'Y.txt'
        kwargs_OE['fname3'] =Surface_name+'Z.txt'       
        beamLine.VKB = Distorted.EllipMirrorDistorted(beamLine, **kwargs_OE)
        
    theta_HKB = theta_KB
    y_HKB = y_VKB + 500
    p_HKB = p_VKB + (y_HKB-y_VKB)/np.cos(2*theta_HKB)
    q_HKB = q_VKB - (y_HKB-y_VKB)/np.cos(2*theta_HKB)
    
    inclination = 2 * theta_VKB
    globalRoll = - np.pi/2
    sourceCenter = [fixedOffset, 0, -p_HKB*2 * theta_VKB] #
    oe_error_H = 0  
    Surface_name=Error_path+'KB38_'
    
    kwargs_OE = dict(
                name='HKB',
                center=[fixedOffset, y_HKB, np.tan(inclination)*(y_HKB-y_VKB)],
                pitch = theta_HKB+inclination*np.cos(globalRoll), 
                positionRoll=globalRoll,
                yaw = inclination*np.sin(globalRoll),                
                rotationSequence='RyRzRx',                                    
#                material = Si,
                limPhysX=[-1, 1], limPhysY=[-150, 150],
                q = q_HKB, p = p_HKB,
                isCylindrical=True)    
    if oe_error_H ==0:        
        beamLine.HKB = raycing.oes.EllipticalMirrorParam(bl=beamLine,**kwargs_OE)
    else:            
        kwargs_OE['get_distorted_surface'] ='error'
        kwargs_OE['fname1'] =Surface_name+'X.txt'
        kwargs_OE['fname2'] =Surface_name+'Y.txt'
        kwargs_OE['fname3'] =Surface_name+'Z.txt'       
        beamLine.HKB = Distorted.EllipMirrorDistorted(beamLine, **kwargs_OE)

""" 观察屏组件 """
Indicator_screen = 1
if Indicator_screen == 1:
    zbins, zppb = 128, 4 
    xbins, xppb = 128, 4

    beamLine.fsm0 = rsc.Screen(beamLine, 'FSM0')
    beamLine.fsm1 = rsc.Screen(beamLine, 'FSM1')
    beamLine.fsm2 = rsc.Screen(beamLine, 'FSM2')

    x_lim, z_lim = 60, 30
    edges = np.linspace(-x_lim, x_lim, xbins+1)*1e-3
    fsmExpX = (edges[:-1] + edges[1:]) * 0.5
    edges = np.linspace(-z_lim, z_lim, zbins+1)*1e-3
    fsmExpZ = (edges[:-1] + edges[1:]) * 0.5  

    dqs = np.linspace(-1000, 1000, 21)
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
    # beamSource = beamLine.source.shine()
    beamLine.fsm0.center = beamLine.slit_WB.center
    beamFSM0 = beamLine.fsm0.expose(beam=beamSource)

    # DCM
    DCM_in = beamSource
    DCM01Global, DCM01Local = beamLine.DCM01.reflect(beam=DCM_in)
    DCM02Global, DCM02Local = beamLine.DCM02.reflect(beam=DCM01Global)
    beamLine.slit_DCM.propagate(beam=DCM02Global)
    beamLine.fsm1.center = beamLine.slit_DCM.center
    beamFSM1 = beamLine.fsm1.expose(beam=DCM02Global)
    DCM_OUT = DCM02Global

    #KB
    KB_in = DCM_OUT
    VKBGlobal, VKBLocal = beamLine.VKB.reflect(KB_in)
    HKBGlobal, HKBLocal = beamLine.HKB.reflect(VKBGlobal)
  
    outDict = {'beamSource': beamSource,
        'beamFSM0':beamFSM0,
        'DCM01Global':DCM01Global,
        'DCM01Local': DCM01Global,
        'DCM02Global': DCM02Global,
        'DCM02Local': DCM02Global,
        'beamFSM1': beamFSM1,
        'VKBGlobal': VKBGlobal,
        'VKBLocal': VKBLocal,
        'HKBLocal': HKBLocal,
        'HKBGlobal': HKBGlobal}

    #1. 样品处 焦点
    t_out_HKB = np.array([-np.sin(2*theta_KB), (np.cos(2*theta_KB))**2, np.sin(4*theta_KB)/2])
    l_HKB = beamLine.HKB.center
    lsamp = l_HKB + t_out_HKB*q_HKB

    beamLine.fsm2.center = lsamp
    beam_samp =  beamLine.fsm2.expose(HKBGlobal)
    outDict['Samp'] = beam_samp

    #2. 多屏
    for i, dq in enumerate(dqs): 
        lsamp = l_HKB + t_out_HKB*(q_HKB+dq)     
        beamLine.fsmn.center = lsamp
        beam_imag = beamLine.fsmn.expose(HKBGlobal)
        outDict['beamFSMn_{0:02d}'.format(i)] = beam_imag
    
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
    beamLine.fsm0.center = beamLine.slit_WB.center
    beamFSM0 = waveOnSlit  

    # DCM
    DCM_in = beamSource
    DCM01Global, DCM01Local = beamLine.DCM01.reflect(beam=DCM_in)
    DCM02Global, DCM02Local = beamLine.DCM02.reflect(beam=DCM01Global)
    beamLine.slit_DCM.propagate(beam=DCM02Global)
    beamLine.fsm1.center = beamLine.slit_DCM.center
    beamFSM1 = beamLine.fsm1.expose(beam=DCM02Global)
    DCM_OUT = DCM02Global

    #单色器到KB镜的传播方式，波动还是几何
    KB_in = DCM_OUT
    VKBGlobal, VKBLocal = beamLine.VKB.reflect(KB_in)
    HKBGlobal, HKBLocal = beamLine.HKB.reflect(VKBGlobal)   

    outDict = {'beamSource': beamSource,
        'beamFSM0':beamFSM0,
        'DCM01Global':DCM01Global,
        'DCM01Local': DCM01Global,
        'DCM02Global': DCM02Global,
        'DCM02Local': DCM02Global,
        'beamFSM1': beamFSM1,
        'VKBGlobal': VKBGlobal,
        'VKBLocal': VKBLocal,
        'HKBLocal': HKBLocal,
        'HKBGlobal': HKBGlobal}

    #1. 样品处 焦点
    #1. Single wave - 样品处-焦点
    t_out_HKB = np.array([-np.sin(2*theta_KB), (np.cos(2*theta_KB))**2, np.sin(4*theta_KB)/2])
    l_HKB = beamLine.HKB.center
    lsamp = l_HKB + t_out_HKB*q_HKB

    beamLine.fsm2.center = lsamp
    waveOnSample = beamLine.fsm2.prepare_wave(beamLine.HKB, fsmExpX, fsmExpZ)
    rw.diffract(HKBLocal, waveOnSample)  
    outDict['Samp'] = waveOnSample
    
    #2. 多屏
    waveOnImags=[]
    for i0, dq in enumerate(dqs):  # prepare the wave
        lsamp = l_HKB + t_out_HKB*(q_HKB+dq)
        beamLine.fsmn.center = lsamp
        waveOnImag = beamLine.fsmn.prepare_wave(beamLine.HKB, fsmExpX, fsmExpZ)
        rw.diffract(HKBLocal, waveOnImag)
        waveOnImags.append(waveOnImag) 
        outDict['beamFSMn_{0:02d}'.format(i0)] = waveOnImag

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
    plot.baseName = plot.title
    plots.append(plot)
    plot.xaxis.limits = [-50, 50]
    plot.yaxis.limits = [-30, 30]  
    
    plot = xrtp.XYCPlot(
        'beamSource', aspect='auto',
        xaxis=xrtp.XYCAxis(r"x'", unit=r"$\mu$rad"),
        yaxis=xrtp.XYCAxis(r"z'", unit=r"$\mu$rad"),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0,
                            limits=[eMin0-1, eMax0+1]),
        title='source angle')
    plot.baseName = plot.title
    plots.append(plot)
    plot.xaxis.limits = [-40, 40]
    plot.yaxis.limits = [-40, 40]   
    #Monitor
    plot = xrtp.XYCPlot(
        'DCM02Global', aspect='auto',
        xaxis=xrtp.XYCAxis(r'$x$', unit=r"$\mu$m"),
        yaxis=xrtp.XYCAxis(r'$z$', unit=r"$\mu$m"),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0, limits=[eMin0-1, eMax0+1]),
        title='DCM')
#    plot.xaxis.limits = [-2, 2]
#    plot.yaxis.limits = [-300, 300]
    plot.baseName = plot.title
    plots.append(plot)
    
    #Monitor
    plot = xrtp.XYCPlot(
        'HKBLocal', aspect='auto',
        xaxis=xrtp.XYCAxis(r'$x$', unit=r"mm"),
        yaxis=xrtp.XYCAxis(r'$y$', unit=r"mm"),
        caxis=xrtp.XYCAxis('energy', 'eV', offset=E0, limits=[eMin0-1, eMax0+1]),
        title='VKBLocal')
    plot.xaxis.limits = [-2, 2]
    plot.yaxis.limits = [-300, 300]
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
    plots.append(plot)
    plot.xaxis.limits = [-x_lim, x_lim]
    plot.yaxis.limits = [-z_lim, z_lim]   


    #Imags @ different position
    for i, dq in enumerate(dqs):
        plot = xrtp.XYCPlot(
            'beamFSMn_{0:02d}'.format(i), aspect='auto',
            xaxis=xrtp.XYCAxis(r'$x$', u'$\mu$m', bins=xbins, ppb=xppb),#
            yaxis=xrtp.XYCAxis(r'$z$', u'$\mu$m', bins=zbins, ppb=zppb),#
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
    print('plot_generator')
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
    """
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine)
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine, 
            generator=plot_generator, generatorArgs=[plots, plots_FSMns, beamLine])    
    """    
    xrtr.run_ray_tracing(plots, repeats=repeats, beamLine=beamLine, 
            generator=plot_generator, generatorArgs=[plots, plots_FSMns, beamLine],
            afterScript=afterScript, afterScriptArgs=[plots])



if __name__ == '__main__':
    main()
