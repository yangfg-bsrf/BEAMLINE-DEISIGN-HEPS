#从xrt的追迹结果，提取光斑的宽度信息，并画图

import matplotlib.pyplot as plt
import numpy as np

def defocus_plot(Plots, dqs): 
    plt.figure(figsize=(7, 5), dpi=72)
    ax1 = plt.subplot(111)
    ax1.set_title(r'FWHM size of beam cross-section near focal position')
    ax1.set_xlabel(r'd$q$ (mm)', fontsize=14)
    ax1.set_ylabel(u'FWHM size (µm)', fontsize=14)    
    qCurve = []
    xCurve = []
    zCurve = []
    for dq, plot in zip(dqs, Plots):
        qCurve.append(dq)            
        xCurve.append(plot.dx)
        zCurve.append(plot.dy)
    ax1.plot(qCurve, xCurve, 'o', label='Horizontal direction')  
    ax1.plot(qCurve, zCurve, '+', label='Vertical direction')
    ax1.set_xlabel(u'Shift (mm)', fontsize=14)
    ax1.set_ylabel(u'FWHM size (µm)', fontsize=14)
    plt.legend()          
    plt.show()
    plt.savefig('Depth_Of_Focus.png', dpi = 200)
    print(plot.title," Defocus Data save :Done")