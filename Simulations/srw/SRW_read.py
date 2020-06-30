# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 19:24:59 2020

@author: yangfg
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
import scipy.signal as signal

filename = 'B8_42.dat'
I = []
with open(filename, "r") as f:
    reader = csv.reader(f, delimiter = ' ')
    for i, row in enumerate(reader):
        if i<10:
            if i == 4:
                x_min = float(row[0].strip('#'))
            if i == 5:
                x_max = float(row[0].strip('#'))            
            if i == 6:
                Nx = int(row[0].strip('#'))  

            if i == 7:
                y_min = float(row[0].strip('#'))
            if i == 8:
                y_max = float(row[0].strip('#'))            
            if i == 9:
                Ny = int(row[0].strip('#'))
        if i>10:
            I.append(float(row[1]))
I0 = np.array(I)
I_2D = I0.reshape(Ny, Nx)           
x0 = np.linspace(x_min, x_max, Nx)*1e6
y0 = np.linspace(y_min, y_max, Ny)*1e6

plt.figure(1)
# plt.plot(x0, I_2D[Ny//2,:])
plt.plot(y0, I_2D[:, Nx//2])




filename = 'B8_45.dat'
I = []
with open(filename, "r") as f:
    reader = csv.reader(f, delimiter = ' ')
    for i, row in enumerate(reader):
        if i<10:
            if i == 4:
                x_min = float(row[0].strip('#'))
            if i == 5:
                x_max = float(row[0].strip('#'))            
            if i == 6:
                Nx = int(row[0].strip('#'))  

            if i == 7:
                y_min = float(row[0].strip('#'))
            if i == 8:
                y_max = float(row[0].strip('#'))            
            if i == 9:
                Ny = int(row[0].strip('#'))
        if i>10:
            I.append(float(row[1]))
I0 = np.array(I)
I_2D = I0.reshape(Ny, Nx)           
x0 = np.linspace(x_min, x_max, Nx)*1e6
y0 = np.linspace(y_min, y_max, Ny)*1e6

# plt.plot(x0, I_2D[Ny//2,:])
plt.plot(y0, I_2D[:, Nx//2])

# plt.figure(1)
# x = x0[abs(x0)<0.5]
# y = y0[abs(y0)<0.5]

# im =plt.imshow(I_2D, interpolation='bicubic', cmap='jet',
#                extent=[min(x0),max(x0),min(y0),max(y0)], aspect='auto')#,vmin=vmin, vmax=vmax)
# plt.axvline(x=0.0,color = 'r', linestyle = '--')
# plt.xlim([-700, 700])
# plt.ylim([-20, 20])
# plt.xlabel('x(um)')
# plt.ylabel('z(um)')
# plt.colorbar(im, shrink = 1)
# plt.show()  
# plt.savefig('2D distribution.png',dpi = 300)

