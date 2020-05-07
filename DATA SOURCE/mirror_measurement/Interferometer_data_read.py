# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import csv

#%%
#读取干涉仪数据文件

x = []
y = []
z = [] 
n = 14
filename = "A0.xyz"
with open(filename, "r") as f:
    reader = csv.reader(f, delimiter = ' ')
    for i, row in enumerate(reader):
        if i<n:
            print(row)
        if i==3:
            nx, ny = int(row[2]),int(row[3])
        if (i>=n)&(row[0]!='#'):
           x.append(int(row[0]))
           y.append(int(row[1]))
           if row[2] == 'No':
               z.append(0)
           else:
               z.append(float(row[2]))
x = np.array(x)
y = np.array(y)
z = np.array(z)  
z_2D = z.reshape(nx,ny)*632;
x_1D = np.unique(x)
y_1D = np.unique(y)


#%%
#画图处理
plt.figure(figsize=(12,5))
plt.subplot(121)
im =plt.imshow(z_2D, interpolation='bicubic', cmap='jet',
                extent=[min(y_1D),max(y_1D),min(x_1D),max(x_1D)], aspect='auto')#,vmin=vmin, vmax=vmax)
plt.colorbar()
plt.title('Total Shape')
plt.subplot(122)
p_start = [264,400]
p_end = [780,652]
z_2D_part = z_2D[p_start[0]:p_end[0],p_start[1]:p_end[1]]
extent = [p_start[1], p_end[1], p_start[0], p_end[0]]
im =plt.imshow(z_2D_part, interpolation='bicubic', cmap='jet',
                extent=extent, aspect='auto')#,vmin=vmin, vmax=vmax)
plt.colorbar()
plt.title('Part Shape')
plt.savefig(filename+'profile_surface.png', dpi =300)