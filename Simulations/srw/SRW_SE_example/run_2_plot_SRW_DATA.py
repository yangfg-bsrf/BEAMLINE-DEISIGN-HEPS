# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pickle

name = '2D_section_B8_SSA'
pickleName = name + '.pickle'
with open(pickleName, 'rb') as f:
    dump = pickle.load(f)   
x, z, y, Ix, Iz = dump

result1 = np.zeros([len(y),len(x)])
for i in range(len(Ix)):
    result1[i,:] =np.array(Ix[i])
    
T1 = np.log(result1.T)
#T1 = result1.T
z = z*1e3 
x = x*1e3
y = y*1e3

#vmin, vmax = 1e13, max(max(Iz))

plt.figure(1)
im =plt.imshow(T1, interpolation='bicubic', cmap='jet',
               extent=[min(y),max(y),min(x),max(x)], aspect='auto')#,vmin=vmin, vmax=vmax)
plt.axvline(x=0.0,color = 'r', linestyle = '--')
#plt.xlim([-1,1])
#plt.ylim([-2,2])
plt.xlabel('y(mm)')
plt.ylabel('x(um)')
plt.colorbar(im, shrink=0.5)
plt.show()  

plt.savefig('ideal_xy_Section.png',dpi = 300)


result2 = np.zeros([len(y),len(z)])
for i in range(len(Iz)):
    result2[i,:] =np.array(Iz[i])
    
T2 = np.log(result2.T)
#T2 = result2.T
plt.figure(2)
im =plt.imshow(T2, interpolation='bicubic', cmap='jet',
               extent=[min(y),max(y),min(z),max(z)], aspect='auto')#,vmin=vmin, vmax=vmax)
plt.axvline(x=0.0,color = 'r', linestyle = '--')

plt.xlabel('y(mm)')
plt.ylabel('z(um)')
plt.colorbar(im, shrink=0.5)
plt.show()  

plt.savefig('ideal_zy_Section.png',dpi = 300)



plt.figure(3)
plt.plot(x,result1[len(y)//2,:],label = 'horizontal')
plt.plot(z,result2[len(y)//2,:],label = 'vertical')
plt.xlabel('position(um)')
plt.ylabel('Intensity(a.u.)')  
#plt.xlim([-0.5,0.5])
plt.legend()
plt.show()
plt.savefig('ideal_yz_Section_1D.png',dpi = 300)
