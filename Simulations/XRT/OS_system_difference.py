#根据linux系统还是windows，选择是否使用交互式画图
import os
import matplotlib
if os.name == 'nt':
    print('the OS is Windows, interative plot is used')

if os.name == 'posix':
    print('the OS is Linux, interative plot is not used')
    matplotlib.use('Agg')
import matplotlib.pyplot as plt    

#test example
import numpy
plt.figure()
x0 = numpy.linspace(0,10,100)
plt.plot(x0,x0)
plt.show()
plt.savefig('test.png', dpi=300)
