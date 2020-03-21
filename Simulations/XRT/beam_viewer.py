#read pickle file exported by XRT tracing file
#for example: crlglobal2.export_beam('crlglobal2.pickle', fformat='pickle')

import pickle
import matplotlib.pyplot as plt
pickleName = 'crlglobal2.pickle'
with open(pickleName, 'rb') as f:
    dump = pickle.load(f) 
plt.figure(2)
x = dump['x']
z = dump['z']
state = dump['state']
plt.scatter(x[state==1],z[state==1])
plt.show()

