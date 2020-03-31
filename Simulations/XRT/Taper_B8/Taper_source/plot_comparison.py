# -*- coding: utf-8 -*-

#%%
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

data_XRT = np.loadtxt('XRT_data.txt',skiprows=2)
data_SPECTRA = np.loadtxt('SPECTRA.dc0', skiprows=2, usecols=[0, 1], unpack=True)
data_SPECTRA = data_SPECTRA.T

fig1 = plt.figure(1, figsize=(7, 5))
ax = plt.subplot(111, label='1')
ax.set_xlabel(u'energy (keV)')
ax.set_ylabel(u'flux (a.u.) ')
ax.plot(data_XRT[:,0]*1e-3, data_XRT[:,1]/max(data_XRT[:,1]), 'g', label='calculated by XRT', lw=2)
ax.plot(data_SPECTRA[:,0]*1e-3, data_SPECTRA[:,1]/max(data_SPECTRA[:,1]), 'b', label='calculated by Spectra', lw=2)
ax.legend(loc='upper center')
fig1.savefig('compareTaper1.png',dpi = 300)
plt.show()

fig2 = plt.figure(2, figsize=(7, 5))
ax = plt.subplot(111, label='1')
ax.set_xlabel(u'energy (keV)')
ax.set_ylabel(u'flux (ph/s/0.1% bw) ')
ax.plot(data_XRT[:,0]*1e-3, data_XRT[:,1], 'g', label='calculated by XRT', lw=2)
ax.plot(data_SPECTRA[:,0]*1e-3, data_SPECTRA[:,1], 'b', label='calculated by Spectra', lw=2)
ax.legend(loc='upper center')
fig1.savefig('compareTaper2.png',dpi = 300)
plt.show()



# %%
