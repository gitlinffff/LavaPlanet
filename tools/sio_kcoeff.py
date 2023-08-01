#!/usr/bin/env python
# coding: utf-8

from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import numpy as np

# read the data
kcoeff_fpath = "/home/linfel/LavaPlanet/kcoeff.100-200-0.01.nc"
kcoeff_dataset = Dataset(kcoeff_fpath,'r')
wavenumber = kcoeff_dataset.variables['Wavenumber'][:]
pressure = kcoeff_dataset.variables['Pressure'][:]
T_SiO = kcoeff_dataset.variables['T_SiO'][:]
temperature = kcoeff_dataset.variables['Temperature'][:]
SiO = kcoeff_dataset.variables['SiO'][:]

kcoeff_TOA = SiO[:,0,1]
kcoeff_MID = SiO[:,39,1]
kcoeff_BOT = SiO[:,79,1]

# create the plot for absorption cross-section of SiO
fig, ax = plt.subplots(3, 1, figsize=(10, 18))

avo = 6.02214076e23

x = 1/wavenumber*10**4
y1 = np.log10(np.exp(kcoeff_TOA)/1000/avo)
y2 = np.log10(np.exp(kcoeff_MID)/1000/avo)
y3 = np.log10(np.exp(kcoeff_BOT)/1000/avo)

ax[0].plot(x,y1, linewidth=1, marker ='*',markersize = 3,label=str(round(pressure[0])) + ' Pa')
ax[1].plot(x,y2, linewidth=1, marker ='*',markersize = 3,label=str(round(pressure[39])) + ' Pa')
ax[2].plot(x,y3, linewidth=1, marker ='*',markersize = 3,label=str(round(pressure[79])) + ' Pa')

for a in ax:
    a.tick_params(axis='both', which='major', labelsize=15)
    a.set_xlabel("wavelength (\u03BCm)",size=15)
    a.set_ylabel("$\log_{10}(\mathrm{m}^2/\mathrm{molecule})$",size=15)
    a.legend(fontsize = 15,markerscale=1.5)

plt.savefig('sio_kcoeff.png',dpi=300)
plt.close()

#plt.title('Spotter vs. WAVEWATCH Ⅲ',fontsize = 14)
