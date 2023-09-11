#!/usr/bin/env python
# coding: utf-8

"""plot absorption cross-section """

from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
import numpy as np

# create the plot for absorption cross-section of SiO
fig, ax = plt.subplots(3, 1, figsize=(10, 18))

avo = 6.02214076e23


for i in range(7):

    # specify input file 
    filename = "kcoeff.lava_planet-B" + str(i+1) + ".nc"
    pwd = "/home/linfel/LavaPlanet/"
    kcoeff_fpath = pwd + filename

    # read the data
    kcoeff_dataset = Dataset(kcoeff_fpath,'r')
    wavenumber = kcoeff_dataset.variables['Wavenumber'][:]
    pressure = kcoeff_dataset.variables['Pressure'][:]
    T_SiO = kcoeff_dataset.variables['T_SiO'][:]
    temperature = kcoeff_dataset.variables['Temperature'][:]
    SiO = kcoeff_dataset.variables['SiO'][:]
    
    # set mask for invalid values
    mask = SiO < np.log(1000 * avo * 1e-50)
    SiO_masked = np.ma.masked_array(SiO, mask=mask)

    # kcoeff at three different pressure(height)
    kcoeff_TOA = SiO_masked[:,0,1]
    kcoeff_MID = SiO_masked[:,39,1]
    kcoeff_BOT = SiO_masked[:,79,1]

    x = 1/wavenumber*10**4
    y1 = np.log10(np.exp(kcoeff_TOA)/1000/avo)
    y2 = np.log10(np.exp(kcoeff_MID)/1000/avo)
    y3 = np.log10(np.exp(kcoeff_BOT)/1000/avo)

    ax[0].plot(x,y1, linewidth=1,label=str(round(pressure[0])) + ' Pa')
    ax[1].plot(x,y2, linewidth=1,label=str(round(pressure[39])) + ' Pa')
    ax[2].plot(x,y3, linewidth=1,label=str(round(pressure[79])) + ' Pa')

# set range of x-axis
#x_range = (12,12.5)

for a in ax:
   
    #a.set_xlim(x_range[0],x_range[1])   # Limit the x-axis range
    a.set_xscale('log')  # use a logarithmic scale for the x-axis
    a.tick_params(axis='both', which='major', labelsize=15)
    a.set_xlabel("wavelength (\u03BCm)",size=15)
    a.set_ylabel("$\log_{10}(\mathrm{m}^2/\mathrm{molecule})$",size=15)
    a.legend(fontsize = 15,markerscale=1.5)


#figname = 'SiO_multibands_' + str(x_range[0]) + '-' + str(x_range[1]) + 'um.png'
figname = 'SiO_multibands.png'
plt.savefig(pwd + 'images/' + figname,dpi=300)
plt.close()
