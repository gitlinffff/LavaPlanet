import numpy as np
import matplotlib.pyplot as plt

pwd = "/home/linfel/LavaPlanet/co_HITRAN2020/"
input_file = pwd + "lb_28Si-16O__SiOUVenIR__100-14285__296K.par"

wvn_list = []
inten_list = []

with open(input_file, 'r') as infile:

    for i,line in enumerate(infile,1):
        wavenumber = float(line[3:15])
        intensity = float(line[15:25])
        wvn_list.append(wavenumber)
        inten_list.append(intensity)


fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
ax.plot(wvn_list,inten_list, linewidth=0.5)

ax.set_yscale('log')    # use a logarithmic scale for the y-axis
ax.set_xlim(2222, 3125)  # Limit the x-axis range
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel("wavenumber $\mathrm{cm}^{-1}$",size=15)
ax.set_ylabel("intensity $\mathrm{cm}^{-1}/(\mathrm{molecule} \mathrm{cm}^{-2})$",size=15)
ax.legend(fontsize = 15,markerscale=1.5)


figname = 'SiO_spectrum_2222-3125.png'
plt.savefig(pwd + figname,dpi=300)
plt.close()

