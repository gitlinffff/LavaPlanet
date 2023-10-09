# plot the histogram for line list

import numpy as np
import matplotlib.pyplot as plt

pwd = "/home/linfel/LavaPlanet/"
input_file = pwd + "co_HITRAN2020/lb_28Si-16O__SiOUVenIR__100-333__296K.par"

inten_list = []

wvn_lower = 100.0
wvn_upper = 250.0
with open(input_file, 'r') as infile:

    for i,line in enumerate(infile,1):
        wavenumber = float(line[3:15])
        intensity = float(line[15:25])
        if (wvn_lower <= wavenumber <= wvn_upper) & (intensity >= 1e-110):
            inten_list.append(intensity)

log10_inten_list = np.log10(inten_list)


fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
ax.hist(log10_inten_list, bins=100, edgecolor='black', alpha=0.7)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_yscale('log')    # use a logarithmic scale for the y-axis
ax.set_ylabel("Frequency",size=15)
ax.set_xlabel("$\log_{10}(intensity) log_{10}(\mathrm{cm}^{-1}/(\mathrm{molecule} \mathrm{cm}^{-2}))$",size=15)


output_fig = pwd + 'images/linelist_hist_'+str(wvn_lower)+'-'+str(wvn_upper)+'.png'
plt.savefig(output_fig,dpi=300)
print(f'# image saved as {output_fig}')
plt.close()

