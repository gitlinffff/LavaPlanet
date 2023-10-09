# plot histogram for line list based on S/S_mean, discard noise data

import numpy as np
import matplotlib.pyplot as plt

pwd = "/home/linfel/LavaPlanet/"
input_file = pwd + "co_HITRAN2020/lb_28Si-16O__SiOUVenIR__100-333__296K.par"


wvn_lower = 100.0
wvn_upper = 250.0
inten_thres = 1e-110

inten_list = []
with open(input_file, 'r') as infile:

    for i,line in enumerate(infile,1):
        wavenumber = float(line[3:15])
        intensity = float(line[15:25])
        if (wvn_lower <= wavenumber <= wvn_upper) & (intensity >= inten_thres):
            inten_list.append(intensity)

S_mean = np.mean(inten_list)
print(f"mean intensity {S_mean}")
S_over_Smean = inten_list/S_mean
log_S_over_Smean = np.log10(S_over_Smean)

fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
ax.hist(log_S_over_Smean, bins=100, edgecolor='black', alpha=0.7)
ax.set_xlim(left=-10)
ax.tick_params(axis='both', which='major', labelsize=15)
#ax.set_xscale('log')    # use a logarithmic scale for the x-axis
ax.set_yscale('log')    # use a logarithmic scale for the y-axis
ax.set_ylabel("Frequency",size=15)
ax.set_xlabel("$log_{10}(\mathrm{S}/\overline{\mathrm{S}})$",size=15)



output_fig = pwd + 'images/linelist_hist_'+str(wvn_lower)+'-'+str(wvn_upper)+'.png'
plt.savefig(output_fig,dpi=300)
print(f'# image saved as {output_fig}')
plt.close()

