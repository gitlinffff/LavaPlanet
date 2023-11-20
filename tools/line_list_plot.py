# plot line list and calculate mean line spacing

import numpy as np
import matplotlib.pyplot as plt

pwd = "/home/linfel/LavaPlanet/"
input_file = pwd + "co_HITRAN2020/lb_28Si-16O__SiOUVenIR__14285-25000__296K.par"

# set wavenumber range
wvn_lower = 14285.0
wvn_upper = 25000.0

# set intensity threshold above which a data is considered as an effective line
inten_thres = 1e-110

wvn_list = []
inten_list = []
cnt_line = 0
with open(input_file, 'r') as infile:

    for i,line in enumerate(infile,1):
        wavenumber = float(line[3:15])
        if wvn_lower <= wavenumber <= wvn_upper:
            intensity = float(line[15:25])
            wvn_list.append(wavenumber)
            inten_list.append(intensity)

            if intensity >= inten_thres:
                cnt_line += 1

# calculate line spacing
print(cnt_line)
d = (wvn_upper - wvn_lower)/cnt_line
print(f"# mean line spacing: {d} cm-1")

# create the plot
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
ax.plot(wvn_list,inten_list, linewidth=0.5)
#ax.plot(wvn_list,inten_list, marker='.', markersize=1, linestyle='None')
ax.set_yscale('log')    # use a logarithmic scale for the y-axis
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel("wavenumber $\mathrm{cm}^{-1}$",size=15)
ax.set_ylabel("intensity $\mathrm{cm}^{-1}/(\mathrm{molecule }\mathrm{cm}^{-2})$",size=15)
#ax.legend(fontsize = 15,markerscale=1.5)


output_fig = pwd + 'images/SiO_spectrum_'+str(wvn_lower)+'-'+str(wvn_upper)+'.png'
plt.savefig(output_fig,dpi=300)
print(f'# image saved as {output_fig}')
plt.close()

