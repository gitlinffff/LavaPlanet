import numpy as np
import matplotlib.pyplot as plt

temp_trk = np.genfromtxt("/home/linfel/LavaPlanet/rdsm/outputs/temp_trk.csv", delimiter=',')
temp_trk = np.flipud(temp_trk)

h_atm = 100.0  # thickness of atmosphere
nlayers = temp_trk.shape[0]
ncol = temp_trk.shape[1]


dh = h_atm/nlayers
dt = 1
h_grid, t_grid = np.mgrid[(dh/2):(h_atm-dh/2)+dh:dh, 0:(ncol-1)+dt:dt]


plt.figure(figsize=(10, 6))
heatmap = plt.pcolor(t_grid,h_grid,temp_trk,cmap='RdBu_r', vmin=-np.max(np.abs(temp_trk)), vmax=np.max(np.abs(temp_trk)))
cb = plt.colorbar(heatmap)
cb.set_label('temperature (K)',size=13)
plt.xlabel('time / $s$',size=13)
plt.ylabel('height / km',size=13)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.savefig("/home/linfel/LavaPlanet/images/temp_trk.png",dpi=300)
