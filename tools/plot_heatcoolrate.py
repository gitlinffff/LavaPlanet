import numpy as np
import matplotlib.pyplot as plt

heat_rate_profile = np.genfromtxt("/home/linfel/LavaPlanet/outputs/heat_rate.csv", delimiter=',')
heat_rate_profile = np.flipud(heat_rate_profile)

h_atm = 100.0  # thickness of atmosphere
nlayers = heat_rate_profile.shape[0]
dh = h_atm/nlayers
dx = 1
h_grid, x_grid = np.mgrid[(dh/2):(h_atm-dh/2)+dh:dh, 1:7+dx:dx]


plt.figure(figsize=(10, 6))
heatmap = plt.pcolor(x_grid,h_grid,heat_rate_profile,cmap='RdBu_r', vmin=-np.max(np.abs(heat_rate_profile)), vmax=np.max(np.abs(heat_rate_profile)))
cb = plt.colorbar(heatmap)
cb.set_label('heating rate (K/s)',size=13)
plt.xlabel('wavenumber / $\mathrm{cm}^{-1}$',size=13)
plt.ylabel('height / km',size=13)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.savefig("/home/linfel/LavaPlanet/images/heating_rate.png",dpi=300)
