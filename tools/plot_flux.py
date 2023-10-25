import numpy as np
import matplotlib.pyplot as plt

h_atm = 100.0  # thickness of atmosphere
layers = 79
dh = h_atm/layers
dx = 1
h, x = np.mgrid[0:h_atm+dh:dh, 1:7+dx:dx]

# plot the downward flux profile
flux_down = np.genfromtxt("/home/linfel/LavaPlanet/outputs/flux_down.csv", delimiter=',')
flux_down = np.flipud(flux_down)

plt.figure(figsize=(10, 6))
heatmap = plt.pcolor(x,h,flux_down,cmap='viridis')
cb = plt.colorbar(heatmap)
cb.set_label('flux',size=13)
plt.xlabel('wavenumber / $\mathrm{cm}^{-1}$',size=13)
plt.ylabel('height / km',size=13)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.savefig("/home/linfel/LavaPlanet/images/flux_down.png",dpi=300)

# plot the upward flux profile
flux_up = np.genfromtxt("/home/linfel/LavaPlanet/outputs/flux_up.csv", delimiter=',')
flux_up = np.flipud(flux_up)

plt.figure(figsize=(10, 6))
heatmap = plt.pcolor(x,h,flux_up,cmap='viridis')
cb = plt.colorbar(heatmap)
cb.set_label('flux',size=13)
plt.xlabel('wavenumber / $\mathrm{cm}^{-1}$',size=13)
plt.ylabel('height / km',size=13)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.savefig("/home/linfel/LavaPlanet/images/flux_up.png",dpi=300)
