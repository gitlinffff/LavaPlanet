from astropy.io import fits
import matplotlib.pyplot as plt
from numpy import exp,pi


def PlanckFunction(T,lmd):
    
    h = 6.62607015e-34   # Js
    c = 299792458        # m/s
    k = 1.380649e-23     # J/K
    
    B = (2*h*c*c)/(lmd**5)/(exp(h*c/k/lmd/T)-1)  # W m-2 sr-1 / m
    
    return B



spec = fits.getdata('hlsp_muscles_multi_multi_hd85512_broadband_v22_const-res-sed.fits',1)

fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.plot(spec['WAVELENGTH'],spec['FLUX'],label='HD85512 Spectra')

r_star = 0.533 * 6.957e8    # 0.533 * radius of the Sun 
d = 36.782 * 299792458 * 365 * 24 * 3600    # 36.782 light year
T = 4800 # temperature K
pl = PlanckFunction(T,spec['WAVELENGTH']*1e-10) * 1e-7 * pi * r_star**2 / d**2   # ergs / s / cm2 / Ang

ax.plot(spec['WAVELENGTH'],pl,label=f'{T} K Planck Function')

ax.set_xlabel('$Wavelength~(Angstroms)$',size=13)
ax.set_ylabel('$Flux~Density~(\mathrm{erg/{cm}^2/s/Ang})$',size=13)
ax.tick_params(axis='both', which='major', labelsize=13)

ax.legend(fontsize=13)

#plt.show()
fig.savefig("/home/linfel/LavaPlanet/images/star_spectra.png",dpi=300)
