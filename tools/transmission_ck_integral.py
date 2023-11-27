from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

def get_band_tau_profile(filename,h):

    # specify input file
    pwd = "/home/linfel/LavaPlanet/2800isothermal/"
    cktable_fpath = pwd + filename
    print(f"# working on {cktable_fpath}")

    # read the data
    ck_dataset = Dataset(cktable_fpath,'r')
    ck_g = ck_dataset.variables['Wavenumber'][:]
    Temperature = ck_dataset.variables['Temperature'][:]
    Pressure = ck_dataset.variables['Pressure'][:]

    nlayers = len(Pressure)

    tau_layers = np.array([])      # create an array to store tau of each layer
    trm_layers = np.array([])

    for ilayer in range(nlayers):
        ck_k_raw = ck_dataset.variables['SiO'][:,ilayer,1]  # k at specific pressure and temperature
        ck_k = np.exp(ck_k_raw)/1000  # m2/mol
        T = Temperature[ilayer]       # K
        p = Pressure[ilayer]          # Pa

        rou = p/T/8.31446261815324    # mol/m3
        d = h/nlayers*1000              # thickness of each layer (m) 
        u = rou * d                   # u is constant within each layer

        # integrate over g space to calculate transmission
        trm = 0.0
        for j in range(len(ck_g)-1):
            k = (ck_k[j] + ck_k[j+1]) / 2
            delta_g = ck_g[j+1] - ck_g[j]
            trm = trm + np.exp(-k*u) * delta_g
        trm_layers = np.append(trm_layers,trm)
        
        # calculate tau for this layer
        itau = -np.log(trm)
        tau_layers = np.append(tau_layers,itau)
    
    return trm_layers, tau_layers, nlayers

# specify the thickness of atmosphere (km)
h_atm = 100.0

# number of bands
nbands = 8

tau_list = []
trm_list = []
# get the tau profile for each band and add to the 2D array
for filenum in list(range(1,nbands+1)):

    filename = f"cktable.lava_planet-B{filenum}.nc"
    band_trm,band_tau,layers = get_band_tau_profile(filename,h_atm)
    tau_list.append(band_tau)
    trm_list.append(band_trm)

tau_all = np.vstack(tau_list).T
trm_all = np.vstack(trm_list).T

# write the tau profile in csv file
np.savetxt("/home/linfel/LavaPlanet/outputs/tau.csv", tau_all, delimiter=',')

# flip the tau and transmission profile to make plot
tau_all = np.flipud(tau_all)
trm_all = np.flipud(trm_all)

# create the grid to plot
dh = h_atm/layers
dx = 1
h_grid, x_grid = np.mgrid[(dh/2):(h_atm-dh/2)+dh:dh, 1:nbands+dx:dx]

# plot the transmission profile
plt.figure(figsize=(10, 6))
heatmap = plt.pcolor(x_grid,h_grid,trm_all,cmap='viridis')
cb = plt.colorbar(heatmap)
cb.set_label('transmission',size=13)
plt.xlabel('wavenumber / $\mathrm{cm}^{-1}$',size=13)
plt.ylabel('height / km',size=13)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.savefig("/home/linfel/LavaPlanet/images/trm_profile.png",dpi=300)

# plot the optical depth profile
plt.figure(figsize=(10, 6))
heatmap = plt.pcolor(x_grid,h_grid,tau_all,cmap='YlOrRd')
cb = plt.colorbar(heatmap)
cb.set_label('optical depth',size=13)
plt.xlabel('wavenumber / $\mathrm{cm}^{-1}$',size=13)
plt.ylabel('height / km',size=13)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.savefig("/home/linfel/LavaPlanet/images/tau_profile.png",dpi=300)
