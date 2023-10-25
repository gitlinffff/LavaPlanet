from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

def get_band_tau_profile(filename,h,layers):

    # specify input file
    pwd = "/home/linfel/LavaPlanet/"
    cktable_fpath = pwd + filename
    print(f"# working on {cktable_fpath}")

    # read the data
    ck_dataset = Dataset(cktable_fpath,'r')
    ck_g = ck_dataset.variables['Wavenumber'][:]
    Temperature = ck_dataset.variables['Temperature'][:]
    Pressure = ck_dataset.variables['Pressure'][:]

    tau_layers = np.array([])      # create an array to store tau of each layer
    trm_layers = np.array([])

    for ilayer in range(layers-1,-1,-1):
        ck_k_raw = ck_dataset.variables['SiO'][:,ilayer,1]  # k at specific pressure and temperature
        ck_k = np.exp(ck_k_raw)/1000  # m2/mol
        T = Temperature[ilayer]       # K
        p = Pressure[ilayer]          # Pa

        rou = p/T/8.31446261815324    # mol/m3
        d = h/(layers-1)*1000              # thickness of each layer (m) 
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
    
    return trm_layers, tau_layers


h_atm = 100.0  # thickness of atmosphere
layers = 80
dh = h_atm/(layers-1)
dx = 1
h, x = np.mgrid[0:100+dh:dh, 1:7+dx:dx]

tau_all = np.empty((layers, 0))
trm_all = np.empty((layers, 0))

# get the tau profile for each band and add to the 2D array
for filenum in list(['1','2','3','4','5','6','7']):

    filename = "cktable.lava_planet-B"+filenum+".nc"
    band_trm,band_tau = get_band_tau_profile(filename,h_atm,layers)
    tau_all = np.hstack((tau_all, band_tau[:, np.newaxis]))
    trm_all = np.hstack((trm_all, band_trm[:, np.newaxis]))

# plot the transmission profile
plt.figure(figsize=(10, 6))
heatmap = plt.pcolor(x,h,trm_all,cmap='viridis')
cb = plt.colorbar(heatmap)
cb.set_label('transmission',size=13)
plt.xlabel('wavenumber / $\mathrm{cm}^{-1}$',size=13)
plt.ylabel('height / km',size=13)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.savefig("/home/linfel/LavaPlanet/images/trm_profile.png",dpi=300)


# plot the optical depth profile
plt.figure(figsize=(10, 6))
heatmap = plt.pcolor(x,h,tau_all,cmap='YlOrRd')
cb = plt.colorbar(heatmap)
cb.set_label('optical depth',size=13)
plt.xlabel('wavenumber / $\mathrm{cm}^{-1}$',size=13)
plt.ylabel('height / km',size=13)
plt.tick_params(axis='both', which='major', labelsize=13)
plt.savefig("/home/linfel/LavaPlanet/images/tau_profile.png",dpi=300)

# flip the tau profile up side down, add a row of 0 at the first row, and write it in csv file
tau_topdown = np.flipud(tau_all)
row_of_zeros = np.zeros(tau_topdown.shape[1], dtype=tau_topdown.dtype)
tau_topdown = np.vstack((row_of_zeros,tau_topdown))
np.savetxt("/home/linfel/LavaPlanet/outputs/tau.csv", tau_topdown, delimiter=',')
