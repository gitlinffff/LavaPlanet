from netCDF4 import Dataset
import numpy as np

# specify input file
filenum = '1'
filename = "cktable.lava_planet-B"+filenum+".nc"
pwd = "/home/linfel/LavaPlanet/"
cktable_fpath = pwd + filename
print(f"# working on {cktable_fpath}")

# read the data
ck_dataset = Dataset(cktable_fpath,'r')
ck_g = ck_dataset.variables['Wavenumber'][:]
Temperature = ck_dataset.variables['Temperature'][:]
Pressure = ck_dataset.variables['Pressure'][:]

layers = 80                       # number of layers
tau_layers = np.array([])         # create an empty array to store tau of each layer

for ilayer in range(layers):
    ck_k_raw = ck_dataset.variables['SiO'][:,ilayer,1]  # k at specific pressure and temperature
    ck_k = np.exp(ck_k_raw)/1000  # m2/mol
    T = Temperature[ilayer]       # K
    p = Pressure[ilayer]          # Pa

    rou = p/T/8.31446261815324    # mol/m3
    d = 10.0/79*1000              # thickness of each layer (m) 
    u = rou * d                   # u is constant within each layer

    # integrate over g space to calculate transmission
    transmission = 0.0
    for j in range(len(ck_g)-1):
        k = (ck_k[j] + ck_k[j+1]) / 2
        delta_g = ck_g[j+1] - ck_g[j]
        transmission = transmission + np.exp(-k*u) * delta_g

    # calculate tau for this layer
    tau = -np.log(transmission)
    tau_layers = np.append(tau_layers,tau)

print(len(tau_layers))
print(tau_layers)
