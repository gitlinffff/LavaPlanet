from netCDF4 import Dataset
import os
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('spcs', type=str, help='atmospheric species')
parser.add_argument('bdnum', type=str, help='band #')
parser.add_argument('cktable', type=str, help='cktable file path')
parser.add_argument('hatm', type=float, help='initial thickness of atmosphere (km)')
args = parser.parse_args()


# read the data
ck_dataset = Dataset(args.cktable,'r')
ck_g = ck_dataset.variables['Wavenumber'][:]
Temperature = ck_dataset.variables['Temperature'][:]
rel_temp = ck_dataset.variables['T_SiO'][:]
Pressure = ck_dataset.variables['Pressure'][:]

nlayers = len(Pressure)
ntemp = len(rel_temp)

# correlated-K method
trm_profile = []                   # create a list to store transmission profile of the atmosphere
for ilayer in range(nlayers):
    T = Temperature[ilayer]                         # K
    p = Pressure[ilayer]                            # Pa
    rou = p/T/8.31446261815324                      # mol/m3
    dh = args.hatm/nlayers*1000                     # thickness of each layer (m) 
    u = rou * dh                                    # u is constant within each layer
    
    trm_lyr = np.array([])         # create an array to store transmissions for one specific layer under each temperature
    for itemp in range(ntemp):
        ck_k_raw = ck_dataset.variables[args.spcs][:,ilayer,itemp]  # k at specific pressure and temperature
        ck_k = np.exp(ck_k_raw)/1000                    # m2/mol

        # integrate over g space to calculate transmission
        trm = 0.0
        for j in range(len(ck_g)-1):
            k = (ck_k[j] + ck_k[j+1]) / 2
            delta_g = ck_g[j+1] - ck_g[j]
            trm = trm + np.exp(-k*u) * delta_g
        trm_lyr = np.append(trm_lyr,trm)
        
    trm_profile.append(trm_lyr)

trm_array = np.vstack(trm_profile)
tau_array = -np.log(trm_array)

# add temperatures as the first row of the tau and trm array
trm_array = np.vstack([rel_temp + Temperature[0], trm_array])
tau_array = np.vstack([rel_temp + Temperature[0], tau_array])

# save the tau and trm array
cwd = os.getcwd()
np.savetxt(f"{cwd}/outputs/trm_B{args.bdnum}.csv", trm_array, delimiter=',')
np.savetxt(f"{cwd}/outputs/tau_B{args.bdnum}.csv", tau_array, delimiter=',')
