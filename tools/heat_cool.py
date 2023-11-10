import numpy as np
from netCDF4 import Dataset

list_hrate = []    
for band in range(7): 
    # read flux data
    flux_dir = np.genfromtxt("/home/linfel/LavaPlanet/outputs/flux_direct.csv", delimiter=',',usecols = band)
    flux_down = np.genfromtxt("/home/linfel/LavaPlanet/outputs/flux_down.csv", delimiter=',',usecols = band)
    flux_up = np.genfromtxt("/home/linfel/LavaPlanet/outputs/flux_up.csv", delimiter=',',usecols = band)

    # read temperature and pressure profile from cktable file
    ck_fname = f"cktable.lava_planet-B{str(band+1)}.nc"
    pwd = "/home/linfel/LavaPlanet/"
    ck_fpath = pwd + ck_fname
    ck_dataset = Dataset(ck_fpath,'r')
    Temperature = ck_dataset.variables['Temperature'][:]
    Pressure = ck_dataset.variables['Pressure'][:]

    # read heat capacity data (J/mol/K)
    cptable = np.genfromtxt("/home/linfel/LavaPlanet/ExoMol_to_HITRAN/28Si-16O/28Si-16O__SiOUVenIR.cp", delimiter='')
    tem_cp = cptable[:,0]
    cp_value = cptable[:,1]

    # check the number of levels
    assert len(flux_down) == len(flux_up), f"Number of levels in downward flux is not equal to number of levels in upward flux."
    assert len(flux_down) == len(Pressure)+1, f"Number of levels in flux is not equal to number of layers in pressure plus 1."

    # number of levels
    nlevel = len(flux_down)

    # calculate net flux
    net_flux = flux_up - flux_down - flux_dir

    # define parameters
    z_atm = 100.0              # height of atmosphere (km)
    dz = z_atm/(nlevel-1)      # height of layer (km)
    molwt = 0.0440849          # molecular weight of SiO  0.0440849 (kg/mol)
    R = 8.31446261815324       # gas constant (J/mol/K)

    heat_rate = np.array([], dtype=float)
    for layer in range(nlevel-1):
        P = Pressure[layer]
        T = Temperature[layer]
        # calculate density
        rou = P*molwt/R/T
        # locate the heat capacity at the specific temperature
        cpidx = np.argmin(np.abs(tem_cp - T))
        cp = cp_value[cpidx] / molwt
        # calculate heating rate
        dTdt = (net_flux[layer] - net_flux[layer+1]) / dz / 1000 / (-rou * cp)

        heat_rate = np.append(heat_rate, dTdt)
    
    list_hrate.append(heat_rate)

heat_rate_profile = np.vstack(list_hrate).T

total_rate = np.sum(heat_rate_profile, axis=1)
heat_rate_profile = np.column_stack((heat_rate_profile, total_rate))
np.savetxt("/home/linfel/LavaPlanet/outputs/heat_rate.csv", heat_rate_profile, delimiter=',')
