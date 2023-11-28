import numpy as np
from netCDF4 import Dataset
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument('bdnum', type=str, help='band #')
parser.add_argument('nfx', type=str, help='net flux file for each level of the atmosphere')
parser.add_argument('cktable', type=str, help='cktable file path of any spectral band')
parser.add_argument('cptable', type=str, help='heat capacity file path')
parser.add_argument('atemp', type=str, help='the file keeping track of the atmospheric temperature')
parser.add_argument('hatm', type=float, help='initial thickness of atmosphere (km)')
parser.add_argument('molwt', type=float, help='molecular weight of the species (kg/mol)')
args = parser.parse_args()


# read net flux data
net_flux = np.genfromtxt(args.nfx, delimiter=',')
#flux_down = np.genfromtxt("/home/linfel/LavaPlanet/outputs/flux_down.csv", delimiter=',',usecols = band)
#flux_up = np.genfromtxt("/home/linfel/LavaPlanet/outputs/flux_up.csv", delimiter=',',usecols = band)

# read pressure profile from initial atmospheric structure file
#Pressure = genfromtxt(args.atm, delimiter='\t', usecols = 3, skip_header=1)
ck_dataset = Dataset(args.cktable,'r')
Pressure = ck_dataset.variables['Pressure'][:]

# read temperature profile
Temperature = np.genfromtxt(args.atemp, delimiter=',', usecols = -1)

# read heat capacity data (J/mol/K)
tem_cp = np.genfromtxt(args.cptable, delimiter='', usecols = 0)
cp_value = np.genfromtxt(args.cptable, delimiter='', usecols = 1)

# check the number of levels
assert len(net_flux) == len(Pressure)+1, f"Number of levels in flux is not equal to number of layers in pressure plus 1."

# number of levels
nlevel = len(net_flux)

# define parameters
dh = args.hatm/(nlevel-1)      # height of layer (km)
R = 8.31446261815324       # gas constant (J/mol/K)

heat_rate = np.array([], dtype=float)
for layer in range(nlevel-1):
    P = Pressure[layer]
    T = Temperature[layer]
    # calculate density
    rou = P * args.molwt / R / T
    # locate the heat capacity at the specific temperature
    cpidx = np.argmin(np.abs(tem_cp - T))
    cp = cp_value[cpidx] / args.molwt
    # calculate heating rate
    dTdt = (net_flux[layer] - net_flux[layer+1]) / dh / 1000 / (-rou * cp)

    heat_rate = np.append(heat_rate, dTdt)
    
#    list_hrate.append(heat_rate)

#heat_rate_profile = np.vstack(list_hrate).T

#total_rate = np.sum(heat_rate_profile, axis=1)
#heat_rate_profile = np.column_stack((heat_rate_profile, total_rate))
#np.savetxt(f"/home/linfel/LavaPlanet/rdsm/outputs/heat_rate_B{args.bdnum}.csv", heat_rate, delimiter=',')
np.savetxt(f"/home/linfel/LavaPlanet/rdsm/outputs/heat_rate.csv", heat_rate, delimiter=',')
