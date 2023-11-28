import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('wdir', type=str, help='working directory')
parser.add_argument('atemp', type=str, help='the file keeping track of the atmospheric temperature')
parser.add_argument('nbands', type=int, help='number of bands')
args = parser.parse_args()


atm_temp = np.genfromtxt(args.atemp, delimiter=',', usecols = -1)
nlayers = len(atm_temp)


for bdnum in range(1,args.nbands+1):
    opd_input = f"{args.wdir}/rdsm/outputs/tau_B{bdnum}.csv"
    opd_grid = np.genfromtxt(opd_input, delimiter=',')
    # opd_grid should have one more row than atm_temp, because the first row of opd_grid is the temperature correponding to the optical depths
    assert len(opd_grid)==nlayers+1, f"layers in temperature file do not match with optical depths in opd input file."
    temp_grid = opd_grid[0,:]
    
    opd_profile = np.array([])
    for ilayer in range(1,nlayers+1):
        itemp = atm_temp[ilayer-1]        # temperature of the ith layer 
        iopd = np.interp(itemp, temp_grid, opd_grid[ilayer,:]) 
        opd_profile = np.append(opd_profile,iopd)

    np.savetxt(f"{args.wdir}/rdsm/outputs/tau_track_B{bdnum}.csv", opd_profile, delimiter=',')    
