import subprocess
from netCDF4 import Dataset
import numpy as np

# working directory
wdir = '/home/linfel/LavaPlanet'

# file architecture
run_opd     = f"{wdir}/rdsm/opd.py"
run_pick_opd = f"{wdir}/rdsm/pkopd.py"
run_flx     = f"{wdir}/rdsm/flux.py"
run_hc      = f"{wdir}/rdsm/heat_cool.py"

bdinfo      = f"{wdir}/lava_planet.inp"
toml        =  "/home/linfel/pydisort/build/bin/sio_isotropic_thermal.toml"
cptable     = f"{wdir}/ExoMol_to_HITRAN/28Si-16O/28Si-16O__SiOUVenIR.cp"
atm         =  "/home/linfel/canoe/data/lava_SiO_atm_isothermal.txt"

# number of bands
nbands = 8

# atmospheric species
spcs = 'SiO'
molwt = 0.0440849   # molecular weight of SiO  0.0440849 (kg/mol)

# specify the initial thickness of atmosphere (km)
hatm = 100.0

# iteration parameters for atmospheric heating and cooling
tstep = 1.0    # second
iterations = 120

# print running headings
print(f"# number of bands: {nbands}")
print(f"# atmospheric species: {spcs}")
print(f"# initial atmospheric thickness: {hatm:.1f} km")
print(f"# time step: {tstep:.1f} s")
print(f"# iterations: {iterations}")

# construct the array of optical depths for each layer and each temperature
for bdnum in range(1,nbands+1):
    print(f"constructing array of optical depths for band {bdnum} ...")
    cktable = f"{wdir}/2800isothermal/cktable.lava_planet-B{bdnum}.nc"
    command = ['python', run_opd, spcs, str(bdnum), cktable, str(hatm)]
    out, err = subprocess.Popen(command,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE).communicate()

# record the initial atmospheric temperature profile
ck_dataset = Dataset(cktable,'r')
temp_trk = ck_dataset.variables['Temperature'][:]
np.savetxt(f"{wdir}/rdsm/outputs/atmtemp.csv", temp_trk, delimiter=',')
#temp_trk = np.genfromtxt(f"{wdir}/rdsm/outputs/temp_trk.csv", delimiter=',').T


for i in range(iterations):

    # obtain optical depths at a given temperature profile of the atmosphere
    atemp = f"{wdir}/rdsm/outputs/atmtemp.csv"
    command = ['python', run_pick_opd, wdir, atemp, str(nbands)]
    out, err = subprocess.Popen(command,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE).communicate()
    print(out.decode(), err.decode())

    # calculate flux
    atemp = f"{wdir}/rdsm/outputs/atmtemp.csv"
    command = ['python', run_flx, wdir, str(nbands), bdinfo, toml, atemp]
    out, err = subprocess.Popen(command,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE).communicate()
    print(out.decode(), err.decode())

    # calculate heating rate
    #nfx = f"{wdir}/rdsm/outputs/flux_net_B{bdnum}.csv"
    nfx = f"{wdir}/rdsm/outputs/total_net_flux.csv"
    atemp = f"{wdir}/rdsm/outputs/atmtemp.csv"
    command = ['python', run_hc, nfx, cktable, cptable, atemp, str(hatm), str(molwt)]
    out, err = subprocess.Popen(command,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE).communicate()

    # obtain new atmospheric temperature after one step of evolution
    hcrate = np.genfromtxt(f"{wdir}/rdsm/outputs/heat_rate.csv", delimiter=',', usecols = -1)
    cur_temp = np.genfromtxt(atemp, delimiter=',', usecols = -1)
    nxt_temp = cur_temp + tstep * hcrate
    np.savetxt(atemp, nxt_temp, delimiter=',')
    assert not np.any(nxt_temp < 0), "invalid negative temperature occurred." 

    temp_trk = np.vstack((temp_trk, nxt_temp)) 


    # percentage progress
    progress = (i+1) / iterations * 100
    print(f"Progress: {progress:.2f}%   ", end="\r",flush=True)


temp_trk = temp_trk.T

np.savetxt(f"{wdir}/rdsm/outputs/temp_trk.csv", temp_trk, delimiter=',')
