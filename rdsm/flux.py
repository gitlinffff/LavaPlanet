#! python3
import subprocess
from numpy import array, genfromtxt
import numpy as np
from pydisort import disort, get_legendre_coefficients, Radiant
import os, sys, unittest
sys.path.append('/home/linfel/LavaPlanet/build/python')
from pyathena import FileMode, io_wrapper, parameter_input
from pyharp import radiation, init_index_map
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument('bdnum', type=str, help='band #')
parser.add_argument('wdir', type=str, help='working directory')
parser.add_argument('bdinfo', type=str, help='spectral bands info  .inp file')
parser.add_argument('toml', type=str, help='toml file for running pydisort')
#parser.add_argument('opd', type=str, help='optical depth file')
parser.add_argument('atemp', type=str, help='the file keeping track of the atmospheric temperature')
#parser.add_argument('wmax', type=float, help='max wavenumber of the band')
#parser.add_argument('wmin', type=float, help='min wavenumber of the band')
args = parser.parse_args()


run_stir = "/home/linfel/LavaPlanet/rdsm/star_irradiance.py"
assert os.path.exists(args.toml), f"{args.toml} does not exist."
#assert os.path.exists(args.opd), f"{args.opd} does not exist."
assert os.path.exists(args.atemp), f"{args.atemp} does not exist."
assert os.path.exists(run_stir), f"{run_stir} does not exist."

ds = disort.from_file(args.toml)
ds.set_header("01. test isotropic scattering + thermal")

# get spectral bands info
os.chdir(args.wdir)
band_infile = io_wrapper()
pin = parameter_input()
band_infile.open(args.bdinfo, FileMode.read)
pin.load_from_file(band_infile)
band_infile.close()
init_index_map(pin)
rad = radiation()
rad.load_all_radiation_bands(pin)
nbands = rad.get_num_bands()
bands_info = []
for i in range(nbands):
    info = {}
    band = rad.get_band(i)
    info['name'] = band.get_name()
    info['wmin'] = band.get_wavenumber_min()
    info['wmax'] = band.get_wavenumber_max()
    bands_info.append(info)



# set level temperature
lyr_temp = genfromtxt(args.atemp, delimiter=',', usecols = -1)
nlayers = len(lyr_temp)
if nlayers == 1:
    temp = np.repeat(lyr_temp,2)
else:
    shifted_temp = np.roll(lyr_temp, shift=1)
    paired_temp = np.stack((shifted_temp, lyr_temp), axis=-1)
    ave_temp = np.mean(paired_temp, axis=1)[1:]
    temp = np.insert(ave_temp, 0, 1.5*lyr_temp[0]-0.5*lyr_temp[1])
    temp = np.append(temp, 1.5*lyr_temp[-1]-0.5*lyr_temp[-2])

# set dimension
ds.set_atmosphere_dimension(
    nlyr=nlayers, nstr=16, nmom=16, nphase=16
).set_intensity_dimension(nuphi=1, nutau=nlayers+1, numu=6).finalize()

# get scattering moments
pmom = get_legendre_coefficients(ds.get_nmom(), "isotropic")

# set boundary conditions
ds.umu0 = 1.0
ds.phi0 = 0.0
ds.albedo = 0.0
ds.fluor = 0.0
ds.fisot = 0.0
ds.btemp = temp[-1]
ds.ttemp = temp[0]
ds.temis = 0.0

# set output polar angles
umu = array([-1.0, -0.5, -0.1, 0.1, 0.5, 1.0])
uphi = array([0.0])

# case No.1
print("==== Case No.1 ====")

# set stellar parameters to calculate irradiance
T_star = 5172.0           # Temperature of 55 Cancri A
r_star = 0.943 * 6.957e8  #radius of the star
d_star = 2309800000       # distance to the star

# set single scattering albedo
ssa = array([0.0])

tot_net_flux = np.zeros(nlayers+1)    # create a zero array to record total net flux
# calculate net flux of each band
for bdnum in range(1,nbands+1):


    # get irradiance at the top of the planet atmosphere
    lmd_low = 10000.0/(bands_info[bdnum-1]['wmax'])
    lmd_up = 10000.0/(bands_info[bdnum-1]['wmin'])
    command = ['python', run_stir, str(T_star), str(r_star), str(d_star), str(lmd_low), str(lmd_up)]
    irradiance, err = subprocess.Popen(command,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE).communicate()
    # set beam flux (boundary condition)
    ds.fbeam = float(irradiance) / ds.umu0
    #ds.fbeam = 0.0
    # set wavenumber range for atmosphere emission
    ds.set_wavenumber_range_invcm(wmin=bands_info[bdnum-1]['wmin'], wmax=bands_info[bdnum-1]['wmax'])

    # set optical depth and output optical depth 
    tau_track_file = f"{args.wdir}/rdsm/outputs/tau_track_B{bdnum}.csv"
    tau = genfromtxt(tau_track_file, delimiter=',', usecols = -1)
    utau = np.insert(np.cumsum(tau),0,0.0)                  # utau is any optical depth where you want results

    result = ds.run_with(
        {
            "temp": temp,
            "tau": tau,
            "ssa": ssa,
            "pmom": pmom,
            "utau": utau,
            "umu": umu,
            "uphi": uphi,
        }
    ).get_intensity()
    #unittest.TestCase.assertEqual(result.shape, (1, nlayers+1, 6))
    assert result.shape == (1, nlayers+1, 6), f"result shape does not match"

    flux = ds.get_flux()[:, [Radiant.RFLDIR, Radiant.FLDN, Radiant.FLUP]]

    # get direct downward flux, diffusive downward flux, diffusive upward flux and net flux
    flux_dir = flux[:,0]
    flux_dn = flux[:,1]
    flux_up = flux[:,2]
    net_flux = flux_up - flux_dn - flux_dir
    #flux_dir = np.vstack(list0).T
    #flux_dn = np.vstack(list1).T
    #flux_up = np.vstack(list2).T

    # add up net flux from each band
    tot_net_flux = tot_net_flux + net_flux
    

#np.savetxt(f"/home/linfel/LavaPlanet/rdsm/outputs/flux_direct_B{args.bdnum}.csv", flux_dir, delimiter=',')
#np.savetxt(f"/home/linfel/LavaPlanet/rdsm/outputs/flux_down_B{args.bdnum}.csv", flux_dn, delimiter=',')
#np.savetxt(f"/home/linfel/LavaPlanet/rdsm/outputs/flux_up_B{args.bdnum}.csv", flux_up, delimiter=',')
#np.savetxt(f"/home/linfel/LavaPlanet/rdsm/outputs/flux_net_B{args.bdnum}.csv", net_flux, delimiter=',')
np.savetxt(f"/home/linfel/LavaPlanet/rdsm/outputs/total_net_flux.csv", tot_net_flux, delimiter=',')
