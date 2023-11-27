import subprocess
import sys
import os
sys.path.append('/home/linfel/LavaPlanet/build/python')
from pyathena import FileMode, io_wrapper, parameter_input
from pyharp import radiation, init_index_map

wdir = '/home/linfel/LavaPlanet'

# file architecture
run_opd     = f"{wdir}/rdsm/opd.py"
run_flx     = f"{wdir}/rdsm/flux.py"

bdnum = 1
bdinfo      = f"{wdir}/lava_planet.inp"
cktable     = f"{wdir}/2800isothermal/cktable.lava_planet-B{bdnum}.nc"
toml        =  "/home/linfel/pydisort/build/bin/sio_isotropic_thermal.toml"

spcs = 'SiO'

# specify the initial thickness of atmosphere (km)
hatm = 100.0

# get spectral bands info
os.chdir(wdir)
band_infile = io_wrapper()
pin = parameter_input()
band_infile.open(bdinfo, FileMode.read)
pin.load_from_file(band_infile)
band_infile.close()
init_index_map(pin)
rad = radiation()
rad.load_all_radiation_bands(pin)
bands_info = []
for i in range(rad.get_num_bands()):
    info = {}
    band = rad.get_band(i)
    info['name'] = band.get_name()
    info['wmin'] = band.get_wavenumber_min()
    info['wmax'] = band.get_wavenumber_max()
    bands_info.append(info)

# calculate optical depth for each layer and each temperature
command = ['python', run_opd, spcs, str(bdnum), cktable, str(hatm)]
out, err = subprocess.Popen(command,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE).communicate()

# calculate flux
tau = f"{wdir}/rdsm/outputs/tau_B{bdnum}.csv"
wmax = bands_info[bdnum-1]['wmax']
wmin = bands_info[bdnum-1]['wmin']

command = ['python', run_flx, toml, tau, str(wmax), str(wmin)]
out, err = subprocess.Popen(command,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE).communicate()
print(out.decode(), err.decode())
