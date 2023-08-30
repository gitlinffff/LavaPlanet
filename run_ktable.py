#! /usr/bin/env python3
import sys, os
sys.path.append('/home/linfel/LavaPlanet/build/python')

from pyathena import FileMode, io_wrapper, parameter_input
from pyharp import radiation, init_index_map
from multiprocessing import Pool, cpu_count
import os, re, subprocess, argparse

cmake_source_dir = '/home/linfel/canoe'
cmake_binary_dir = '/home/linfel/LavaPlanet/build'

# prepare parse and add argument
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
                    help = 'athena input file'
                    )
args = vars(parser.parse_args())

# file architecture
hitbin  = f"{cmake_binary_dir}/bin/hitbin.release"
rfm     = f"{cmake_binary_dir}/bin/rfm.release"
run_rfm = f"{cmake_binary_dir}/bin/run_rfm.py"
kcoeff  = f"{cmake_binary_dir}/bin/kcoeff.release"
cktable = f"{cmake_binary_dir}/bin/cktable.py"

parfile = ""
hitfile = "/home/linfel/LavaPlanet/co_HITRAN2020/lb_28Si-16O__SiOUVenIR__100-14285__296K.par"

base_name = os.path.basename(args['input'])
path_name = os.path.dirname(args['input'])
inpfile = base_name.split('.')[0]

# number of threads
max_threads = cpu_count() // 2

# ktable specifics
atm     = f"{cmake_source_dir}/data/lava_SiO_atm_isothermal.txt"

# temperature grid
temp    = "-50 50 3"

# flags
generate_tab = True
generate_nc  = True
generate_cktable = True

# spectral bands
file = io_wrapper()
pin = parameter_input()

if path_name != '':
    file.open(f"{path_name}/{inpfile}.inp", FileMode.read)
else:
    file.open(f"{inpfile}.inp", FileMode.read)
pin.load_from_file(file)
file.close()

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
    info['wres'] = band.get_wavenumber_res()
    cia_list, hitran_list = [], []
    for j in range(band.get_num_absorbers()):
        if band.get_absorber(j).get_category() == "cia":
            cia_list.append(band.get_absorber(j).get_name())
        if band.get_absorber(j).get_category() == "hitran":
            hitran_list.append(band.get_absorber(j).get_name())
    info['cia'] = '\"' + ' '.join(cia_list) + '\"'
    info['hitran'] = ' '.join(hitran_list)
    bands_info.append(info)

# number of parallel threads
nthreads = min(max_threads, len(bands_info))

# run ktable in single thread
def RunSingleKtable(info):
  bname = info['name']
  print(f"working on band {bname} ...")
  wmin, wmax, wres = info['wmin'], info['wmax'], info['wres']
  wave = f'{wmin} {wmax} {wres}'
  tab_folder = f'{wmin}-{wmax}'
  kinp = f'kcoeff.inp-{bname}'
  kncfile = f'kcoeff.{inpfile}-{bname}.nc'
  # create tab files and kcoeff.inp
  if generate_tab:
    script = ['python', run_rfm,'--hitbin', hitbin,'--rfm', rfm,
              '--par', parfile, '--hit', hitfile,
              '--atm',atm,'--wave',wave,'--temp',temp,
              '--molecule',info['hitran'],'--output',kinp,'--rundir',tab_folder]
    out, err = subprocess.Popen(script,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE).communicate()
    print(out.decode(), err.decode())
    #if not ('Successful' in out.decode()):
    #  raise RuntimeError("Error in generating tab files.")

  # run kcoeff
  if generate_nc:
    script = [kcoeff,'-i',kinp,'-o',kncfile]
    out, err = subprocess.Popen(script,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE).communicate()
    print(out.decode(), err.decode())
    #if err != b'':
    #    raise RuntimeError("Error in generating kcoeff.**.nc.")
  print("band %s finishes." % wave)

  # create correlated-K table
  if generate_cktable:
    script = [cktable,'--kcoeff',kncfile,'--atm',atm,'--cia',info['cia'],
              '--input', kinp, '--output', f'{inpfile}-{bname}',
              f' > log.cktable-{inpfile}-{bname} &']
    with open('run_cktable.sh', 'a') as f:
      f.write(' '.join(script))
      f.write('\n\n')
  print("band %s finishes." % wave)

if __name__ == '__main__':
    if generate_cktable:
        with open('run_cktable.sh', 'w') as f:
            f.write('#!/bin/bash\n')

    print(bands_info)
# parallel on spectral bands
    pool = Pool(nthreads)
    pool.map(RunSingleKtable, bands_info)
