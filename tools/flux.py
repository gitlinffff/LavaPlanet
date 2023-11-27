#! python3
import subprocess
from numpy import array, genfromtxt
import numpy as np
from pydisort import disort, get_legendre_coefficients, Radiant
import os, unittest
import sys
sys.path.append('/home/linfel/LavaPlanet/build/python')
from pyathena import FileMode, io_wrapper, parameter_input
from pyharp import radiation, init_index_map


# cdisort test01
class PyDisortTests(unittest.TestCase):
    def setUp(self):
        self.toml_path = "/home/linfel/pydisort/build/bin/sio_isotropic_thermal.toml"
        self.tau_path = "/home/linfel/LavaPlanet/outputs/tau.csv"
        self.atm_path = "/home/linfel/canoe/data/lava_SiO_atm_isothermal.txt"
        self.stir_path = "/home/linfel/LavaPlanet/tools/star_irradiance.py"
        assert os.path.exists(self.toml_path), f"{self.toml_path} does not exist."
        assert os.path.exists(self.tau_path), f"{self.tau_path} does not exist."
        assert os.path.exists(self.atm_path), f"{self.atm_path} does not exist."
        assert os.path.exists(self.stir_path), f"{self.stir_path} does not exist."

    def test_isotropic_scattering(self):
        ds = disort.from_file(self.toml_path)
        ds.set_header("01. test isotropic scattering + thermal")

        # get spectral bands info
        os.chdir("/home/linfel/LavaPlanet")
        band_infile = io_wrapper()
        pin = parameter_input()
        band_infile.open("/home/linfel/LavaPlanet/lava_planet.inp", FileMode.read)
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
        
        # get number of layers
        nlayers = len(genfromtxt(self.tau_path, delimiter=',', usecols = 0))
        
        # set dimension
        ds.set_atmosphere_dimension(
            nlyr=nlayers, nstr=16, nmom=16, nphase=16
        ).set_intensity_dimension(nuphi=1, nutau=nlayers+1, numu=6).finalize()

        # set level temperature
        lyr_temp = genfromtxt(self.atm_path, delimiter='\t', usecols = 3, skip_header=1)
        if len(lyr_temp) == 1:
            temp = np.repeat(lyr_temp,2)
        else:
            shifted_temp = np.roll(lyr_temp, shift=1)
            paired_temp = np.stack((shifted_temp, lyr_temp), axis=-1)
            ave_temp = np.mean(paired_temp, axis=1)[1:]
            temp = np.insert(ave_temp, 0, 1.5*lyr_temp[0]-0.5*lyr_temp[1])
            temp = np.append(temp, 1.5*lyr_temp[-1]-0.5*lyr_temp[-2])
        
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
        
        list0 = []
        list1 = []
        list2 = []

        for band in range(rad.get_num_bands()):
            
            # get irradiance at the top of the planet atmosphere
            lmd_low = 10000.0/(bands_info[band]['wmax'])
            lmd_up = 10000.0/(bands_info[band]['wmin'])
            command = ['python', self.stir_path, str(T_star), str(r_star), str(d_star), str(lmd_low), str(lmd_up)]
            irradiance, err = subprocess.Popen(command,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE).communicate()
            
            # set beam flux (boundary condition)
            ds.fbeam = float(irradiance) / ds.umu0
            #ds.fbeam = 0.0
            # set wavenumber range for atmosphere emission
            ds.set_wavenumber_range_invcm(wmin=bands_info[band]['wmin'], wmax=bands_info[band]['wmax'])
           
            # set optical depth and output optical depth 
            tau = genfromtxt(self.tau_path, delimiter=',', usecols = band)
            utau = np.insert(np.cumsum(tau),0,0.0)  # utau is any optical depth where you want results
            
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
            self.assertEqual(result.shape, (1, nlayers+1, 6))

            flux = ds.get_flux()[:, [Radiant.RFLDIR, Radiant.FLDN, Radiant.FLUP]]
            
            list0.append(flux[:,0])
            list1.append(flux[:,1])
            list2.append(flux[:,2])

        flux_dir = np.vstack(list0).T
        flux_dn = np.vstack(list1).T
        flux_up = np.vstack(list2).T
        
        np.savetxt("/home/linfel/LavaPlanet/outputs/flux_direct.csv", flux_dir, delimiter=',')
        np.savetxt("/home/linfel/LavaPlanet/outputs/flux_down.csv", flux_dn, delimiter=',')
        np.savetxt("/home/linfel/LavaPlanet/outputs/flux_up.csv", flux_up, delimiter=',')

if __name__ == "__main__":
    unittest.main()
