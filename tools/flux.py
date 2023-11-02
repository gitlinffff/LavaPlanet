#! python3
import subprocess
from numpy import array, pi, genfromtxt
import numpy as np
from pydisort import disort, get_legendre_coefficients, Radiant
import os, unittest


# cdisort test01
class PyDisortTests(unittest.TestCase):
    def setUp(self):
        self.toml_path = "/home/linfel/pydisort/build/bin/sio_isotropic_thermal.toml"
        self.tau_path = "/home/linfel/LavaPlanet/outputs/tau.csv"
        self.atm_path = "/home/linfel/canoe/data/lava_SiO_atm_isothermal.txt"
        assert os.path.exists(self.toml_path), f"{self.toml_path} does not exist."
        assert os.path.exists(self.tau_path), f"{self.tau_path} does not exist."
        assert os.path.exists(self.atm_path), f"{self.atm_path} does not exist."

    def test_isotropic_scattering(self):
        ds = disort.from_file(self.toml_path)
        ds.set_header("01. test isotropic scattering")

        # wavelength of bands (um)
        bwv = { 0:[40.0,100.0],
                1:[7.462686567164179,12.34567901234568],
                2:[4.0,5.376344086021505],
                3:[2.699784017278618,3.490401396160558],
                4:[2.052123948286477,2.59000259000259],
                5:[1.652892561983471,2.052123948286477],
                6:[0.7,1.652892561983471]}
        
        # get number of layers
        nlayers = len(genfromtxt(self.tau_path, delimiter=',', usecols = 0))
        
        # set dimension
        ds.set_atmosphere_dimension(
            nlyr=nlayers, nstr=16, nmom=16, nphase=16
        ).set_intensity_dimension(nuphi=1, nutau=nlayers+1, numu=6).finalize()

        # set level temperature
        temp = genfromtxt(self.atm_path, delimiter='\t', usecols = 3, skip_header=1)
        temp = np.insert(temp,0,2800.0)
        
        # get scattering moments
        pmom = get_legendre_coefficients(ds.get_nmom(), "isotropic")

        # set boundary conditions
        ds.umu0 = 1.0
        ds.phi0 = 0.0
        ds.albedo = 0.0
        ds.fluor = 0.0
        ds.btemp = 2800.0
        ds.ttemp = 2800.0
        ds.temis = 0.1

        # set output optical depth and polar angles
        umu = array([-1.0, -0.5, -0.1, 0.1, 0.5, 1.0])
        uphi = array([0.0])

        # case No.1
        print("==== Case No.1 ====")
        
        # set stellar parameters to calculate irradiance
        T_star = 5172.0           # Temperature of 55 Cancri A
        r_star = 0.943 * 6.957e8  #radius of the star
        d_star = 2309800000       # distance to the star
        
        ds.fisot = 0.0
        ssa = array([0.0])
        
        list0 = []
        list1 = []
        list2 = []
        

        for band in range(1):
            lmd_low = bwv[band][0]
            lmd_up = bwv[band][1]
            command = ['python', 'star_irradiance.py', str(T_star), str(r_star), str(d_star), str(lmd_low), str(lmd_up)]
            irradiance, err = subprocess.Popen(command,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE).communicate()
            print(irradiance)
            ds.fbeam = float(irradiance) / ds.umu0
 
            ds.set_wavenumber_range_invcm(wmin=100.0, wmax=250.0)
            
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
