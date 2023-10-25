#! python3
from numpy import array, pi, genfromtxt
from pydisort import disort, get_legendre_coefficients, Radiant
from numpy.testing import assert_allclose
import os, unittest


# cdisort test01
class PyDisortTests(unittest.TestCase):
    def setUp(self):
        self.toml_path = "/home/linfel/pydisort/build/bin/isotropic_scattering.toml"
        self.tau_path = "/home/linfel/LavaPlanet/outputs/tau.csv"
        assert os.path.exists(self.toml_path), f"{self.toml_path} does not exist."
        assert os.path.exists(self.tau_path), f"{self.tau_path} does not exist."

    def test_isotropic_scattering(self):
        ds = disort.from_file(self.toml_path)
        ds.set_header("01. test isotropic scattering")

        # set dimension
        ds.set_atmosphere_dimension(
            nlyr=81, nstr=16, nmom=16, nphase=16
        ).set_intensity_dimension(nuphi=1, nutau=81, numu=6).finalize()

        # get scattering moments
        pmom = get_legendre_coefficients(ds.get_nmom(), "isotropic")

        # set boundary conditions
        ds.umu0 = 0.1
        ds.phi0 = 0.0
        ds.albedo = 0.0
        ds.fluor = 0.0

        # set output optical depth and polar angles
        umu = array([-1.0, -0.5, -0.1, 0.1, 0.5, 1.0])
        uphi = array([0.0])

        # case No.1
        print("==== Case No.1 ====")
        ds.fbeam = pi / ds.umu0
        ds.fisot = 0.0
        utau = genfromtxt(self.tau_path, delimiter=',', usecols = 0)
        ssa = array([0.2,0.3,0.1,0.2])

        tau = utau
        result = ds.run_with(
            {
                "tau": tau,
                "ssa": ssa,
                "pmom": pmom,
                "utau": utau,
                "umu": umu,
                "uphi": uphi,
            }
        ).get_intensity()
        print(result.shape)
        self.assertEqual(result.shape, (1, 81, 6))
        #print(result)


        flux = ds.get_flux()[:, [Radiant.RFLDIR, Radiant.FLDN, Radiant.FLUP]]
        print(flux.shape)
        print(flux)



if __name__ == "__main__":
    unittest.main()
