from scipy.integrate import quad
from numpy import exp,pi
import argparse

def PlanckFunction(lmd,T):

    h = 6.62607015e-34   # Js
    c = 299792458        # m/s
    k = 1.380649e-23     # J/K

    B = (2*h*c*c)/(lmd**5)/(exp(h*c/k/lmd/T)-1)  # W m-2 sr-1 / m

    return B

parser = argparse.ArgumentParser()
parser.add_argument('T', type=float, help='Temperature of the star (K)')
parser.add_argument('lmd_low', type=float, help='wavelength lower limit of integral of Planck Function (um)')
parser.add_argument('lmd_up', type=float, help='wavelength upper limit of integral of Planck Function (um)')
args = parser.parse_args()

I = quad(PlanckFunction, args.lmd_low*1e-6, args.lmd_up*1e-6, args=args.T)
print(I[0])
