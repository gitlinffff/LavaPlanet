from netCDF4 import Dataset
import math
import numpy as np
import matplotlib.pyplot as plt

# define the cumulative density distribution using Malkmus band model
def g_malkmus(k,B,S):
    a = 0.5 * (np.pi * B * S)**0.5
    b = 0.5 * (np.pi * B / S)**0.5
    # calculate the cumulative density distribution
    cdd = 0.5 * (1.0 - math.erf(a/k**0.5-b*k**0.5)) + 0.5 * (1.0 - math.erf(a/k**0.5+b*k**0.5)) * np.exp(np.pi * B)

    return cdd

# generate an uneven list of independent variable
def generate_k_list(low_limit,up_limit):
    i = low_limit
    k_list = []
    step_order = 3   # step relative to the range that it belongs to
    while i<up_limit:
        list_i = list(np.arange(10**i, 10**(i+1), 10**(i-step_order)))
        k_list = k_list + list_i
        i += 1
    return k_list

# specify input file
filename = "cktable.lava_planet-B5.nc"
pwd = "/home/linfel/LavaPlanet/"
cktable_fpath = pwd + filename

# read the data
ck_dataset = Dataset(cktable_fpath,'r')
ck_g = ck_dataset.variables['Wavenumber'][:]
ck_k_raw = ck_dataset.variables['SiO'][:,79,1]  # k at specific pressure and temperature

# adjust unit of k
molwt = 0.0440849 # molecular weight of SiO  0.0440849 kg/mol
ck_k = np.exp(ck_k_raw)/1000/molwt
log_ck_k = np.log10(ck_k)

# built the cumulative density distribution function of Malkmus model
k_model = generate_k_list(np.floor(log_ck_k.min()),np.ceil(log_ck_k.max()))
B = 1e-7
S = 1.0
g_model = [g_malkmus(ki,B,S) for ki in k_model]

# calculate the sum of squares of errors (SSE)
SSE = 0
for i in range(len(ck_g)):
    g = ck_g[i]
    k_obsr = ck_k[i]

    # get the Malkmus k value based on the given g value 
    for j in range(len(g_model)-1):
        if (g - g_model[j]) * (g - g_model[j+1]) <= 0:
            k_pred = k_model[j+1]
            break
        #if ((g - g_model[j]) > 0) & ((g - g_model[j+1]) > 0):
        #    print('g_model is invalid')

    SSE = SSE + (k_obsr - k_pred)**2

# create the plot
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.plot(g_model,k_model, linewidth=0.5)                           # plot Malkmus model k-g curve
ax.plot(ck_g, ck_k, marker='o', markersize=1, linestyle='None')   # plot cktable points

ax.set_yscale('log')    # use a logarithmic scale for the x-axis
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel("g",size=15)
ax.set_ylabel("k(g)",size=15)
plt.title('SSE %f'%SSE,fontsize=15)
#ax.legend(fontsize = 15,markerscale=1.5)

output_fig = 'images/test5.png'
plt.savefig(pwd+output_fig,dpi=300)
