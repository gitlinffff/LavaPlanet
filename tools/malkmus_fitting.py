# method: grid scanning, using error of log(k) to calculate SSE

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
def generate_uneven_seq(low_limit,up_limit,step_order):
    """ example:
        generate_uneven_seq(0,2,1):   [1, 1.1, 1.2, 1.3, 1.4,..., 9.9, 10, 11, 12, 13, 14,..., 99]
        generate_uneven_seq(-1,1,2):  [0.1, 0.101, 0.102, 0.103,..., 0.999, 1, 1.01, 1.02, 1.03,..., 9.99]
    """
    i = low_limit
    k_list = []
    while i<up_limit:
        list_i = list(np.arange(10**i, 10**(i+1), 10**(i-step_order)))
        k_list = k_list + list_i
        i += 1
    return k_list

# specify input file
filename = "cktable.lava_planet-B1.nc"
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

# build a series of cumulative density distribution functions of Malkmus model using different S and B parameters to determine the optimal B and S
k_model = generate_uneven_seq(np.floor(log_ck_k.min())-1,np.ceil(log_ck_k.max())+4,3)
B_list = generate_uneven_seq(-4,-3,0)#; B_list.reverse()
S_list = generate_uneven_seq(1,2,0)#; S_list.reverse()

# open a file to export (B,S,SSE)
outfile = open("/home/linfel/LavaPlanet/outputs/B1_BS_SSE.csv", 'w')
outfile.writelines(f"B,S,SSE\n")

# initialize objective variables
SSE_min = 10e20
B_optm = 0
S_optm = 0

# search best model that can minimize SSE
iteration = len(B_list)*len(S_list); ite = 0  # total number of iterations
for B in B_list:
    for S in S_list:
        ite += 1
        print('B %.10f  S %.10f, B_optm %.10f  S_optm %.10f  SSE_min %f,  progress %.2f%%'%(B,S,B_optm,S_optm,SSE_min,ite/iteration*100))
        
        try_next_parameter = False  # Initialize it as False. It might be set as True in the following codes in some conditions.
        g_model = [g_malkmus(ki,B,S) for ki in k_model]
        # check the model
        if ck_g.max() > max(g_model): print("ck_g.max() > max(g_model)")
        if ck_g.min() < min(g_model): print("ck_g.min() < min(g_model)")

        # calculate the sum of squares of errors (SSE)
        SSE = 0
        for i in range(len(ck_g)):
            g = ck_g[i]
            k_obsr = log_ck_k[i]
            # get the Malkmus k value based on the given g value 
            k_pred = ''
            for j in range(len(g_model)-1):
                if (g - g_model[j]) * (g - g_model[j+1]) <= 0:
                    k_pred = np.log10(k_model[j+1])
                    break
            if k_pred != '':
                SSE = SSE + (k_obsr - k_pred)**2
            else:
                try_next_parameter = True # Set it True when k_pred cannot be found due to inappropriate B or S value.
                break

        if try_next_parameter == False:
            if SSE < SSE_min:
                SSE_min = SSE; B_optm = B; S_optm = S
            
            # write the two parameters B and S as well as corresponding SSE to the file
            outfile.writelines(f"{B},{S},{SSE}\n")

outfile.close()
print('# %d iterations '%iteration)
print(f"# optimal B:{B_optm:.2e}    optimal S:{S_optm:.2e}   SSE_min:{SSE_min}")

# get the optimal cumulative density distribution functions of Malkmus model
g_optm_model = [g_malkmus(ki,B_optm,S_optm) for ki in k_model]

# create the plot
print('# creating the plot ...')
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.plot(g_optm_model,k_model, linewidth=0.5)                           # plot optimal Malkmus model k-g curve
ax.plot(ck_g, ck_k, marker='o', markersize=1, linestyle='None')        # plot cktable points

ax.set_yscale('log')    # use a logarithmic scale for the x-axis
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel("g",size=15)
ax.set_ylabel("k(g)",size=15)
plt.title(f"SSE:{SSE_min}\noptimal B:{B_optm:.2e}\noptimal S:{S_optm:.2e}",fontsize=8)
#ax.legend(fontsize = 15,markerscale=1.5)

output_fig = 'images/test1.png'
plt.savefig(pwd+output_fig,dpi=300)
print('# image saved  %s'%(pwd+output_fig))

print('# completed!')
