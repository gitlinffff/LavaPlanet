# method: gradient descend, using error of log(k) to calculate SSE

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

# calculate the sum of squares of errors (SSE)
def calculate_SSE(x_obs,y_obs,x_model,y_model):
    j_key = 0
    y_pred = []
    for i in range(len(x_obs)):
        x = x_obs[i]
        # get the Malkmus k value based on the given g value 
        for j in range(j_key,len(x_model)-1):
            if (x - x_model[j]) * (x - x_model[j+1]) <= 0:
                y_pred.append(y_model[j+1])
                j_key = j
                break
    #print(len(y_pred),len(y_obs))
    assert len(y_pred) == len(y_obs), f"lengths of y_pred and y_obs not equal" 
    SSE = np.sum((y_pred - y_obs)**2)
    
    return SSE

def get_n_neighbors(n,l,x0,y0):
    dtheta = 2*np.pi/n
    for i in range(n):
        x1 = x0 + l * np.cos(i*dtheta)
        y1 = y0 + l * np.sin(i*dtheta)
        if (x1<=0)|(y1<=0):
            continue
        yield (x1,y1) 

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
log_k_model = np.log10(k_model)

# open a file to export (B,S,SSE)
outfname = "outputs/B1_BS_SSE.csv"
outfile = open(pwd+outfname, 'w')
outfile.writelines(f"B,S,SSE\n")

# search best model that can minimize SSE
B_iterate = 6.0041e-03         # set the initial value of B when gradient descend starts 
S_iterate = 6.3089e1         # set the initial value of S when gradient descend starts
SSE_min = calculate_SSE(ck_g,log_ck_k,[g_malkmus(ki,B_iterate,S_iterate) for ki in k_model],k_model)
SSE_now = SSE_min       # SSE at the central (B,S) location in each iteration
num_nbr = 25            # set the number of directions for searching largest descend
step = 1e-1             # set the step length for each moving
diff = 1e10             # difference in SSE between two successive iterations
threshold = 1e-4        # stop iteration if the difference in SSE between two successive iterations falls below this threshold
num_loop = 0            # used for keeping recording how many iterations has completed

while diff > threshold:
    count = 0
    for B,S in get_n_neighbors(num_nbr,step,B_iterate,S_iterate):
        g_model = [g_malkmus(ki,B,S) for ki in k_model]
        # check the model
        #if num_loop == 250:
        #    #print(ck_g.max(),max(g_model),ck_g.min(),min(g_model))
        #    print(f"B: {B:.5e}  S: {S:.5e}")

        if (ck_g.max() > max(g_model))|(ck_g.min() < min(g_model)): 
            continue
        
        SSE_i = calculate_SSE(ck_g,log_ck_k,g_model,log_k_model)

        if SSE_i <= SSE_min:
            SSE_min = SSE_i; B_iterate = B; S_iterate = S
            count += 1
    
        num_loop += 1
    
    if count == 0:      # it means that SSE at all the sorrounding (B,S) locations are larger than SSE at the center, then decreasing the step is needed
        step = step/5
    else:               # it means that a lower SSE can be found among the sorrounding (B,S) locations
        diff = abs(SSE_min - SSE_now)
        SSE_now = SSE_min

    # write the two parameters B and S as well as corresponding SSE to the file
    outfile.writelines(f"{B_iterate:.6e},{S_iterate:.6e},{SSE_min}\n")
    print(f"B: {B_iterate:.4e}  S: {S_iterate:.4e}  SSE: {SSE_min}  loops: {num_loop}")

outfile.close()
print("# gradient descent track saved  %s"%(pwd+outfname))

# get the optimal cumulative density distribution functions of Malkmus model
g_optm_model = [g_malkmus(ki,B_iterate,S_iterate) for ki in k_model]

# create the plot
print('# creating the plot ...')
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.plot(g_optm_model,k_model, linewidth=0.5)                           # plot optimal Malkmus model k-g curve
ax.plot(ck_g, ck_k, marker='o', markersize=1, linestyle='None')        # plot cktable points

ax.set_yscale('log')    # use a logarithmic scale for the x-axis
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel("g",size=15)
ax.set_ylabel("k(g)",size=15)
plt.title(f"SSE:{SSE_min}\noptimal B:{B_iterate:.4e}\noptimal S:{S_iterate:.4e}",fontsize=8)
#ax.legend(fontsize = 15,markerscale=1.5)

output_fig = 'images/test1.png'
plt.savefig(pwd+output_fig,dpi=300)
print('# image saved  %s'%(pwd+output_fig))

print('# completed!')
