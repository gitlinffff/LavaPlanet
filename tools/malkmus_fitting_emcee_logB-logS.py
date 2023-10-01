from netCDF4 import Dataset
import math
import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt

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

# define the cumulative density distribution using Malkmus band model
def g_malkmus(k,B,S):
    a = 0.5 * (np.pi * B * S)**0.5
    b = 0.5 * (np.pi * B / S)**0.5
    # calculate the cumulative density distribution
    cdd = 0.5 * (1.0 - math.erf(a/k**0.5-b*k**0.5)) + 0.5 * (1.0 - math.erf(a/k**0.5+b*k**0.5)) * np.exp(np.pi * B)

    return cdd

# get the model prediction value of k
def approach_prediction_value(g_o,k_o,B,S):
    g_error_thres = 1e-4
    k_step = 1.0
    k_p = k_o
    g_1 = g_malkmus(k_p,B,S)
    error_before = g_1 - g_o

    while (abs(error_before) > g_error_thres) & (abs(k_step) > 1e-12):
        k_p = 10**(np.log10(k_p) + k_step)
        g_1 = g_malkmus(k_p,B,S)
        error_now = g_1 - g_o

        if error_before * error_now < 0:
            k_step = k_step * (-0.5)
        if (error_before * error_now > 0) & (abs(error_before) <= abs(error_now)):
            k_step = k_step * (-2)
        
        error_before = error_now

    return k_p

# calculate the likelihood function (sum of squares of errors) (SSE)
def ln_likelihood(params,g_obs,k_obs):
    B,S = params
    k_pred = []
    for i in range(len(g_obs)):
        g_o = g_obs[i]
        k_o = k_obs[i]
        k_p = approach_prediction_value(g_o,k_o,B,S)
        k_pred.append(k_p)
    assert len(k_pred) == len(k_obs), f"lengths of k_pred and k_obs not equal"
    SSE = np.sum((np.log10(k_pred) - np.log10(k_obs))**2)
    return -0.5 * SSE

# define the prior function
def ln_prior(log_params):
    logB,logS = log_params
    if (logB < -10) or (logB > 0) or (logS < -5) or (logS > 5):  # It is not good for logB to be larger than 0
        return -np.inf
    return 0.0

# define the Posterior Probability
def ln_posterior(log_params, x, y):
    lp = ln_prior(log_params)
    params = 10**log_params
    if not np.isfinite(lp):
        return -np.inf  
    return lp + ln_likelihood(params, x, y)


# specify input file
filename = "cktable.lava_planet-B3.nc"
pwd = "/home/linfel/LavaPlanet/"
cktable_fpath = pwd + filename
print(f"# working on {cktable_fpath}")

# read the data
ck_dataset = Dataset(cktable_fpath,'r')
ck_g = ck_dataset.variables['Wavenumber'][:]
ck_k_raw = ck_dataset.variables['SiO'][:,79,1]  # k at specific pressure and temperature

# adjust unit of k
molwt = 0.0440849 # molecular weight of SiO  0.0440849 kg/mol
ck_k = np.exp(ck_k_raw)/1000/molwt
log_ck_k = np.log10(ck_k)

# set the initial guess
ndim = 2  # Number of parameters (logB, logS)
nwalkers = 10  # Number of walkers (chains)
initial_logB_guess = -2
initial_logS_guess = -1
starting_guess = np.array([initial_logB_guess, initial_logS_guess])  # Provide your initial guesses
p0 = np.random.randn(nwalkers, ndim) * 0.1+ starting_guess       # np.random.randn:  used to generate random number in standard normal distribution

# Run the MCMC Sampler
sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_posterior, args=(ck_g, ck_k))
n_steps = 3000
sampler.run_mcmc(p0, n_steps, progress=True)

# final results
samples = sampler.chain[:, :, :].reshape((-1, ndim))
optm_params = np.median(samples, axis=0)
parameter_stddevs = np.std(samples, axis=0)
optm_SSE = -2*ln_likelihood(10**optm_params,ck_g, ck_k)
print(f"# optimal parameters [logB,logS]: {optm_params}")
print(f"# SSE at the optimal parameters: {optm_SSE}")
print(f"# parameters stddevs: {parameter_stddevs}")

# integrated autocorrelation time
try:
    tau = sampler.get_autocorr_time()
    print(tau)
except:
    print("# The chain is shorter than 50 times the integrated autocorrelation time for 2 parameter(s). Use this estimate with caution and run a longer chain!")
    
# get flat samples to make corner plot for the MCMC results
print('# creating the corner plot for MCMC ...')
flat_samples = sampler.get_chain(discard=200,thin=20,flat=True)
# discard: discard the first 200 steps in the chain
# thin: take only every thin steps from the chain.
# flat: flatten the chain across the ensemble
print(f"# (samples,dimension): {flat_samples.shape}")
fig_cornerplot = corner.corner(flat_samples, labels=["logB","logS"])
output_fig1 = "images/test3_corner.png"
fig_cornerplot.savefig(pwd+output_fig1,dpi=300)
print('# image saved  %s'%(pwd+output_fig1))

# get the optimal cumulative density distribution functions of Malkmus model
k_model = generate_uneven_seq(np.floor(log_ck_k.min())-1,np.ceil(log_ck_k.max())+1,3)
optm_B = 10**optm_params[0]
optm_S = 10**optm_params[1]
g_optm_model = [g_malkmus(ki,optm_B,optm_S) for ki in k_model]

# create the plot
print('# creating the model fit plot ...')
fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ax.plot(g_optm_model,k_model, linewidth=0.5)                           # plot optimal Malkmus model k-g curve
ax.plot(ck_g, ck_k, marker='o', markersize=1, linestyle='None')        # plot cktable points
ax.set_yscale('log')    # use a logarithmic scale for the x-axis
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel("g",size=15)
ax.set_ylabel("k(g)",size=15)
plt.title(f"SSE: {optm_SSE}\noptimal B: {optm_B:.4e}\noptimal S: {optm_S:.4e}",fontsize=8)
#ax.legend(fontsize = 15,markerscale=1.5)
output_fig2 = 'images/test3_model_fitting.png'
plt.savefig(pwd+output_fig2,dpi=300)
print('# image saved  %s'%(pwd+output_fig2))

print('# completed!')

