from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


band = {'1':'100-250',
        '2':'810-1340',
        '3':'1860-2500',
        '4':'2865-3704',
        '5':'3861-4873',
        '6':'4873-6050',
        '7':'6050-14285.7'}

for filenum in band.keys():
    filename = "cktable.lava_planet-B"+filenum+".nc"
    pwd = "/home/linfel/LavaPlanet/"
    cktable_fpath = pwd + filename

    # read the data
    ck_dataset = Dataset(cktable_fpath,'r')
    ck_g = ck_dataset.variables['Wavenumber'][:]
    ck_k_raw = ck_dataset.variables['SiO'][:,79,1]  # k at specific pressure and temperature

    # adjust unit of k
    molwt = 0.0440849 # molecular weight of SiO  0.0440849 kg/mol
    ck_k = np.exp(ck_k_raw)/1000/molwt


    # create the plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(ck_g, ck_k, marker='o', markersize=1, linestyle='None',label=f"{band[filenum]} "+"$\mathrm{cm}^{-1}$")        # plot cktable points
    ax.set_yscale('log')    # use a logarithmic scale for the x-axis
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel("g",size=15)
    ax.set_ylabel("k(g)",size=15)
    ax.legend(fontsize = 15,markerscale=1)
    output_fig = f'images/band{filenum}_ck.png'
    plt.savefig(pwd+output_fig,dpi=300)
    print('# image saved  %s'%(pwd+output_fig))
