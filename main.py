"""
Created on May 19 2019
author     : Ahmad Borzou, PhD
email      : ahmad_borzou[at]baylor.edu  
affiliation: Baylor University
"""

from Stability import *


## solve the stability equations
lns_arr, lnsp_arr, lnspp_arr,\
y_arr, yp_arr, ypp_arr,\
r_arr, rho_arr, P_arr,\
M_arr = SolveStability()


## dr in SI units
dr = r_arr[1]-r_arr[0]


###################
## plotting section
###################	
print ("Working on the plots ...")
font = {
	'family': 'serif',
	'size': 14,
	'style': 'normal',
	'weight': 'medium',
	'fontname': 'Arial'
}
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')


## distances in kilo-parsec
r_arr_kpc = np.array(r_arr)/(co.parsec*1000)

## calculate the chemical potential
mu_arr = co.k*ins.T0*y_arr*(gal.lnz0+lns_arr)

## plot
fig, axs = plt.subplots(1,5, figsize=[16, 4])
ax = axs[0]
ax.plot(r_arr_kpc,rho_arr,linewidth=2,color='black',label='numeric')
ax.set_ylabel(r'$ \rho \quad\left( \mathrm{kg}\cdot\mathrm{m}^{-3} \right) $', fontdict=font)
ax.set_xlabel(r'$ \mathrm{r} \quad \left( \mathrm{kpc} \right) $', fontdict=font)
ax.set_yscale('log')
ax = axs[1]
ax.plot(r_arr_kpc,np.array(y_arr)*gal.T0,linewidth=2,color='black')
ax.set_ylabel(r'$ \mathrm{T} \quad \left( \mathrm{K} \right) $', fontdict=font)
ax.set_xlabel(r'$ \mathrm{r} \quad \left( \mathrm{kpc} \right) $', fontdict=font)
ax.set_xscale('log')
ax = axs[2]
lnz = gal.lnz0+lns_arr
ax.plot(r_arr_kpc,lnz,linewidth=2,color='black')
ax.set_ylabel(r'$ \mathrm{Ln}(\mathrm{z}) $', fontdict=font)
ax.set_xlabel(r'$ \mathrm{r} \quad \left( \mathrm{kpc} \right) $', fontdict=font)
ax = axs[3]
ax.plot(r_arr_kpc,np.array(M_arr)/co.MSun,linewidth=2,color='black',label='numeric')
ax.set_ylabel(r'$ \mathrm{M} \quad \left( \mathrm{M}_{\odot} \right)  $', fontdict=font)
ax.set_xlabel(r'$ \mathrm{r} \quad \left( \mathrm{kpc} \right) $', fontdict=font)
#ax.set_xscale('log')
ax.set_yscale('log')
ax = axs[4]
ax.plot(r_arr_kpc,mu_arr,linewidth=2,color='black')
ax.set_ylabel(r'$ \mu \quad\left( \mathrm{Joule} \right) $', fontdict=font)
ax.set_xlabel(r'$ \mathrm{r} \quad \left( \mathrm{kpc} \right) $', fontdict=font)
ax.set_xscale('log')
plt.subplots_adjust(hspace=0.4,wspace=0.4)
plt.tight_layout()
plt.show()












