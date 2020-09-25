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
M_arr, v02_vr2_arr, temphi_arr  = SolveStability()


## dr in SI units
dr = r_arr[1]-r_arr[0]

## Gravity Potential Energy
U_g  = -4.*np.pi*co.G*dr*sum(M_arr*rho_arr*r_arr)
## Kinetic energy d(U_k) = 3./2.*P*dV
U_k  = 3./2.*4.*np.pi*dr*sum(r_arr**2*P_arr)


## gravitational potential
phi0  = -4.*np.pi*co.G*dr*sum(rho_arr*r_arr)
phi_minus_phi0_arr = -co.G*M_arr/r_arr + 4.*np.pi*co.G*temphi_arr
phi_arr = phi_minus_phi0_arr + phi0
print ("Gravitational Potential Energy (U_g): %1.1e\nAverage Kinetic Energy (U_k): %1.1e"%(U_g,U_k))

## free falling velocity starting at rest from the edge
v02   = v02_vr2_arr[-1]
v_arr = np.sqrt(v02 - v02_vr2_arr)
v_    = v_arr[0:-1]## exclude v(R) which is zero
t_freefall = dr*sum(1/v_)
print ("dynamical time: %1.1e (yr)"%(t_freefall/(365.*24.*60.*60.)))


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
ax.plot(r_arr_kpc,rho_arr/co.rho_c,linewidth=2,color='black',label='numeric')
ax.set_ylabel(r'$ \rho \quad\left( \rho_{\mathrm{c}} \right) $', fontdict=font)
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












