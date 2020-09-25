"""
Created on May 19 2019
author     : Ahmad Borzou, PhD
email      : ahmad_borzou[at]baylor.edu  
affiliation: Baylor University
"""

## Fermi-Dirac integral\n,
## this function is useful in case ln(z)<100 which is still very degenerate case
## if ln(z)> 100, use the approximation by Sommerfeld

import numpy as np
from scipy.special import gamma
from scipy.integrate import quad
from scipy.optimize import root_scalar
import warnings
## the numerical calculation of the integral 
## returns warnings. But, we have compared the 
## results with numerical tables and validated
## them. Therefore, we suppress the warnings 
warnings.filterwarnings("ignore")

def f(nu,z):
	def integrand(x,z,nu):
		return x**(nu-1.)/(np.exp(x)/z+1.)
	ans   = quad(integrand,0,np.inf,args=(z,nu))/gamma(nu) 
	if abs(ans[0])>0. and ans[1]/ans[0] > 0.01:
		print ("very high error")
	return ans[0]
	
	
	
## pickle the partial degenerate values to speed up the code
#N = 0
#dz = 0.001
#z  = 0.01
#z_arr   = []
#f52_arr = []
#f32_arr = []
#f12_arr = []
#while z <= np.exp(10.):
#	if z > 0.1:
#		dz = 0.1
#	z_arr.append(z)
#	f52_arr.append(f(5./2.,z))
#	f32_arr.append(f(3./2.,z))
#	f12_arr.append(f(1./2.,z))
#	z += dz
#	if N%10000 == 0:
#		print("z: %1.2e"%(z))
#	N+=1
	
	
import pickle
## write a file
#with open("FermiDirac.pickle", "wb") as fil:
#	pickle.dump([np.array(z_arr),np.array(f52_arr),np.array(f32_arr),np.array(f12_arr)], fil)
#

#with open("FermiDirac.pickle", "rb") as fil:
#z_arr,f52_arr,f32_arr,f12_arr = pickle.load(fil)
	
	
