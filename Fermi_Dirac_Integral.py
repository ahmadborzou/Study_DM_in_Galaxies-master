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
	
	
	
