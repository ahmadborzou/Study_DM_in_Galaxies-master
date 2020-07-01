"""
Created on May 19 2019
author     : Ahmad Borzou, PhD
email      : ahmad_borzou[at]baylor.edu  
affiliation: Baylor University
"""

import Constants as co
import numpy as np



###########################
## Everythin in SI units ##
###########################




##################################################################################
## m: the mass of Dark matter   ## co.m = 1 eV in SI units
## n0: the number density of dark matter
## T0: the temperature at the center
## n and b: parameters of temperature profile in the form of y = (1 + b \xi^n)^{-1} with n > 1. ## set b=0 for isothermal 
##################################################################################
m  = co.m*100.
n0 = 1.e-22/m 
T0 = 0.0001
n = 2.	
b = 100.


##################################################################################
## the temperature profile and its first and second derivative with respect to \xi
##################################################################################
def y(xi):
	return 1/(1 + b*xi**n)
	
def y_p(xi):
	return -(b*n*xi**(-1 + n))/(1 + b*xi**n)**2

def y_pp(xi):
	return  (b*n*xi**(-2 + n)*(1 + b*xi**n + n*(-1 + b*xi**n)))/(1 + b*xi**n)**3

## make sure that the temperature profile has correct initial conditions
if abs(y(0)-1.)>1.e-5 or abs(y_p(0))>1.e-5:
	print("y0: %g yp0: %g"%(y(0),y_p(0)))
	raise(Exception("Wrong temperature profile. Initial values have to be y=0, y_p=0."))
##################################################################################
##################################################################################


