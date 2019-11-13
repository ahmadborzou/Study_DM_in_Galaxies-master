"""
Created on May 19 2019
author     : Ahmad Borzou, PhD
email      : ahmad_borzou[at]baylor.edu  
affiliation: Baylor University
"""

import Constants as co
########################################################################
## Set the initial values 
## m = mass
## n0 = number density at r=0
## T0 = temperature at r=0
## a, b are used to set the temperature profile as d^y/d(xi)^2 = -a xi^b
########################################################################

##############
## In SI units 
##############




################################
## Dwarfs with DM mass of 200 eV
#m  = co.m*200.
#n0 = 1.e-20/m 
#T0 = 0.0001
#n0 = 1.e-21/m 
#T0 = 0.000003
#a  = 0. 
#b  = 0.
##-------
#a  = 10.
#b  = 0.



##############################
## Dwarfs with DM mass of 2 eV
m  = co.m*2.
n0 = co.rho0/m*1.3
T0 = 0.001
a  = 2.2e6
b  = 0.
##-------
#a  = 0.
#b  = 0.


############################
## classic isothermal galaxy
#m  = co.m*1.e6
#n0 = 8.e-20/m #co.rho0/m 
#T0 = 1300.
#a  = 0.
#b  = 0.


#############
#m  = co.m*2.
#n0 = co.rho0/m*1.3
#T0 = 1.e-3
#a  = ?
#b  = ?
#y_pp = lambda xi, a, b: a*(b + xi)**(-2)


###################################
## second derivative of T/T0 wrt xi
## d^2(y)/d(xi)^2 
def y_pp(xi,a=a,b=b):
	#return <whatever other form of the temperature profile>
	return -a*xi**b
	

