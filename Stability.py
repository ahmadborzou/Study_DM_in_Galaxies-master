"""
Created on May 19 2019
author     : Ahmad Borzou, PhD
email      : ahmad_borzou[at]baylor.edu  
affiliation: Baylor University
"""

import numpy as np
from Galaxies import profile
import Constants as co
import Inputs as ins



## a galaxy(DM mass, num density, initial temperature)
gal = profile(ins.m,ins.n0,ins.T0)

def SolveStability():
	
	## temperature setup
	y    = 1. 
	y_p  = 0. 
	## fugacity setup
	lns = 0. ## s0=1 =>lns0 = 0
	lns_p = 0.
	## increments of xi
	dxi = gal.xi(co.parsec*0.1)
	## start xi from dxi rather than 0 to avoid singularity
	xi = dxi
	## number of loops in the while loop
	N = 0
	## mass of galaxy 
	M = 0.

		
	## containers to store the outputs
	xi_arr     = []
	r_arr      = []
	lns_arr    = []
	lnsp_arr   = []
	lnspp_arr  = []
	y_arr      = []
	yp_arr     = []
	ypp_arr    = []
	rho_arr    = []
	P_arr      = []
	M_arr      = []

	
	## append the initial values
	xi_arr.append(xi)
	r_arr.append(gal.r(xi))
	lns_arr.append(lns)
	lnsp_arr.append(lns_p)
	lnspp_arr.append(0.)
	y_arr.append(y)
	yp_arr.append(y_p)
	ypp_arr.append(0.)
	rho_arr.append(gal.m*gal.n(lns,y))
	P_arr.append(gal.P(lns,y))
	M_arr.append(M)	
	
	
	## stop the loop when density is 1/1000 of the initial value
	while gal.n(lns,y) > 0.001*gal.n0:
		## second derivative of y and Ln(s)
		y_pp   = ins.y_pp(xi)
		lns_pp = gal.lnspp(lns_p,lns,y_pp,y_p,y,xi)

		## determine Ln(s) and y 
		if N == 0:  ## Newton method in the first loop
			lns += lns_p*dxi
			y   += y_p*dxi
		else: ## Verlet method
			lns = -lns_arr[-2] + 2.*lns + lns_pp*dxi**2
			y   = -y_arr[-2]   + 2.*y   + y_pp*dxi**2
		## for unrealistic temperature profiles, temperature can turn negative
		## break the loop
		if y < 0.:
			break
		## determine first derivatives of Ln(s) and y
		if N == 0:  ## Newton method in the first loop
			lns_p += lns_pp*dxi
			y_p   += y_pp*dxi
		else: ## Verlet method
			lns_p = lnsp_arr[-2] + 2.*lns_pp*dxi
			y_p   = yp_arr[-2]   + 2.*y_pp*dxi
		## determine xi
		xi  += dxi   
		## mass at radius r in SI units
		M   += 4.*np.pi*(gal.r(xi))**2*gal.m*gal.n(lns,y)*gal.r(dxi)
		
		
		## break the loop and inform the user if fugacity or its derivatives are infinite
		if lns == np.inf:
			print "***\n***\nlog of z/z0 is infinite.\nBreaking the loop ...\n***\n***"
			break
		if lns_p == np.inf:
			print "***\n***\n1st derivative of log of z/z0 is infinite.\nBreaking the loop ...\n***\n***"
			break
		if lns_pp == np.inf:
			print "***\n***\n2nd derivative of log of z/z0 is infinite.\nBreaking the loop ...\n***\n***"
			break
			
		## append the quantities to the containers
		xi_arr.append(xi)
		r_arr.append(gal.r(xi))
		lns_arr.append(lns)
		lnsp_arr.append(lns_p)
		lnspp_arr.append(lns_pp)
		y_arr.append(y)
		yp_arr.append(y_p)
		ypp_arr.append(y_pp)
		rho_arr.append(gal.m*gal.n(lns,y))
		P_arr.append(gal.P(lns,y))
		M_arr.append(M)

		## print out
		try:
			if N%1000 == 0 and N>0:
				print "r: %1.1e (kpc) xi: %1.1e"%(gal.r(xi)/(co.parsec*1000.),xi)
				print "Ln(z): %1.1e Ln(s): %1.1e Ln(s)': %1.1e Ln(s)'': %1.1e"%(gal.lnz0+lns,lns,lns_p,lns_pp)
				print "T: %1.1e (K) y: %1.1e y': %1.1e y'': %1.1e"%(gal.T0*y,y,y_p,y_pp)
				print "n: %1.1e (1/m^3) n/n0: %1.1e"%(gal.n(lns,y),gal.n(lns,y)/gal.n0)
				print "M: %1.1e (sun)"%(M/co.MSun)
				print "--"
		except:
			## this is just a print. no action needed
			pass
			
		## change the counter
		N+=1
		
		


	###################################
	## convert the lists to numpy array
	lns_arr   = np.array(lns_arr)
	lnsp_arr  = np.array(lnsp_arr)
	r_arr     = np.array(r_arr)
	rho_arr   = np.array(rho_arr)
	P_arr     = np.array(P_arr)
	M_arr     = np.array(M_arr)	
	
	
	## return the containers
	return lns_arr, lnsp_arr, lnspp_arr,\
           y_arr, yp_arr, ypp_arr,\
		   r_arr, rho_arr, P_arr,\
		   M_arr
