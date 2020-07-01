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
import warnings



## a galaxy(DM mass,n0,T0)
gal = profile(ins.m,ins.n0,ins.T0)

def SolveStability():

	## increments of xi
	dxi = gal.xi(co.parsec*0.1)
	## start xi from dxi rather than 0 to avoid singularity
	xi = np.copy(dxi)
	
	## temperature setup
	y    = ins.y(xi) 
	y_p  = ins.y_p(xi)

	## fugacity setup
	lns = 0. ## s0=1 =>lns0 = 0
	lns_p = 0.
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
	rhop_arr   = []
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
	rhop_arr.append(gal.m*gal.n_p(lns_p,lns,y_p,y))
	P_arr.append(gal.P(lns,y))
	M_arr.append(M)	
	
	
	## stop the loop when density is 1/1000 of the initial value
	while gal.n(lns,y) > 1.e-3*gal.n0:
		#if gal.r(xi)/(1000.*co.parsec) < 0.04:
		#	dxi = gal.xi(co.parsec*0.001)
		## second derivative of y and s
		y_pp   = ins.y_pp(xi)
		lns_pp = gal.lnspp(lns_p,lns,y_pp,y_p,y,xi)
		
		## set the temperature and its derivative
		y    = ins.y(xi) 
		y_p  = ins.y_p(xi)
		
		## determine s and y 
		if N == 0:  ## Newton method in the first loop
			lns += lns_p*dxi
		else: ## Verlet method
			lns = -lns_arr[-2] + 2.*lns + lns_pp*dxi**2
		## for unrealistic temperature profiles temperature can turn negative
		## break the loop and inform the user
		if y < 0.:
			print ("***\n***\nNegative temperature. Breaking ...\n***\n***")
			raise(Exception("Negative temperature. Unacceptable solution"))
			break
		## determine first derivative of s
		if N == 0:  ## Newton method in the first loop
			lns_p += lns_pp*dxi
		else: ## Verlet method
			lns_p = lnsp_arr[-2] + 2.*lns_pp*dxi
		## determine xi
		xi  += dxi   
		## mass at radius r in SI units using M = M + dM
		M   += 4.*np.pi*(gal.r(xi))**2*gal.m*gal.n(lns,y)*gal.r(dxi)
				
		## break the loop and inform the user if derivative of density is infinite  
		if gal.n_p(lns_p,lns,y_p,y) == np.inf:
			print ("***\n***\nderivative of number density is infinite.\nnumber density will sharply fall to zero.\nthis is usually the edge of the system.\nBreaking the while loop ...\n***\n***")
			break
		## break the loop and inform the user if fugacity or its derivatives are infinite
		if lns == np.inf:
			print ("***\n***\nlog of z/z0 is infinite.\nThis is full degeneracy limit. We can't handle this at this point.\nBreaking the loop ...\n***\n***")
			break
		if lns_p == np.inf:
			print ("***\n***\n1st derivative of log of z/z0 is infinite.\nWe can't handle this at this point.\nBreaking the loop ...\n***\n***")
			break
		if lns_pp == np.inf:
			print ("***\n***\n2nd derivative of log of z/z0 is infinite.\nWe can't handle this at this point.\nBreaking the loop ...\n***\n***")
			break
		## add the calculated quantities into the containers
		xi_arr.append(xi)
		r_arr.append(gal.r(xi))
		lns_arr.append(lns)
		lnsp_arr.append(lns_p)
		lnspp_arr.append(lns_pp)
		y_arr.append(y)
		yp_arr.append(y_p)
		ypp_arr.append(y_pp)
		rho_arr.append(gal.m*gal.n(lns,y))
		rhop_arr.append(gal.m*gal.n_p(lns_p,lns,y_p,y))
		P_arr.append(gal.P(lns,y))
		M_arr.append(M)

		try:
			if N%1000 == 0 and N>0:
				print ("r: %1.1e (kpc) xi: %1.1e "%(gal.r(xi)/(co.parsec*1000.),xi))
				print ("lnz: %1.1e lns: %1.1e lnsp: %1.1e lnspp: %1.1e"%(gal.lnz0+lns,lns,lns_p,lns_pp))
				print ("T: %1.1e (Kelvin) y: %1.1e yp: %1.1e ypp: %1.1e"%(gal.T0*y,y,y_p,y_pp))
				print ("n: %1.1e (1/m^3) rho: %1.1e (kg/m^3)"%(gal.n(lns,y),gal.n(lns,y)*gal.m))
				print ("M: %1.1e (sun)"%(M/co.MSun))
				print ("--")
		except:
			## this is just a print. no action needed
			pass
			
		## change the counter
		N+=1
		
		


	###################################
	## convert the lists to numpy array
	y_arr   = np.array(y_arr)
	yp_arr  = np.array(yp_arr)
	ypp_arr = np.array(ypp_arr)
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
