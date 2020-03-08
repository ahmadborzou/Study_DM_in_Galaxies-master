"""
Created on May 19 2019
author     : Ahmad Borzou, PhD
email      : ahmad_borzou[at]baylor.edu  
affiliation: Baylor University
"""

import numpy as np
import Constants as co
from Fermi_Dirac_Integral import f as f_nu
from scipy.optimize import root_scalar


## for partial degeneracy use the pickled files
import cPickle as pickle
with open("FermiDirac.pickle", "rb") as fil:
	z_arr,f52_arr,f32_arr,f12_arr = pickle.load(fil)

class profile:
	## a constructor
	## inputs are dark matter mass, number density, and Temperature all at the center
	def __init__(self,m,n0,T0):
		## the mass of DM
		self.m    = m
		## the number density of DM at r = 0
		self.n0   = n0
		## the temperature of DM at r = 0
		self.T0   = T0
		## alpha parameter 
		self.alpha= co.h/np.sqrt(2.*np.pi*self.m)
		## Fermi_Dirac integral at xi=0 with nu=3/2
		## using the EOS
		f320 = self.n0*self.alpha**3/(2.*(co.k*self.T0)**(3./2.))
		## if z << 1, all Fermi-Dirac integrals are equal to z
		if f320 < 0.1:
			self.z0 = f320
			self.lnz0 = np.log(self.z0)
		else:
			## at this point there are 2 possibilities. 1. High degeneracy 2. partial degeneracy
			## assume that it is high degenerate. Use Sommerfeld approximation and EOS
			temp_lnz0 = (3.*np.sqrt(np.pi)*self.n0*self.alpha**3/(8.*(co.k*self.T0)**(3./2.)))**(2./3.)
			## test the assumption
			if temp_lnz0 > 10.:
				self.lnz0 = temp_lnz0
				self.z0   = np.exp(self.lnz0) 
			else:
				## if code arrives here z ~ 1 and system is partially degenerate
				## this is the slowest part because no perturbation exist to be used
				## z0 is the root of this function
				def RootFunc(z):
					return f_nu(3./2.,z)-f320
				## derivative of this function wrt z
				def RootFuncPrime(z):
					return f_nu(1./2.,z)/z
				## do the optimization to find z0
				ans = root_scalar(f=RootFunc,method='newton',fprime=RootFuncPrime,x0=1.)
				## make sure that it converged
				if ans.converged:
					self.z0   = ans.root
					self.lnz0 = np.log(self.z0)
				else:
					raise Exception("\n***\n\ncould not determine the value of z0. Exiting ..\n\n***\n")
		
	## a function to convert xi to r in SI units
	def r(self,xi):
		return xi*np.sqrt(self.alpha**3/(8.*np.pi*co.G*(self.m)**2*np.sqrt(co.k*self.T0)))

	## a function to convert r in SI units to dimensionless xi	
	def xi(self,r):
		return r*np.sqrt((8.*np.pi*co.G*(self.m)**2*np.sqrt(co.k*self.T0))/self.alpha**3)
		
	## Ln(z) = Ln(z0*s) = Ln(z0) + Ln(s)
	## this form prevents overflow
	def lnz(self,lns):
		ans = self.lnz0 + lns
		return ans

	## Dirac-Fermi function f_5/2
	def f_52(self,lns):
		## High degeneracy
		if self.lnz(lns) > 10:
			return 8./(15.*np.sqrt(np.pi))*(self.lnz(lns))**(5./2.)
		## classical
		elif self.z0*np.exp(lns) < 0.1:
			return self.z0*np.exp(lns)
		else:
			## partial degeneracy. Numerically solve the integral
			return f_nu(5./2.,self.z0*np.exp(lns))
			## use the pickled values. z increments are 0.1
			#return f52_arr[int(self.z0*np.exp(lns)/0.1)-1]
			
	## Dirac-Fermi function f_3/2
	def f_32(self,lns):
		if self.lnz(lns) > 10:
			return 4./(3.*np.sqrt(np.pi))*(self.lnz(lns))**(3./2.)
		elif self.z0*np.exp(lns) < 0.1:
			return self.z0*np.exp(lns)
		else:
			return f_nu(3./2.,self.z0*np.exp(lns))
			#return f32_arr[int(self.z0*np.exp(lns)/0.1)-1]
						
	## Dirac-Fermi function f_1/2
	def f_12(self,lns):
		if self.lnz(lns) > 10:
			return 2./np.sqrt(np.pi)*(self.lnz(lns))**(1./2.)
		elif self.z0*np.exp(lns) < 0.1:
			return self.z0*np.exp(lns)
		else:
			return f_nu(1./2.,self.z0*np.exp(lns))		
			#return f12_arr[int(self.z0*np.exp(lns)/0.1)-1]
		
	## the second derivative of Ln(s) (s = z/z0)
	def lnspp(self,lnsp,lns,ypp,yp,y,xi):
		## high degeneracy
		if self.lnz(lns) > 10:
			## this is f52/f32
			h  = self.f_52(lns)/self.f_32(lns)
			hp = lnsp*(1. - self.f_52(lns)*self.f_12(lns)/(self.f_32(lns))**2)			
		## classic case
		elif self.z0*np.exp(lns) < 0.1:
			## this is f52/f32
			h  = 1.
			## its derivative
			hp = 0.
		## partial degeneracy
		else:
			## this is f52/f32
			h  = self.f_52(lns)/self.f_32(lns)
			hp = lnsp*(1. - self.f_52(lns)*self.f_12(lns)/(self.f_32(lns))**2)
		## at xi ~ 0 sp, and yp are closer to 0 than xi.
		if self.r(xi) < co.parsec*0.0001:
			return -np.sqrt(y)*self.f_32(lns) - yp*lnsp/y             - 5.*h*ypp/(2.*y) - 5.*hp*yp/(2.*y) 
		else:
			## write sp^2/s as sp/s*sp to avoid overflow
			return -np.sqrt(y)*self.f_32(lns) - yp*lnsp/y -2.*lnsp/xi - 5.*h*ypp/(2.*y) - 5.*hp*yp/(2.*y) -5.*h*yp/(xi*y)

		
	## derivative of number density n wrt xi   
	def n_p(self,lnsp,lns,yp,y):
		return 2.*(co.k*self.T0)**(3./2.)/self.alpha**3*(  \
				3./2.*np.sqrt(y)*yp*self.f_32(lns)+  \
				y**(3./2.)*self.f_12(lns)*lnsp  \
				)
				
	## number density n 
	def n(self,lns,y):
		return 2.*(co.k*self.T0*y)**(3./2.)/self.alpha**3*self.f_32(lns)
	
	
	## pressure
	def P(self,lns,y):
		return 2.*(co.k*self.T0*y)**(5./2.)/self.alpha**3*self.f_52(lns)
		
