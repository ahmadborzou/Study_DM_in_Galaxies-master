# Study DM in Galaxies
 
This is a software to study spherical dark matter halos of galaxies with any temperature profile at any degeneracy level.

## How to run the software:

1. Set the initial values in Inputs.py
2. From the command line execute 
	$ python main.py
	
The software can also be used online at    
https://mybinder.org/v2/gh/ahmadborzou/Study_DM_in_Galaxies-master.git/master

here is a short tutorial   
https://youtu.be/j4G_38baj_w





## Here is the architecture of the software:
The software is written in python 3.7 and depends on numpy 1.17.3, scipy 1.3.1, and matplotlib 3.1.1.

	1. main.py:
	----------
	Get the solutions of the stability equations and plot them. To run the code, one should call this file. 
	
	2. Stability.py:
	---------------
	Here the Verlet method is used to find the solutions of the stability equations. First a class of galaxy is constructed. 
	The second derivatives of the temperature and logarithm of fugacity are passed by the constructed galaxy. The numeric solutions
	are found by moving from the center in steps of 0.1 pc. To increase the accuracy of the solutions, decrease the value of 
	dxi to for example 0.01 or 0.001 pc. 
	Valid solutions are those that do not change if the accuracy is improved. 
	
	3. Galaxies.py:
	--------------
	Here the inputs from the user are used to construct a galactic profile. The parameters like \alpha, \xi, 
	Fermi-Dirac integrals, density, pressure, derivatives of the logarithm of fugacity using differential equations, 
	all corresponding to the User's inputs, are calculated here. 
	
	4. Fermi_Dirac_Integral.py:
	--------------------------
	The full definition of the Fermi-Dirac integrals are defined here. The integrals can be calculated numerically. 
	To find the degeneracy level z from the integrals, we solve an optimization problem. 
    When fugacity is less than 0.01 or larger than exp(30), the approximations of the integrals are known and this exact forms will not be used. 
	
	5. Inputs.py:
	------------
	This is the file to be modified by User before running the code. The user needs to comment the default values and insert his/her
	DM mass, number density at the center, and temperature at the center. The inputs are all in the SI units. 
	Also the temperature profile, and its first and second derivatives are taken in this file. The default form is y=(1+b\xi^n)^{-1} with n>1. If the user would like to work with this form, 
	the values of n and b should be set only. If other forms of temperature profile are desired, the return of functions 
	should be overwritten
	
	6. Constants.py:
	---------------
	This file contains the SI values of 1. the mass of the sun 2. parsec 3. speed of light 4. a typical mass density 
	5. Planck constant 6. Boltzmann constant 7. Newton constant 8. the mass equivalent to 1 eV. 
	These values are often used in the software and can be used by the user when inserting inputs.


