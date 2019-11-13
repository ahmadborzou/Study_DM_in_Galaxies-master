# Study DM in Galaxies
 
This is a software to study degenerate and non-degenerate spherical dark matter halos of galaxies with any temperature profile.

## How to run the software:

1. Set the initial values in Inputs.py
2. From the command line execute 
	$ python main.py
	
The software can also be used online at    
https://mybinder.org/v2/gh/ahmadborzou/Study_DM_in_Galaxies/master 

here is a short tutorial   
https://youtu.be/m-HqiTLKfOA





## Here is the architecture of the software:
The software is written in python 2.7 and depends on numpy 1.16.2, scipy 1.2.1, and matplotlib 2.2.4.

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
	The values of these integrals are computed and stored(pickled) in tables. In the run time the tables are used instead 
	of the functions of this file. This optimizes the speed of the program. When fugacity is less than 0.1 or larger than 
	exp(10), the approximations of the integrals are known and this exact forms will not be used. 
	
	5. Inputs.py:
	------------
	This is the file to be modified by User befor running the code. The user needs comment the default values and insert his/her
	DM mass, number density at the center, and temperature at the center. The inputs are all in the SI units. 
	Also \frac{d^y}{d\xi^2} is taken in this file. The default form is -a*\xi^b. If the user would like to work with this form, 
	the values of a and b should be given. If other forms of temperature profile are desired, the return of function y_pp() 
	should be overwritten
	
	6. Constants.py:
	---------------
	This file contains the SI values of 1. the mass of the sun 2. parsec 3. speed of light 4. a typical mass density 
	5. Planck constant 6. Boltzmann constant 7. Newton constant 8. the mass equivalent to 1 eV. 
	These values are often used in the software and can be used by the user when inserting inputs.


