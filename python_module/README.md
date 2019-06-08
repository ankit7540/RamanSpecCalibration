# Python module for Intensity calibration 

This module includes two files, `wavelength_sensivity_fit.py`, `compute_spectra.py`. 

`wavelength_sensivity_fit.py` : It accepts the data of experimental rotational Raman intensity (from common rotational state) for and performinf a fit to determine the coefficients of a polynomial representing the wavelength dependent sensitivity.

`compute_spectra.py` : This is an independent module which imports energy state data (for H<sub>2</sub>, HD and D<sub>2</sub>) along with ratio of polarizability anisotropy matrix elements to compute the theoretical rotational Raman spectra of these gases. Temperature is included for Boltzmann population using which observable Raman spectra can be modeled. For obtaining the ratio of theoretical Raman intensities originating from a common rotational state (required for the `wavelength_sensivity_fit.py` module ), a function is included  here which simply takes the ratio of the intensities of the individual bands.


Requirements
----------------
Python 2.7 or Python 3.x with Numpy and math modules

Usage
----------------
Following commands are run under the Python interpreter environment. (alternatively use any Python IDE like Spyder, IDLE or PyCharm)

1. After cloning the repository and moving in the `python-module` directory, add the current folder to path allowing to import the module in the current folder (this is required for Python 2.7). 
    > import sys
    
    > sys.path.append("..")
     
2. Import the `rovibME` which should be in your current folder. (Directly execute the following command when using Python3)
    > import rovibME
    
    
