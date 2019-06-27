# Python module for Raman calibration

This module includes three files, `wavenumber_cal.py`,  `wavelength_sensivity_fit.py`, `compute_spectra.py`.

`wavenumber_cal.py` : It requires data on the band positions in pixels with the uncertainties, as vectors in a 1D txt file. (NumPy array can be directly used as well). Core program is just fitting the pixel with the reference wavelength data. The reference data in already included in the file with the error.

`wavelength_sensivity_fit.py` : It accepts the data of experimental rotational Raman intensity (from common rotational state) for and performing a fit to determine the coefficients of a polynomial representing the wavelength dependent sensitivity.

`compute_spectra.py` : This is an independent module which imports energy state data (for H<sub>2</sub>, HD and D<sub>2</sub>) along with ratio of polarizability anisotropy matrix (at 532.2 nm) elements to compute the theoretical rotational Raman spectra of these gases. Temperature is included for Boltzmann population using which observable Raman spectra can be modeled. For obtaining the ratio of theoretical Raman intensities originating from a common rotational state (required for the `wavelength_sensivity_fit.py` module ), a function is included  here which simply takes the ratio of the intensities of the individual bands. (For matrix elements of polarizability anisotropy at other wavelengths refer to this GitHub repository <https://github.com/ankit7540/H2-PolarizabilityMatrixElements> and research paper <https://aip.scitation.org/doi/abs/10.1063/1.5011433>.)


Requirements
----------------
Python 2.7 or Python 3.x with NumPy, SciPy and math modules

Usage
----------------
Following commands are run under the Python interpreter environment. (alternatively use any Python IDE like Spyder, IDLE or PyCharm)

1. After cloning the repository and moving in the `python-module` directory, and refer to the redme. (add the current folder to path allowing to import the module in the current folder; this is required for Python 2.7).
    > import sys

    > sys.path.append("..")

2. Import the module `rovibME` which should be in your current folder. (Directly execute the following command when using Python3)
    > import rovibME
