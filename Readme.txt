# Raman_InstrResponseFunction
----------------------------

This is a set of functions in Python and IgorPro's scripting language for the determination of the Instrument Response Function (IRF) or the wavelength dependent sensitivity of the Raman spectrometer.    This repository requires the data on the Ji, frequency, intensity information for Stokes and anti-Stokes bands originating from a common rotational state in separate 2D data files/waves as the input (details later). Programs in Python and IgorPro are independent and perform the same job. The main scheme of the code is for the non-linear weighted minimization to obtain coefficients for a polynomial which represents the IRF. Independent validation of the obtained IRF should be done for a measure of accuracy.

## Input data required
----------------------------
 
 - Stokes_data (2D data matrix, with Ji, frequency, intensity as the columns)
 - antiStokes_data (2D data matrix, with Ji, frequency, intensity as the columns)
 

## Available programs
--------------------------- 



## Requirements for running program
---------------------------

 - Python 3 or 2.7 for the Pyhton module
 - IgorPro 6.36+ / 7+ for the Igor functions
 
 
## Usage
--------------------------- 
Clone the repository or download the zip file. As per your choice of the language ( Python or IgorPro) go through the specific README and proceed. 
 
## Comments
--------------------------- 
 - On convergence of the minimization
 
 - On accuracy
 

## Credits
--------------------------- 
 