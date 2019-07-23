# Python module for obtaining wavelength sensitivity for Raman intensity calibration

This module includes `wavelength_sensitivity.py` for performing fit using experimental band ratios to obtain the wavelength sensitivity modelled as a polynomial.

Usage
----------------
Required input data should be structured as following.

2D waves with columns  containing the data as follows.

|   J   |  Expt  ratio |  theoretical ratio |   J+2  transition wavenumber |  J-2  transition wavenumber | weight  |

For example,   (header  given for  clarity). Data wave as example is included in the present directory  as file, `dataD2.txt`

| J | Expt Ratio | Theortcl  Ratio | J+2, freq | J-2, freq | Weight |
|---|     ------------|-----------------      |-----------   |-----------   |--------    |
| 2 | 2.4872     | 2.3827               | 414.65    | -179.07   | 11.67   |
| 3 | 1.7910     | 1.6652               | 529.90    | -297.53   | 9.48    |
| 4 | 1.4674     | 1.3889               | 642.81    | -414.65   | 10.72  |
| 5 | 1.3904     | 1.2331               | 752.92    | -529.90   | 4.49    |
| 6 | 1.2955     | 1.1287               | 859.83    | -642.81   | 3.26    |


There are three functions for computing the residuals, `residual_linear`, `residual_quadratic` and `residual_cubic`. The  path to the data waves in the Igor experiment file is  to be modified in the function for the residual.


Requirements
----------------
Python 2.7 or Python 3.x with NumPy, SciPy and math modules

Usage
----------------
Following commands are run under the Python interpreter environment. (alternatively use any Python IDE like Spyder, IDLE or PyCharm)



1. After cloning the repository and moving in the `python-module` directory, and refer to the readme.  Prepare the required data as mentioned above and modify the path ....   (add the current folder to path allowing to import the module in the current folder; this is required for Python 2.7).
    > import sys

    > sys.path.append("..")

2. Import the module `wavelength_sensitivity` which should be in your current folder. (Directly execute the following command when using Python3)
    > import wavelength_sensitivity
