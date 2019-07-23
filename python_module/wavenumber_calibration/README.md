# Python module for wavenumber calibration

This module includes `wavenum_calbr.py` for performing fit of  reference transition wavenumbers against the band positions in pixels.  

Usage
----------------
Four 1D numpy arrays are required.
  - reference  transition wavenumbers
  - error in  the reference (sigma)

  - band positions in pixels
  - error in band  positions  

In the present exmaple, these are included as `freq.txt`, `freq_sigma.txt`, `pixel.txt` and `pixel_sigma.txt`  



Requirements
----------------
Python 2.7 or Python 3.x with NumPy  and SciPy modules

Usage
----------------
Following commands are run under the Python interpreter environment. (Alternatively use any Python IDE like Spyder, IDLE or PyCharm). *Spyder3 is has been used while writing and debugging.*

***When using Python interpreter in terminal***

1. After cloning the repository and moving in the `python-module` directory,  refer to the readme.  Prepare the required data as mentioned above which will be loaded in the module  as numpy array. If required, change the path to the data files in the code.  

2. Import the python module. If  using python 2.7 add the current folder to path allowing to import the module in the current folder.
    > import sys

    > sys.path.append("..")

If using Python3, directly import as
    > import wavenum_calbr

***When using Python IDE like Spyder***

1. After cloning the repository and moving in the `python-module` directory,  refer to the readme.  Prepare the required data as mentioned above which will be loaded in the module  as numpy array. Open the  file in the IDE and make changes  to the file path if required. Run the code.

Output
----------------

Fit  output includes coefficients for the polynomial. Using number of pixels on the CCD the wavenumber can be  generated. Required commands are included in the code.

Output for cubic polynomial is shown below:

```
Beta: [-1.04565629e+03  2.01636951e+00 -2.12196552e-04  1.85233049e-08]
Beta Std Error: [8.78769213e-02 3.69790901e-04 4.60070077e-07 1.69164597e-10]
Beta Covariance: [[ 7.57752176e-03 -3.13681344e-05  3.80574949e-08 -1.36759109e-11]
 [-3.13681344e-05  1.34180673e-07 -1.65904859e-10  6.02941948e-14]
 [ 3.80574949e-08 -1.65904859e-10  2.07694741e-13 -7.61125428e-17]
 [-1.36759109e-11  6.02941948e-14 -7.61125428e-17  2.80799598e-20]]
Residual Variance: 1.0191133117725222
Inverse Condition #: 9.292532597737295e-05
Reason(s) for Halting:
  Sum of squares convergence

Fit coefs [-1.04565629e+03  2.01636951e+00 -2.12196552e-04  1.85233049e-08]
```

Matplotlib is used to generate the plots of fit and  residuals.

<p align="center">
  <img src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/wavenum_calbr_fit_py.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/wavenum_calbr_fit_py.png" width="392" height="265" />
</p>

<p align="center">
  <img src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/wavenum_calbr_fitResdy_py.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/wavenum_calbr_fitResdy_py.png" width="392" height="265" />
</p>

<p align="center">
  <img src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/wavenum_calbr_fitResdx_py.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/wavenum_calbr_fitResdx_py.png" width="392" height="265" />
</p>
