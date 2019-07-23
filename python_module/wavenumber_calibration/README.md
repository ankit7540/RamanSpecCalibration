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
