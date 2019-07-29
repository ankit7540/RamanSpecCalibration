# Python module for obtaining wavelength sensitivity for Raman intensity calibration

This module includes `wavelength_sensitivity.py` for performing fit using experimental band ratios to obtain the wavelength sensitivity modelled as a polynomial.

Usage
----------------
Required input data should be structured as following.

2D waves with columns  containing the data as follows.

|   J   |  Expt  ratio |  theoretical ratio |   J+2  transition wavenumber  (cm<sup>-1</sup>) |  J-2  transition wavenumber (cm<sup>-1</sup>) | weight  |

For example,   (header  given for  clarity). Data wave as example is included in the present directory  as file, `dataD2.txt`

| J | Expt Ratio | Theoretical  Ratio | J+2, freq | J-2, freq | Weight |
|---|     ------------|-----------------      |-----------   |-----------   |--------    |
| 2 | 2.4872     | 2.3827               | 414.65    | -179.07   | 11.67   |
| 3 | 1.7910     | 1.6652               | 529.90    | -297.53   | 9.48    |
| 4 | 1.4674     | 1.3889               | 642.81    | -414.65   | 10.72  |
| 5 | 1.3904     | 1.2331               | 752.92    | -529.90   | 4.49    |
| 6 | 1.2955     | 1.1287               | 859.83    | -642.81   | 3.26    |


There are three functions for computing the residuals, `residual_linear`, `residual_quadratic` and `residual_cubic`. The  path to the data waves in the specific usage is to be modified. See following code block,
```
#********************************************************************
# Load data files
# Modify this as user case
dataH2 = np.loadtxt("./dataH2.txt")
dataHD = np.loadtxt("./dataHD.txt")
dataD2 = np.loadtxt("./dataD2.txt")
dataO2 = np.loadtxt("./dataO2.txt")
dataO2_p = np.loadtxt("./dataO2_p.txt")
xaxis = np.loadtxt("./Ramanshift_axis.txt")

#********************************************************************
```


Requirements
----------------
Python 2.7 or Python 3.x with NumPy, SciPy and math modules

Usage
----------------
Following commands are run under the Python interpreter environment. (Alternatively use any Python IDE like Spyder, IDLE or PyCharm). *Spyder3 is has been used while writing and debugging  python  codes given  here.*

***When using Python interpreter in terminal***

1. After cloning the repository and moving in the `python-module` directory,  refer to the readme.  Prepare the required data as mentioned above which will be loaded in the module  as numpy array. If required, change the path to the data files in the code.  

2. Import the python module. If  using python 2.7 add the current folder to path allowing to import the module in the current folder.

```
   import sys

   sys.path.append("..")
```

If using Python3, directly import as

  ```
     import wavelength_sensitivity
  ```

***When using Python IDE like Spyder***

1. After cloning the repository and moving in the `python-module` directory,  refer to the readme.  Prepare the required data as mentioned above which will be loaded in the module  as NumPy array. Open the  file in the IDE and make changes  to the file path if required. Run the code.
