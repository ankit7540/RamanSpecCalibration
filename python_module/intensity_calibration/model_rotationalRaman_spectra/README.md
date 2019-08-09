# Theoretical intensities  for rotational Raman  bands of H<sub>2</sub>, HD and D<sub>2</sub>

This section includes one file `compute_spectra.py` for the computation of the rotational Raman intensities of the three gases. Temperature is required  for computing the Boltzmann population. Wavelength dependent polarizability anisotropy matrix elements are required for the transition probability of the rotational Raman transitions. As an example, the data for the 532.199323 nm, specific to our experiment, was obtained  from  this repository (<https://github.com/ankit7540/H2-PolarizabilityMatrixElements>),  will be used here.


Requirements
----------------
Python 2.7 or Python 3.x with NumPy  and  matplotlib

Usage
----------------
Following commands are run under the Python interpreter environment. (Alternatively use any Python IDE like Spyder, IDLE or PyCharm). *Spyder3 is has been used while writing and debugging.*

***When using Python interpreter in terminal***

1. After cloning the repository and moving in the `python-module` directory,  refer to the readme.

2. Import the python module. If  using python 2.7 add the current folder to path allowing to import the module in the current folder.
    > import sys

    > sys.path.append("..")

If using Python3, directly import as
    > import compute_spectra

***When using Python IDE like Spyder***

1. After cloning the repository and moving in the `python-module` directory,  refer to the readme. Open the file in the IDE and make changes to the file path if required. Run the code.

Available  Functions
-------------------
`sumofstate_H2(  T  )` and similar for HD and D<sub>2</sub> :  The sum of state or the partition function is defined as <br>

<p align="center">
  <img   src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/partition_function_defn.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/partition_function_defn.png" width="581" height="85" />
</p>

`spectra_H2(T, J_{Stokes}, J_{anti-Stokes})` and similar for HD and D<sub>2</sub> : For computing the relative Raman intensities (for pure rotation) for the specific gas.
Here,  temperature T  is in  Kelvin. `J_{Stokes}` and `J_{anti-Stokes})` represents  the maximum rotational J whose intensity is  to be  computed. Intensities are  computed for `J=J_{anti-Stokes},.., J_{Stokes}`. Computed  intensities are relative normalized intensities.

Functions can be used directly for computing the sum of states for the three gases for given temperature (for example,  `sumofstate_D2(298)`). When computing spectra the sum of states  are computed within the function (for example, `spectra_H2(298, 6, 6)`).

Output
-------------


```
 Sum of state for  H2 at 333 K :  8.547700236069472
 Sum of state for  HD at 375 K :  6.220035416830953
 Sum of state for  D2 at 298 K :  32.86935389137424


Normalized spectra H2:
 [[1.68807110e-05]
 [1.00177858e-03]
 [3.81730673e-03]
 [7.08815800e-02]
 [6.79526166e-02]
 [3.35825806e-01]
 [1.00000000e+00]
 [1.50122643e-01]
 [1.06331157e-01]
 [4.64706587e-03]
 [1.05514387e-03]
 [1.58920099e-05]]


Normalized spectra HD:
 [[4.82808002e-05]
 [6.71015239e-04]
 [6.22339903e-03]
 [3.72083472e-02]
 [1.35642025e-01]
 [2.60840366e-01]
 [8.69625413e-01]
 [1.00000000e+00]
 [5.98423539e-01]
 [2.14344369e-01]
 [4.83621415e-02]
 [7.08669001e-03]
 [6.90790930e-04]
 [4.57832852e-05]]


Normalized spectra D2 :
 [[2.76095983e-04]
 [1.03622463e-03]
 [1.18052776e-02]
 [2.50459482e-02]
 [1.54194158e-01]
 [1.64009643e-01]
 [4.19691073e-01]
 [9.40953417e-01]
 [6.27306694e-01]
 [1.00000000e+00]
 [2.73110763e-01]
 [2.14165350e-01]
 [3.08832831e-02]
 [1.33245422e-02]]
```
<br>

Using  matplotlib  spectra can be plotted using the band position and the intensities. For example,<br>

<p align="center">
  <img src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/compute_spectra_python.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/compute_spectra_python.png" width="392" height="265" />
</p>
