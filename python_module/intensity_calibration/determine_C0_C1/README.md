# Generating C<sub>0</sub>/C<sub>1</sub> correction

The module `gen_correction.py` allows one to generate the correction (C<sub>0</sub>/C<sub>1</sub>)
using the function `gen_C0_C1`.

Here, the C<sub>0</sub> correction is obtained from the wavenumber axis supplied by the user. The C<sub>1</sub> correction is determined from a broadband white light spectrum, which is assumed to be a Black-Body emitter (standard tungsten lamps work well for this purpose).

---

For a similar program written in IgorPro see [gen_correction.ipf](https://github.com/ankit7540/RamanSpec_BasicOperations/blob/master/intensity_corr/) in my other repository.


## Requirements

Python 3.6 or higher with numpy, scipy and matplotlib modules.


## Function

```
>>> import gen_correction
        **********************************************************

         This module is for generating the wavenumber-dependent
         intensity correction curve.

         Main function is listed below.

        **********************************************************
        gen_C0_C1 ( Ramanshift,  laser_nm, wl_spectra, norm_pnt,
                        mask = None, set_mask_nan = None, export = None)
        **********************************************************

         REQUIRED PARAMETERS
                         Ramanshift = vector, the x-axis in relative wavenumbers
                         laser_nm = scalar, the laser wavelength in nanometers
                         wl_spectra = broadband whitelight spectra (1D or 2D)
                         norm_pnt =  scalar, normalization point (corrections will be set
                                        to unity at this point
          OPTIONAL PARAMETERS
                         mask = vector, mask wave for selecting specific region to fit
                         set_mask_nan= boolean, 1 will set the masked region in
                                        the output correction to nan, 0 will not do so.
                         export = 0 or 1, setting to 1 will export the correction as a txt
                                     file with name intensity_correction.txt
                          ------------------------------------------
                          All vectors required here should be numpy arrays.
                          See line 14 to 18 to define/load the numpy arrays
                                              before execution.
        **********************************************************

```

## Input parameters

 - Ramanshift =  numpy array, 1D
 - laser_nm = scalar
 - wl_spectra =  numpy array, 1D or 2D
 - norm_pnt =  scalar (index where the normalization is performed)
               Should be a number between 0 and (n-1) where n is the number of rows
               in the wavenumber axis.
 - mask = numpy array, 1D, boolean ( 0 for points which are to be fit and 1
                                     for the points to be masked )
 - set_mask_nan = scalar 0 or 1
 - export = 0 or 1


## Examples

Check this jupyter notebook showing an example with data. [Examples/Example_for_C0_C1.ipynb](https://github.com/ankit7540/IntensityCalbr/blob/master/PythonModule/determine_C0_C1_correction/Examples/Example_for_C0_C1.ipynb)

### loading files

```

wavenumber = np.loadtxt("wavenumber.txt")

wl = np.loadtxt("wl.txt")

mask = np.loadtxt("mask.txt")

mask = mask.astype(bool)    # conversion to boolean

mask_array =  np.invert(mask)   # inverting bool, scipy requires masked region to be true

```

### Running the main program.

```

corr = gen_C0_C1 (wavenumber, 532.2, wl, 582, mask = mask_array, set_mask_nan = 1, export = 1)
	 Mask is available. Using mask and fitting.
	 Optimized coefs : [3.33443440e-19 3.76715151e+03]
	 Correction will be exported as intensity_correction.txt

```

---
