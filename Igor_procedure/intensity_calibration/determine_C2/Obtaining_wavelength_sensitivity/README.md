# Wavelength sensitivity  from the rotational Raman spectra of H<sub>2</sub>, HD and D<sub>2</sub>

This section includes one file `Weighted_fit_wavelength_sensitivity.ipf` for obtaining the wavelength dependent sensitivity modelled as a polynomial. Ratio  of intensities from common rotational states are compared to the corresponding theoretical ratio to obtain the wavelength dependent sensitivity curve.

Requirements
----------------
IgorPro 6.36 or  7 and higher.

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


There are three functions for computing the residuals, `Residual_Linear`, `Residual_Quadratic` and `Residual_Cubic`. The  path to the data waves in the Igor experiment file is  to be modified in the function for the residual.

See following code block  from  the residual  functions,
```
// ---------------------- waves containing  data   for  fit------------------------------
//  Modify the location of  the  data waves in your experiment
wave dataH2   = root:dataH2
wave dataD2  = root:dataD2
wave dataHD  = root:dataHD
wave dataO2_pR  = root:DataO2_p
wave dataO2_o1s1  = root:DataO2
// --------------------------------------------------------------------------------------

```

Next, create a temporary wave (required  for passing to the optimize function)  and  then  run the  function, `run_Opt_Linear`, `run_Opt_Quadratic` or  `run_Opt_Cubic` with guess  coefficients.

```
make  /n=1 temp=0 // this wave is passed to the optimize function, not used in present implementation
run_Opt_Linear  ( temp , 0.115   )
run_Opt_Quadratic  ( temp , 0.115 , 1e-2  )
run_Opt_Cubic  ( temp , 0.115 , 1e-2,  1e-4 )
```
Using the obtained coefficients  and the wavenumber axis  vector (1D  wave for  the xaxis)  generate the sensitivity curve.
