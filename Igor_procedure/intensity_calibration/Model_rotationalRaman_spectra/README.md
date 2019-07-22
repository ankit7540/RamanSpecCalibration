# Theoretical intensities  for rotational Raman  bands of H<sub>2</sub>, HD and D<sub>2</sub>

This section includes one file `Compute_spectra.ipf` for the computation of the rotational Raman intensities of the three gases. Temperature is required  for computing the Boltzmann population. Wavelength dependent polarizability anisotropy matrix elements are required for the transition probability of the rotational Raman transitions. As an example, the data for the 532.199323 nm, specific to our experiment, was obtained  from  this repository (<https://github.com/ankit7540/H2-PolarizabilityMatrixElements>),  will be used here.


Requirements
----------------
Tested with  Igor Pro  7.

Usage
----------------
Following commands is run in  the  command browser to initialize  the procedure on first run.

```
rotational_analysis()
  Folder (root:Params:RotnRamanIns)  generated. Please load the data files in the following folder.
  root:Params:RotnRamanIns
```

Above  makes  new  folder  for keeping data for the energy of rovibrational states and the Polarizability anisotropy Matrix Elements.

Next,  load the data  waves to  the folder, `root:Params:RotnRamanIns`
All the  required  data is available in the sub-directory called `energy_levels_and_gamma`

After loading the  waves, use the commmands available for computing the sumofstates (Q) and relative intensities  of  rotational Raman bands of  H<sub>2</sub>, HD and D<sub>2</sub>

Available  Functions
-------------------
`sumofstate_H2(  T  )` Computing  partition function for the gases. Definition is as follows,
<img src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/partition_function_defn.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/partition_function_defn.png" width="854" height="125" />

`spectra_H2(T, J_{Stokes}, J_{anti-Stokes})`, and  similar for HD and D<sub>2</sub> Here,  temperature T  is in  Kelvin. `J_{Stokes}` and `J_{anti-Stokes})` represents  the maximum rotational J whose intensity is  to be  computed. Intensities are  computed for `J=J_{anti-Stokes},.., J_{Stokes}`. Computed  intensities are relative normalized intensities.

Functions can be used directly for computing the sum of states for the three gases for given temperature (for example,  `sumofstate_D2(298)`). When computing spectra the sum of states  are computed within the function (for example, `spectra_H2(298, 6, 6)`).

Output
-------------





```
•print SumOfstates_H2 ( 298 )
  7.71851
•print SumOfstates_D2 ( 458 )
  49.8149

•spectra_H2 ( 298, 5  , 5 )
•spectra_HD ( 298, 5  , 5 )
•spectra_D2 ( 298, 7  , 7 )



```
<br>

Using  matplotlib  spectra can be plotted using the band position and the intensities. For example,<br>
<img src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/python_module/intensity_calibration/model_rotationalRaman_spectra/spectra.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/python_module/intensity_calibration/model_rotationalRaman_spectra/spectra.png" width="392" height="265" />
