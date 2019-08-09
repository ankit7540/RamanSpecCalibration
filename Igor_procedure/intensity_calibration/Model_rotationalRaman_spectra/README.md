# Theoretical intensities  for rotational Raman  bands of H<sub>2</sub>, HD and D<sub>2</sub>

This section includes one file `Compute_spectra.ipf` for the computation of the rotational Raman intensities of the three gases. Temperature is required  for computing the Boltzmann population. Wavelength dependent polarizability anisotropy matrix elements are required for the transition probability of the rotational Raman transitions. As an example, the data for the 532.199323 nm, specific to our experiment, was obtained  from  this repository (<https://github.com/ankit7540/H2-PolarizabilityMatrixElements>),  will be used here.


Requirements
----------------
IgorPro 6.36 or  7 and higher.

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
All the  required  data is available in the sub-directory called `energy_levels_and_gamma` as  Igor text files (.itx)

After loading the  waves, use the commands available for computing the Sumofstates (Q) and relative intensities  of  rotational Raman bands of  H<sub>2</sub>, HD and D<sub>2</sub>

Available  Functions
-------------------
`sumofstate_H2(  T  )` Computing  partition function for the gases. Definition is as follows,<br>

<p align="center">
  <img   src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/partition_function_defn.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/partition_function_defn.png" width="581" height="85" />
</p>


`spectra_H2(T, J_{Stokes}, J_{anti-Stokes})`, and  similar for HD and D<sub>2</sub> : For computing spectra (normalized relative Raman  intensities for  pure rotation in the three gases). Here,  temperature T  is in  Kelvin. `J_{Stokes}` and `J_{anti-Stokes})` represents  the maximum rotational J whose intensity is  to be  computed. Intensities are  computed for `J=J_{anti-Stokes},.., J_{Stokes}`. Computed  intensities are relative normalized intensities.

Functions can be used directly for computing the sum of states for the three gases for given temperature (for example,  `sumofstate_D2(298)`). When computing spectra the sum of states  are computed within the function (for example, `spectra_H2(298, 6, 6)`) for obtaining the Boltzmann  populations.

Output
-------------
`SumOfstates_H2 ( 298 )` outputs a scalar.

`spectra_H2 ( 298, 5  , 5 )`  outputs  three  waves, posnH2, specH2, spectraH2.
posnH2 :  1D wave containing  the Raman transition  wavenumbers.
specH2 :  2D wave containing  rotational  state, Raman transition  wavenumbers, relative Raman intensities  and the  transition wavenumber in  absolute  wavenumbers.
spectraH2 :  1D wave containing  the relative Raman intensities.


```
•print SumOfstates_H2 ( 298 )
  7.71851
•print SumOfstates_D2 ( 458 )
  49.8149

•spectra_H2 ( 298, 5  , 5 )
•print  posnH2
  posnH2[0]= {-1034.67,-814.424,-587.032,-354.373,354.373,587.032,814.424,1034.67,1246.1,1447.28}

•print  specH2
  specH2[0][0]= {5,4,3,2,0,1,2,3,4,5}
  specH2[0][1]= {-1034.67,-814.424,-587.032,-354.373,354.373,587.032,814.424,1034.67,1246.1,1447.28}
  specH2[0][2]= {0.00100178,0.00381731,0.0708816,0.0679526,0.335826,1,0.150123,0.106331,0.00464706,0.00105514}
  specH2[0][3]= {19824.6,19604.4,19377,19144.3,18435.6,18202.9,17975.5,17755.3,17543.8,17342.7}

•print  spectraH2
    spectraH2[0]= {0.00100178,0.00381731,0.0708816,0.0679526,0.335826,1,0.150123,0.106331,0.00464706,0.00105514}

```

Similarly  for HD and D<sub>2</sub>.
<br>

Plot  showing modelled spectra of the three gases.</br>
<p align="center">
  <img   src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/Igor_plot.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/Igor_plot.png" width="598" height="385" />
</p>
