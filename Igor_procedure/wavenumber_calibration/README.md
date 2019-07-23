# Igor procedure for  wavenumber calibration

 `wavenumber_calibration.ipf` for performing fit for obtaining wavenumber  axis. Data on the band positions in pixels with the uncertainties, as waves in the Igor experiment. Core program is just fitting the pixel with the reference wavelength data. It includes  the function `wavenumber_fit` which  performs the fit.



`data_waves.itx`  in the present  folder includes  4 waves, `freq`, `freq_sigma`, `pixel` and `pixel_sigma` which  serve as example for using the procedure.

Usage
----------------
```
wavenumber_fit ( wavenum_reference,
                 error_in_reference,  
                 band_position_in_pixel,
                 error_band_posn ,
                 num_pixel_on_CCD )

```
`wavenum_reference`       =   1D wave, list  of  reference Raman transition wavenumbers

`error_in_reference`      =   1D wave, error in the above Raman transition wavenumbers

`band_position_in_pixel`  =   1D wave, corresponding band  positions on  the CCD (in  pixel)

`error_band_posn`         =   1D wave, error in the band positions

`num_pixel_on_CCD`        =   variable,  number of pixels  on the CCD,  for example, 1024, 1340  or 1600. This is used for generating  the wavenumber (relative) axis.


Example
-------------------
```
•wavenumber_fit ( freq, freq_sigma,  pixel, pixel_sigma , 1600 )
  Fit converged properly
  fit_freq= poly(W_coef,x)
  residual_y= freq[p] - (poly(W_coef,pixel[p]))
  W_coef={-1045.7,2.0164,-0.00021219,1.8518e-008}
  V_chisq= 47.7758;V_npnts= 50;V_numNaNs= 0;V_numINFs= 0;
  V_startRow= 0;V_endRow= 50;
  W_sigma={0.0871,0.000367,4.57e-07,1.68e-10}
  Coefficient values ± one standard deviation
  	K0	=-1045.7 ± 0.0871
  	K1	=2.0164 ± 0.000367
  	K2	=-0.00021219 ± 4.57e-007
  	K3	=1.8518e-008 ± 1.68e-010
```


A  plot  is generated in the end as shown  below.

<p align="center">
  <img src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/wavenum_calbr_fit_ig.png" data-canonical-src="https://github.com/ankit7540/RamanSpecCalibration/blob/master/img/wavenum_calbr_fit_ig.png" width="505" height="325" />
</p>


Error  analysis
-----------------------
For  error  analysis refer to the manuscript and supplementary material.
