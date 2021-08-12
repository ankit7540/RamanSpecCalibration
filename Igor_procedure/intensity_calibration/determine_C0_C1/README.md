# Obtain C<sub>0</sub>/C<sub>1</sub>

Intensity calibration of Raman spectra based on a three step procedure has been established in our research. (See [10.1002/jrs.5955](https://onlinelibrary.wiley.com/doi/full/10.1002/jrs.5955) for more details)

These three corrections for the intensity calibration steps are labelled as C<sub>0</sub>, C<sub>1</sub> and C<sub>2</sub>.

 - C<sub>0</sub> is derived from the vector representing the x-axis.
 - C<sub>1</sub> comes from the fit of the broadband white-light spectrum ( which is corrected for C<sub>0</sub>) to the number of photons from a black-body emitter.
 - and C<sub>2</sub> comes from the Raman intensities (which are already corrected using C<sub>0</sub> and C<sub>1</sub>. This correction rectifies any inconsistencies not covered in the previous steps of intensity calibration. See the above paper for more details.  
The total multiplicative correction is given as  C<sub>0</sub>/(C<sub>1</sub> * C<sub>2</sub>).


The procedure file in this folder, `gen_correction.ipf` has the function defined as `gen_C0_C1 ` for generating the first two corrections, and which returns : (C<sub>0</sub>/C<sub>1</sub>). For common Raman measurements, where very high accuracy in the relative Raman intensities is not required, this might suffice.
