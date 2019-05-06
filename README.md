# RamanSpec : IntensityCalibration

Set of functions in Python and IgorPro's scripting language for the determination of the the wavelength dependent sensitivity of the Raman spectrometer for intensity calibration. This repository requires the data on the Ji, frequency, and the measured rotational Raman intensities from H<sub>2</sub>, HD and D<sub>2</sub> (and O<sub>2</sub>) if available. Programs in Python and IgorPro are independent and perform the same job. The main scheme of the code is for the non-linear weighted minimization to obtain coefficients for a polynomial which represents the wavelength dependent sensitivity. An independent validation of the obtained sensitivity should be done for a measure of accuracy. 

## Principle

## Input data required
 
 - Stokes and anti-Stokes intensity data from measured rotational Raman spectra of H<sub>2</sub>, HD and D<sub>2</sub> (and O<sub>2</sub>) if available.

| Ji | Expt_Stokes/antiStokes intensity | Theory_Stokes/antiStokes intensity | Stokes_Freq | anti-Stokes freq | weight |
|----|----------------------------------|-------|-----------------|-------------|-----------------|
|    |                                  |       |                 |             |                 |
|    |                                  |       |                 |             |                 |


## Available programs
 - Set of Igor Procedures
 - A Python module 

 
## Usage

Clone the repository or download the zip file. As per your choice of the programming environment ( Python or IgorPro) refer to the specific README inside the folders and proceed. 
 
## Comments
 
 - On convergence of the minimization scheme : The convergence of the optimization has been tested with artificial and actual data giving expected results. However, in certain cases convergence in the minimization may not be achieved based on the specific data set and the error in the intensity.
 
 - Accuracy of the calibration : It is suggested to perform an independent validation of the intensity calibration. This validation can be using anti-Stokes to Stokes intensity for determining the sample temperature and calculating the depolarization ratio from spectra which is corrected for the wavelength dependent sensitivity.
 

## Credits
 
 

| First Header  | Second Header |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  

