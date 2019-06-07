# RamanSpec : IntensityCalibration

Set of functions in Python and IgorPro's scripting language for the determination of the the wavelength dependent sensitivity of the Raman spectrometer for intensity calibration. This repository requires the data on the rotational state (J), frequency, and the measured rotational Raman intensities from H<sub>2</sub>, HD, D<sub>2</sub>, N<sub>2</sub> and O<sub>2</sub>. Programs in Python and IgorPro are independent and perform the same job. The main scheme of the code is for the non-linear weighted minimization to obtain coefficients for a polynomial which represents the wavelength dependent sensitivity. An independent validation of the obtained sensitivity should be done for a measure of accuracy. 

## Principle
Refer to the following research papers for the principle. 

## Input data required
 - List of all data required : rotational state (J), experimental band area ratio (Stokes/ anti-Stokes), theoretical band area ratio (Stokes/anti-Stokes), transition frequency (Stokes) in cm<sup>-1</sup>, transition frequency (anti-Stokes) in cm<sup>-1</sup> and the weight (used for fit)
 - All of the above correspond to pair of observed bands originating from a common rotational state. For several such pairs the data is to be structured in the following format in a txt file (without header).

| Ji | Experimental (Stokes/antiStokes intensity) | Theoretical (Stokes/antiStokes intensity) | Stokes_Freq | anti-Stokes freq | weight |
|----|----------------------------------|-------|-----------------|-------------|-----------------|
|    |                                  |       |                 |             |                 |
|    |                                  |       |                 |             |                 |

*For  O<sub>2</sub>, when using the vibration-rotation transitions (S1- and O1-branch), include the data and the frequencies for these transitions.*

See specific program's readme regarding the use of the above strucutured data in the program for fit. 

## Available programs
 - Set of Igor Procedures
 - A Python module 
for performing non-linear fit on the above mentioned data set to obtain the wavelength dependent sensitivity. Additionaly, programs to compute the theoretical pure rotational Raman spectra (for H<sub>2</sub>, HD and D<sub>2</sub>) are also included.

 
## Usage

Clone the repository or download the zip file. As per your choice of the programming environment ( Python or IgorPro) refer to the specific README inside the folders and proceed. 
 
## Comments
 
 - On convergence of the minimization scheme : The convergence of the optimization has been tested with artificial and actual data giving expected results. However, in certain cases convergence in the minimization may not be achieved based on the specific data set and the error in the intensity.
 
 - Accuracy of the calibration : It is suggested to perform an independent validation of the intensity calibration. This validation can be using anti-Stokes to Stokes intensity for determining the sample temperature and calculating the depolarization ratio from spectra which is corrected for the wavelength dependent sensitivity.
 

## Credits
Non-linear optimization in SciPy : 
 - Travis E. Oliphant. Python for Scientific Computing, Computing in Science & Engineering, 9, 10-20 (2007), DOI:10.1109/MCSE.2007.58 

