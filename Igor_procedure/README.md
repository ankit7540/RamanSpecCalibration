# IgorPro module for Intensity calibration

This module includes two file having several functions each, `RamanCalibration.ipf`, `compute_spectra.ipf`

`RamanCalibration.ipf` : Includes two major functions for wavenumber and intensity calibration. For wavenumber calibration : Data on the band positions in pixels with the uncertainties, as waves in the Igor experiment. Core program is just fitting the pixel with the reference wavelength data. The reference data in already included in the procedure with the error. For intensity calibration : Data on the band intensities and frequencies from common rotational states for different gas species arranged in a  predefined way as 2D waves. See **Usage** section for details.

`compute_spectra.ipf` : This is an independent Igor procedure which uses energy state data (for H<sub>2</sub>, HD and D<sub>2</sub>) along with ratio of polarizability anisotropy matrix (at 532.2 nm) elements to compute the theoretical rotational Raman spectra of these gases. Temperature is included for Boltzmann population using which observable Raman spectra can be modeled. For obtaining the ratio of theoretical Raman intensities originating from a common rotational state (required for the `RamanCalibration.ipf` module ), a function is included  here which simply takes the ratio of the intensities of the individual bands.

For matrix elements of polarizability anisotropy for H<sub>2</sub>, HD and D<sub>2</sub>, at other wavelengths refer to this GitHub repository (<https://github.com/ankit7540/H2-PolarizabilityMatrixElements>) and research paper (<https://aip.scitation.org/doi/abs/10.1063/1.5011433>).


Requirements
----------------
Tetsed with IgorPro 6.36 and 7 or higher.

Usage
----------------
Data required is assumed to be included in the Igor experiment as waves.

For wavenumber calibration, two waves `band position (in pixel)` with `error in band position` are required.

For intensity calibration, dataset for each specific gas is required as a 2D wave containing rotational state (J), experimental band area ratio (Stokes/ anti-Stokes), theoretical band area ratio (Stokes/anti-Stokes), transition frequency (Stokes) in cm<sup>-1</sup>, transition frequency (anti-Stokes) in cm<sup>-1</sup> and the weight (used for fit). For  O<sub>2</sub>, when using the vibration-rotation transitions (S1- and O1-branch), include the data and the frequencies for these transitions. (All of the above correspond to pair of observed bands originating from a common rotational state.)
