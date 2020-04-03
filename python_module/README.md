# Python module for Raman calibration

This module is broadly divided  into two sections , intensity and wavenumber calibration.

**Intensity  calibration :**  Includes following files,

`compute_spectra.py` : This is an independent module which imports energy state data (for H<sub>2</sub>, HD and D<sub>2</sub>) along with ratio of polarizability anisotropy matrix  elements (at 532.2 nm, for our experiment) to compute the theoretical rotational Raman spectra of these gases. Temperature is included for Boltzmann population using which observable Raman spectra can be modeled. For obtaining the ratio of theoretical Raman intensities originating from a common rotational state (required for the `wavelength_sensivity.py` module ), taking a ratio of the computed intensities from common  rotational  states is required. This ratio is independent of the  Boltzmann  populations.

For matrix elements of polarizability anisotropy for H<sub>2</sub>, HD and D<sub>2</sub>, at other wavelengths refer to this GitHub repository (<https://github.com/ankit7540/H2-PolarizabilityMatrixElements>) and research paper (<https://aip.scitation.org/doi/abs/10.1063/1.5011433>).

`wavelength_sensivity.py` : It accepts the data of experimental rotational Raman intensity (from common rotational state) for performing a fit to determine the coefficients of a polynomial representing the wavelength dependent sensitivity.

**Wavenumber calibration :**  Includes following file,

`wavenumber_cal.py` : It requires data on the band positions in pixels with the uncertainties, as vectors in a 1D txt file. (NumPy array can be directly used as well). Corresponding reference data is also required as 1D arrays containing values and the errors. Core program is just fitting the pixel with the reference wavelength data. Example data sets are included.


Requirements
----------------
Python 2.7 or Python 3.x with NumPy, SciPy, math and matplotlib modules.

Usage
----------------
Following commands are run under the Python interpreter environment. Alternatively use any Python IDE like Spyder, IDLE or PyCharm). IDE are best suited since the obtained results can be immediately viewed using the matplotlib plots.
