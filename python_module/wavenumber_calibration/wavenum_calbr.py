# -*- coding: utf-8 -*-

"""Module describing the weighted fit to
determine the wavenumber calibration of the spectrometer using a  polynomial
as a model function"""

# Load necessary modules
import numpy as np
import scipy.odr
import matplotlib.pyplot as plt

#********************************************************************

# Load data files : MODIFY  THE FOLLOWING FOR THE INPUT DATA

#   ref_freq    : reference  wavenumber, 1D array
#   ref_error   : error in the reference, 1D array
#   pixel       : pixel position  of bands, 1D array
#   pixel_error : error in band psoition

# Above are loaded from txt files in the present case as example

ref_freq = np.loadtxt("./freq.txt")
ref_error = np.loadtxt("./freq_sigma.txt")
pixel = np.loadtxt("./pixel.txt")
pixel_error = np.loadtxt("./pixel_sigma.txt")

#********************************************************************
# Fit function defined here
def fit_func_quadratic(B, x):
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0] + B[1]*x  + B[2]*(x**2)

#********************************************************************


#********************************************************************
# Fit function defined here
def fit_func_cubic(B, x):
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0] + B[1]*x  + B[2]*(x**2)  + B[3]*(x**3)

#********************************************************************

quadraticPolynomial = scipy.odr.Model(fit_func_quadratic)
cubicPolynomial = scipy.odr.Model(fit_func_cubic)

Data = scipy.odr.RealData(pixel, ref_freq ,  sx=pixel_error, sy=ref_error)

# Change initial coefs defined below as required
Fit_q = scipy.odr.ODR(Data, quadraticPolynomial, [-1000, 2,  1e-4 ])
Fit_c = scipy.odr.ODR(Data, cubicPolynomial, [-1000, 2,  1e-4, 1e-5])


# Cubic fit is shown here. Uncomment the line for quadratic if required.
# Quadratic fit  needs change on line numbers: 65,66 ; 80,81 ; 83,84

#output = Fit_q.run()
output = Fit_c.run()

output.pprint()

ODR_result =  output.beta
ODR_std = output.sd_beta
ODR_resd_x = output.delta
ODR_resd_y = output.eps

print ("\nFit coefs", ODR_result)

#  For 1600  pixels on  the detector. Change 1600 to  your case if required.
xpoints = np.arange(0, 1600, dtype=float)

#fit = fit_func_quadratic( ODR_result, xpoints)
fit = fit_func_cubic( ODR_result, xpoints)

#wavenumber_axis = fit_func_quadratic( ODR_result, xpoints)
wavenumber_axis = fit_func_cubic( ODR_result, xpoints)

#print ("ODR residual in x", ODR_resd_x)
#print ("ODR residual in y", ODR_resd_y)

#********************************************************************

# Plotting the data :  MATPLOTLIB REQUIRED

txt = ("*Generated from 'wavenum_calbr.py' on the\
      \nGitHub Repository: RamanSpecCalibration ")

# FIGURE 0 INITIALIZED

plt.figure(0)
ax0 = plt.axes()
plt.title('Fitting result', fontsize=20)
plt.plot( pixel,  ref_freq, 'ro',  label='wavenumber (ref)')
plt.plot( xpoints, fit, 'b', linewidth=2, label='cubic_fit')

plt.xlabel('pixel', fontsize=16)
plt.ylabel('Wavenumber', fontsize=16)
plt.grid(True)
ax0.tick_params(axis='both', labelsize =15)

ax0.minorticks_on()
ax0.tick_params(which='minor', right='on')
ax0.tick_params(axis='y', labelleft='on', labelright='on')
plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey",\
         transform=plt.gcf().transFigure)
plt.legend(loc='upper left', fontsize=15)

#plt.savefig('wavenum_calbr_fit_py.png', bbox_inches='tight', dpi=300)
#********************************************************************

ymax = np.amax(np.absolute(ODR_resd_y))
xmax = np.amax(np.absolute(ODR_resd_x))

#********************************************************************

# FIGURE 1 INITIALIZED
plt.figure(1)
ax1 = plt.axes()
plt.title('Residual of the fit')
plt.plot( pixel, ODR_resd_y, "ro" , label="residual in y")
plt.xlabel('pixel', fontsize=16)
plt.ylabel('residual in  y')
plt.grid(True)
plt.legend(loc='upper left')
ax1.minorticks_on()
ax1.set_ylim([-1*ymax-0.025, ymax+0.025])
plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey", transform=plt.gcf().transFigure)
ax1.minorticks_on()
ax1.tick_params(which='minor', right='on')
ax1.tick_params(axis='y', labelleft='on', labelright='on')
ax1.tick_params(axis='both', labelsize =15)
ax1.tick_params(axis='both', labelsize =15)

#plt.savefig('wavenum_calbr_fitResdy_py.png', bbox_inches='tight', dpi=300)
#********************************************************************


# FIGURE 2 INITIALIZED
plt.figure(2)
ax1 = plt.axes()
plt.title('Residual of the fit')
plt.plot( pixel, ODR_resd_x, "bo", label="residual in x")
plt.xlabel('pixel', fontsize=16)
plt.ylabel('residual in x')
plt.grid(True)
plt.legend(loc='upper left')
ax1.minorticks_on()
ax1.set_ylim([-1*xmax-0.1, xmax+0.1])
plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey", transform=plt.gcf().transFigure)
ax1.minorticks_on()
ax1.tick_params(which='minor', right='on')
ax1.tick_params(axis='y', labelleft='on', labelright='on')
ax1.tick_params(axis='both', labelsize =15)
ax1.tick_params(axis='both', labelsize =16)

#plt.savefig('wavenum_calbr_fitResdx_py.png', bbox_inches='tight', dpi=300)
#********************************************************************
