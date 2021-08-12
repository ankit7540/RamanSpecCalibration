#############################################################################

import package_util

import numpy as np
from numpy.polynomial import Polynomial
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#plt.rcParams["font.family"] = "Arial"

#############################################################################

# define input arrays here

#axis=np.loadtxt('axis.txt')  # xaxis, in relative wavenumbers
#wl=np.loadtxt('wl.txt')      # white light spectra, 1D or 2D
#mask=np.loadtxt('mask.txt')  # mask, integers, set 1 for masked region
#mask = mask.astype(bool)     # convert mask to boolean
#mask_array =  np.invert(mask)  # use this if mask is required to be inverted

#############################################################################

print('\t**********************************************************')

print('\t ')
print('\t This module is for generating the wavenumber-dependent')
print('\t intensity correction curve.')

print('\n\t Main function is listed below.')
print('\t ')
print('\t**********************************************************')
print('\tgen_C0_C1 ( Ramanshift,  laser_nm, wl_spectra, norm_pnt,  ')
print ( '\t\t        mask = None, set_mask_nan = None, export = None) ')
print('\t**********************************************************')
print('\n\t REQUIRED PARAMETERS')
print('\t\t\t Ramanshift = vector, the x-axis in relative wavenumbers')
print('\t\t\t laser_nm = scalar, the laser wavelength in nanometers')
print('\t\t\t wl_spectra = broadband whitelight spectra (1D or 2D)')
print('\t\t\t norm_pnt =  scalar, normalization point (corrections will be set')
print('\t\t\t                to unity at this point')
print('\t  OPTIONAL PARAMETERS')
print('\t\t\t mask = vector, mask wave for selecting specific region to fit')
print('\t\t\t set_mask_nan= boolean, 1 will set the masked region in')
print('\t\t\t                the output correction to nan, 0 will not do so.')
print('\t\t\t export = 0 or 1, setting to 1 will export the correction as a txt')
print('\t\t\t             file with name intensity_correction.txt')

print('\t\t\t  ------------------------------------------')
print('\t\t\t  All vectors required here should be numpy arrays.')
print('\t\t\t  See line 14 to 18 in the file to define/load the numpy arrays')
print('\t\t\t                      before execution.')

print('\t**********************************************************')

print('\t\t\t  Example:')

print('\t\t\t  gen_C0_C1 (Ramanshift, 532, wl_spectra, 500 )')
print('\t\t\t  gen_C0_C1 (Ramanshift, 532, wl_spectra, 500, mask=maskw, ')
print('\t\t\t\t\t\t                  set_mask_nan=0, export=1 )')

print('\t**********************************************************')

#############################################################################


#############################################################################
# This file has function(s) for determination of the C0 and C1 corrections

def gen_C0_C1 (Ramanshift, laser_nm, wl_spectra, norm_pnt, mask = None,
               set_mask_nan = None, export = None):

    '''Ramanshift = vector, the x-axis in relative wavenumbers
       laser_nm = scalar, the laser wavelength in nanometers,
       wl_spectra = broadband whitelight spectra (1D or 2D),
       norm_pnt =  normalization point (corrections will be set
                                        to unity at this point),

       OPTIONAL PARAMETERS
       mask = vector, mask wave for selecting specific region to fit ,
       set_mask_nan= boolean, 1 will set the masked region in
       the output correction to
          nan, 0 will not '''


    C0 = gen_C0 (Ramanshift, norm_pnt)

    # check the whitelight spectrum provided
    dim = wl_spectra.shape
    if(dim[0] != Ramanshift.shape[0]):
        print ("\t Error : Dimension mismatch for wl spectra and the xaxis")

    if (wl_spectra.ndim > 1    ):
        print ("\t wl spectra is 2D")
        wl_spectra = np.mean(wl_spectra, axis=1)
        #print(wl_spectra.shape)

    # normalize the wl
    wl_norm = wl_spectra / np.amax(wl_spectra)

    # correct wl with C0
    wl_norm_C0 = np.multiply(wl_norm , C0)

    # check if mask supplied --------
    if (isinstance(mask, np.ndarray) != 1):
        print ("\t Mask not available. Proceeding with fit..")
        C1 = gen_C1 (Ramanshift, laser_nm ,  wl_norm_C0 ,   norm_pnt)

    elif (isinstance(mask, np.ndarray) == 1):
        print ("\t Mask is available. Using mask and fitting.")
        C1 = gen_C1_with_mask (Ramanshift, laser_nm ,  wl_norm_C0 , mask,
                               norm_pnt)

        if (set_mask_nan==1 ):
            C1 [mask] = np.nan

    correction = (C0/C1)
    #--------------------------------
    # export output
    if (export == 1):
        print('\t Correction will be exported as intensity_correction.txt')

        # export the correction as txt
        filename = 'intensity_correction.txt'
        np.savetxt(filename, correction, fmt='%3.7f', newline='\n',
                   header='intensity_corr')

    return correction
    #----------------------------------------------------------


#############################################################################

def gen_C0 (Ramanshift, norm_pnt):
    '''Ramanshift = vector, the x-axis in wavenumbers
       norm_pnt =  normalization point (corrections will be set
                                        to unity at this point) '''

    spacing = np.diff(Ramanshift)

    # normalization
    wavenum_corr  = spacing / spacing [norm_pnt]
    temp_xaxis = np.arange(spacing.shape[0])


    # fitting the wavenumber diff
    c = Polynomial.fit(temp_xaxis , wavenum_corr, deg=3)

    # extrpolate the last point
    wavenum_corr = np.append(wavenum_corr, c(spacing.shape[0]+1))

    return wavenum_corr


#############################################################################

def photons_per_unit_wavenum_abs(x,a,T) :
    return (a*599584916*(x**2))/(np.exp(0.1438776877e-1*x/T)-1)

#############################################################################
#############################################################################

def gen_C1 (Ramanshift, laser_nm ,  wl_spectra ,   norm_pnt):
    '''Ramanshift = vector, the x-axis in wavenumbers
       norm_pnt =  normalization point (corrections will be set
                                        to unity at this point) '''

    abs_wavenumber = ((1e7/laser_nm)-Ramanshift)


    init_guess=np.array([1e-18, 2799 ])
    # perform fit
    popt, pcov = curve_fit(photons_per_unit_wavenum_abs,
    abs_wavenumber, wl_spectra, p0=init_guess, bounds=([1e-28,1e-8], [900., 9000.]))

    print("\t Optimized coefs :", popt)

    # generate fit
    fit = photons_per_unit_wavenum_abs(abs_wavenumber, *popt)

    plt.plot(abs_wavenumber, wl_spectra,'o',abs_wavenumber, fit)
    #plt.plot(abs_wavenumber, wl_spectra,'o' )
    plt.grid()
    plt.ylim([0, 1.2])
    plt.title('Fit of broadband white light spectrum with black-body emission' )
    plt.xlabel('Wavenumber / cm$^{-1}$ (absolute)')
    plt.ylabel('Relative intensity')
    plt.show()

    #---------------------------
    C1 = wl_spectra / fit


    return C1

#############################################################################

def gen_C1_with_mask (Ramanshift, laser_nm ,  wl_spectra , mask ,  norm_pnt):

    '''Ramanshift = vector, the x-axis in wavenumbers
       norm_pnt =  normalization point (corrections will be set
                                        to unity at this point) '''

    abs_wavenumber = ((1e7/laser_nm)-Ramanshift)

    masked_wl  =   np.ma.masked_array(wl_spectra, mask=mask)
    masked_xaxis = np.ma.masked_array(abs_wavenumber, mask=mask)

    init_guess=np.array([1e-19, 2899 ])
    # perform fit
    popt, pcov = curve_fit(photons_per_unit_wavenum_abs,
    masked_xaxis, masked_wl, p0=init_guess,
    method='lm')

    print("\t Optimized coefs :", popt)

    # generate fit
    fit = photons_per_unit_wavenum_abs(abs_wavenumber, *popt)

    plt.plot(abs_wavenumber, masked_wl,'o',abs_wavenumber, fit,'--')
    #plt.plot(abs_wavenumber, wl_spectra,'o' )
    plt.grid()
    plt.ylim([0, 1.2])
    plt.title('Fit of broadband white light spectrum with black-body emission' )
    plt.xlabel('Wavenumber / cm$^{-1}$ (absolute)')
    plt.ylabel('Relative intensity')

    plt.show()

    #---------------------------
    C1 = wl_spectra / fit


    return C1
#############################################################################
