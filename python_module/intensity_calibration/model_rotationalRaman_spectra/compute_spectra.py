#!/usr/bin/python
'''Module for computing the pure rotational Raman spectra from H2, HD and D2'''

import math
import numpy as np
import matplotlib.pyplot as plt

# Constants ------------------------------
K = np.float64(1.38064852e-23) # J/K
H = np.float64(6.626070040e-34)  # J.s
C = np.float64(2.99792458e+10)   # cm/s
# ----------------------------------------

# Laser properties------------------------
omega = 18789.9850    #absolute cm-1
omega_sc = omega/1e4  # scaled frequency (for better numerical accuracy)
# ----------------------------------------

# Load data on the energy levels and the polarizability anisotropy

# Data on the rovibrational energy levels has been extracted from the calculated dissociation
# energy data published in the following works :
#   a)  J. Komasa, K. Piszczatowski, G.  Lach, M. Przybytek, B. Jeziorski,
#        and K. Pachucki, J. Chem. Theory and Comput. 7, 3105 (2011).
#
#   b) K. Pachucki and J. Komasa, Phys. Chem. Chem. Phys. 12, 9188 (2010).

eJH2v0 = np.loadtxt("./energy_levels_and_gamma/H2eV0.txt")
eJH2v1 = np.loadtxt("./energy_levels_and_gamma/H2eV1.txt")
eJHDv0 = np.loadtxt("./energy_levels_and_gamma/HDeV0.txt")
eJHDv1 = np.loadtxt("./energy_levels_and_gamma/HDeV1.txt")
eJD2v0 = np.loadtxt("./energy_levels_and_gamma/D2eV0.txt")
eJD2v1 = np.loadtxt("./energy_levels_and_gamma/D2eV1.txt")

#   Data on the matrix elements of polarizability anisotropy has been taken from
#   our previous work.
#   c) A. Raj, H. Hamaguchi, and H. A. Witek, J. Chem. Phys. 148, 104308 (2018).

ME_H2_532 = np.loadtxt("./energy_levels_and_gamma/ME_gamma_532.199323_H2.txt")
ME_HD_532 = np.loadtxt("./energy_levels_and_gamma/ME_gamma_532.199323_HD.txt")
ME_D2_532 = np.loadtxt("./energy_levels_and_gamma/ME_gamma_532.199323_D2.txt")

#print(eJH2v0)
#print(eJH2v1)
#print(eJD2v0)
#print(eJD2v1)
#print(eJHDv0)
#print(eJHDv1)

#print(ME_H2_532)
#print(ME_HD_532)
#print(ME_D2_532)

#********************************************************************
#********************************************************************

# set of functions for computing the Sum of states of molecular hydrogen and
#   its isotopologues at given temperature.
# Data on energy levels is needed for the specific molecule

def sumofstate_H2(T):
    """calculate the sum of state for H2 molecule at T """

    Q = np.float64(0.0)

    #--- nuclear spin statistics ------------
    g_even = 1 	# hydrogen
    g_odd = 3
    # ---------------------------------------

    jmaxv0 = len(eJH2v0)
    jmaxv1 = len(eJH2v1)

    # contribution from the ground vibrational state
    for i in range(0, jmaxv0):
        E = eJH2v0[i]
        energy = (-1*E*H*C)
        factor = (2*i+1)*math.exp(energy/(K*T))
        if i % 2 == 0:
            factor = factor*g_even
        else:
            factor = factor*g_odd
        Q = Q+factor

    # contribution from the first vibrational state
    for i in range(0, jmaxv1):
        E = eJH2v1[i]
        energy = (-1*E*H*C)
        factor = (2*i+1)*math.exp(energy/(K*T))
        if i % 2 == 0:
            factor = factor*g_even
        else:
            factor = factor*g_odd
        Q = Q+factor

#   return the sum of states for H2
    return Q

#********************************************************************
#********************************************************************

# compute the temperature dependent sum of state for HD which includes contributions
# from the ground and first vibrational state of electronic ground state.

def sumofstate_HD(T):
    """calculate the sum of state for HD molecule at T """

    Q = np.float64(0.0)

    #--- nuclear spin statistics ------------
    g_even = 1        # hydrogen deuteride
    g_odd = 1
    # ---------------------------------------

    jmaxv0 = len(eJHDv0)
    jmaxv1 = len(eJHDv1)

    # contribution from the ground vibrational state
    for i in range(0, jmaxv0):
        E = eJHDv0[i]
        energy = (-1*E*H*C)
        factor = (2*i+1)*math.exp(energy/(K*T))
        if i % 2 == 0:
            factor = factor*g_even
        else:
            factor = factor*g_odd
        Q = Q+factor

    # contribution from the first vibrational state
    for i in range(0, jmaxv1):
        E = eJHDv1[i]
        energy = (-1*E*H*C)
        factor = (2*i+1)*math.exp(energy/(K*T))
        if i % 2 == 0:
            factor = factor*g_even
        else:
            factor = factor*g_odd
        Q = Q+factor

#   return the sum of states for HD
    return Q

#********************************************************************
#********************************************************************

# compute the temperature dependent sum of state for D2 which includes contributions
# from the ground and first vibrational state of electronic ground state.

def sumofstate_D2(T):
    """calculate the sum of state for D2 molecule at T """

    Q = np.float64(0.0)

    #--- nuclear spin statistics ------------
    g_even = 6        # deuterium
    g_odd = 3
    # ---------------------------------------

    jmaxv0 = len(eJD2v0)
    jmaxv1 = len(eJD2v1)

    # contribution from the ground vibrational state
    for i in range(0, jmaxv0):
        E = eJD2v0[i]
        energy = (-1*E*H*C)
        factor = (2*i+1)*math.exp(energy/(K*T))
        if i % 2 == 0:
            factor = factor*g_even
        else:
            factor = factor*g_odd
        Q = Q+factor

    # contribution from the first vibrational state
    for i in range(0, jmaxv1):
        E = eJD2v1[i]
        energy = (-1*E*H*C)
        factor = (2*i+1)*math.exp(energy/(K*T))
        if i % 2 == 0:
            factor = factor*g_even
        else:
            factor = factor*g_odd
        Q = Q+factor

#   return the sum of states for D2
    return Q

#********************************************************************
#********************************************************************

def normalize1d(array_name):
    """Normalize a 1D array using the max value"""
    max_val = np.max(array_name)
    size = len(array_name)
    for i in range(0, size, 1):
        array_name[i] = array_name[i]/max_val

#********************************************************************
#********************************************************************

def spectra_H2(T, Js, Jas):
    """Compute in intensities and position for rotational Raman bands of H2 """

    #--- nuclear spin statistics ------------
    g_even = 1        # hydrogen
    g_odd = 3
    # ---------------------------------------

    sos = sumofstate_H2(T)

    # define arrays for intensity, position and J_{i}
    specH2 = np.zeros(shape=(Js+Jas, 4))
    posnH2 = np.zeros(shape=(Js+Jas, 1))
    spectraH2 = np.zeros(shape=(Js+Jas, 1))

    # Stokes lines
    for i in range(0, Js+1):
        E = eJH2v0[i]
        energy = (-1*E*H*C)
        popn = (2*i+1)*math.exp(energy/(K*T))
        bj = (3*(i+1)*(i+2))/(2*(2*i+1)*(2*i+3))
        position = (eJH2v0[i+2]-eJH2v0[i])
        gamma = ME_H2_532[i+1][2]

        if i % 2 == 0:
            factor = popn*g_even*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos
        else:
            factor = popn*g_odd*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos

        specH2[(Jas-1)+i][0] = i
        specH2[(Jas-1)+i][1] = position
        specH2[(Jas-1)+i][2] = factor  # unnormalized intensity, arbitrary unit
        specH2[(Jas-1)+i][3] = (omega-position)

        posnH2[(Jas-1)+i] = position
        spectraH2[(Jas-1)+i] = factor

    # Anti-Stokes lines
    i = Jas
    for i in range(Jas, 1, -1):
        E = eJH2v0[i]
        energy = (-1*E*H*C)
        popn = (2*i+1)*math.exp(energy/(K*T))
        bj = (3*(i)*(i-1))/(2*(2*i-1)*(2*i+1))
        position = -1*(eJH2v0[i] - eJH2v0[i-2])
        gamma = ME_H2_532[i-1][2]

        if i % 2 == 0:
            factor = popn*g_even*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos
        else:
            factor = popn*g_odd*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos

        specH2[Jas-i][0] = i
        specH2[Jas-i][1] = position
        specH2[Jas-i][2] = factor  # unnormalized intensity, arbitrary unit
        specH2[Jas-i][3] = (omega-position)

        posnH2[Jas-i] = position
        spectraH2[Jas-i] = factor

    # return the output

    normalize1d(spectraH2)
    specH2[:, 2] = spectraH2[:, 0]
    print("Normalized spectra H2: \n", spectraH2)
    #print(specH2)  # assign normalized intensity
    np.savetxt('posnH2.txt', (posnH2), fmt='%+6.6f') #save the band position.
    np.savetxt('spectraH2.txt', (spectraH2), fmt='%+5.11f') #save the band intensity
    np.savetxt('specH2.txt', (specH2), fmt='%+5.12f') #save the 2D array which has data on
    #  initial J number, band position, band intensity and the position in abs. wavenumbers

#********************************************************************
#********************************************************************

def spectra_HD(T, Js, Jas):
    """Compute in intensities and position for rotational Raman bands of HD """

    #--- nuclear spin statistics ------------
    g_even = 1        # hydrogen deuteride
    g_odd = 1
    # ---------------------------------------

    sos = sumofstate_HD(T)

    # define arrays for intensity, position and J_{i}
    specHD = np.zeros(shape=(Js+Jas, 4))
    posnHD = np.zeros(shape=(Js+Jas, 1))
    spectraHD = np.zeros(shape=(Js+Jas, 1))

    # Stokes lines
    for i in range(0, Js+1):
        E = eJHDv0[i]
        energy = (-1*E*H*C)
        popn = (2*i+1)*math.exp(energy/(K*T))
        bj = (3*(i+1)*(i+2))/(2*(2*i+1)*(2*i+3))
        position = (eJHDv0[i+2]-eJHDv0[i])
        gamma = ME_HD_532[i+1][2]

        if i % 2 == 0:
            factor = popn*g_even*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos
        else:
            factor = popn*g_odd*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos

        specHD[(Jas-1)+i][0] = i
        specHD[(Jas-1)+i][1] = position
        specHD[(Jas-1)+i][2] = factor  # unnormalized intensity, arbitrary unit
        specHD[(Jas-1)+i][3] = omega-position
        posnHD[(Jas-1)+i] = position
        spectraHD[(Jas-1)+i] = factor


    # Anti-Stokes lines
    i = Jas
    for i in range(Jas, 1, -1):
        E = eJHDv0[i]
        energy = (-1*E*H*C)
        popn = (2*i+1)*math.exp(energy/(K*T))
        bj = (3*(i)*(i-1))/(2*(2*i-1)*(2*i+1))
        position = -1*(eJHDv0[i]-eJHDv0[i-2])
        gamma = ME_HD_532[i-1][2]

        if i % 2 == 0:
            factor = popn*g_even*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos
        else:
            factor = popn*g_odd*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos

        specHD[Jas-i][0] = i
        specHD[Jas-i][1] = position
        specHD[Jas-i][2] = factor  # unnormalized intensity, arbitrary unit
        specHD[Jas-i][3] = omega-position
        posnHD[Jas-i] = position
        spectraHD[Jas-i] = factor


#   return the output
    normalize1d(spectraHD)
    specHD[:, 2] = spectraHD[:, 0]  # assign normalized intensity
    print("Normalized spectra HD: \n", spectraHD)

    np.savetxt('posnHD.txt', (posnHD), fmt='%+6.6f') #save the band position.
    np.savetxt('spectraHD.txt', (spectraHD), fmt='%+5.11f') #save the band intensity
    np.savetxt('specHD.txt', (specHD), fmt='%+5.12f') #save the 2D array which has data on
    #  initial J number, band position, band intensity and the position in abs. wavenumbers

#********************************************************************
#********************************************************************

def spectra_D2(T, Js, Jas):
    """Compute in intensities and position for rotational Raman bands of D2 """

    #--- nuclear spin statistics ------------
    g_even = 6        # deuterium
    g_odd = 3
    # ---------------------------------------

    sos = sumofstate_D2(T)

    # define arrays for intensity, position and J_{i}
    specD2 = np.zeros(shape=(Js+Jas, 4))
    posnD2 = np.zeros(shape=(Js+Jas, 1))
    spectraD2 = np.zeros(shape=(Js+Jas, 1))

    # Stokes lines
    for i in range(0, Js+1):
        E = eJD2v0[i]
        energy = (-1*E*H*C)
        popn = (2*i+1)*math.exp(energy/(K*T))
        bj = (3*(i+1)*(i+2))/(2*(2*i+1)*(2*i+3))
        position = (eJD2v0[i+2]-eJD2v0[i])
        gamma = ME_D2_532[i+1][2]

        if i % 2 == 0:
            factor = popn*g_even*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos
        else:
            factor = popn*g_odd*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos

        specD2[(Jas-1)+i][0] = i
        specD2[(Jas-1)+i][1] = position
        specD2[(Jas-1)+i][2] = factor  # unnormalized intensity, arbitrary unit
        specD2[(Jas-1)+i][3] = omega-position
        posnD2[(Jas-1)+i] = position
        spectraD2[(Jas-1)+i] = factor


    # Anti-Stokes lines
    i = Jas
    for i in range(Jas, 1, -1):
        E = eJD2v0[i]
        energy = (-1*E*H*C)
        popn = (2*i+1)*math.exp(energy/(K*T))
        bj = (3*(i)*(i-1))/(2*(2*i-1)*(2*i+1))
        position = -1*(eJD2v0[i]-eJD2v0[i-2])
        gamma = ME_D2_532[i-1][2]

        if i % 2 == 0:
            factor = popn*g_even*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos
        else:
            factor = popn*g_odd*bj*omega_sc*((omega_sc-position/1e4)**3) *(gamma**2)/sos

        specD2[Jas-i][0] = i
        specD2[Jas-i][1] = position
        specD2[Jas-i][2] = factor  # unnormalized intensity, arbitrary unit
        specD2[Jas-i][3] = omega-position
        posnD2[Jas-i] = position
        spectraD2[Jas-i] = factor

#   return the output
    normalize1d(spectraD2)
    specD2[:, 2] = spectraD2[:, 0]  # assign normalized intensity
    print("Normalized spectra D2 :\n", spectraD2)

    np.savetxt('posnD2.txt', posnD2, fmt='%+6.6f') #save the band position.
    np.savetxt('spectraD2.txt', spectraD2, fmt='%+5.11f') #save the band intensity
    np.savetxt('specD2.txt', (specD2), fmt='%+5.12f') #save the 2D array which has data on
    #  initial J number, band position, band intensity and the position in abs. wavenumbers

#********************************************************************
#********************************************************************

print(" Sum of state for  H2 at 333 K : ", sumofstate_H2(333))

print(" Sum of state for  HD at 375 K : ", sumofstate_HD(375))

print(" Sum of state for  D2 at 298 K : ", sumofstate_D2(298))

#********************************************************************
print ("\n")
spectra_H2(298, 6, 6)
print ("\n")
spectra_HD(298, 7, 7)
print ("\n")
spectra_D2(298, 6, 8)

#  Spectra can  be plotted using matplotlib simply using the band posn and the spectra

#********************************************************************

# Plotting the data :  MATPLOTLIB REQUIRED

txt = ("*Generated from 'compute_spectra.py' on the\
      \nGitHub Repository: RamanSpecCalibration ")

# FIGURE 0 INITIALIZED

wavenumH2= np.loadtxt("./posnH2.txt")
wavenumHD= np.loadtxt("./posnHD.txt")
wavenumD2= np.loadtxt("./posnD2.txt")

spectraH2= np.loadtxt("./spectraH2.txt")
spectraHD= np.loadtxt("./spectraHD.txt")
spectraD2= np.loadtxt("./spectraD2.txt")

plt.figure(0)
ax0 = plt.axes()
plt.title('Calculated  pure rotational Raman spectra', fontsize=17)

plt.stem( wavenumH2,  spectraH2, 'r', label='$H_2$' ,  markerfmt='ro' , basefmt='k-')
plt.stem( wavenumHD,  spectraHD, 'g', label='$HD$',  markerfmt='go', basefmt='k-')
plt.stem( wavenumD2,  spectraD2, 'b', label='$D_2$',  markerfmt='bo', basefmt='k-')


#plt.plot( wavenumH2,  spectraH2, 'r|',  label='wavenumber (ref)')
#plt.plot( wavenumHD,  spectraHD, 'go',  label='wavenumber (ref)')
#plt.plot( wavenumD2,  spectraD2, 'bo',  label='wavenumber (ref)')

plt.xlabel('Wavenumber / cm-1', fontsize=16)
plt.ylabel('Relative intensity', fontsize=16)
plt.grid(True)
ax0.tick_params(axis='both', labelsize =15)

ax0.set_xlim([1700, -1065])

ax0.minorticks_on()
ax0.tick_params(which='minor', right='on')
ax0.tick_params(axis='y', labelleft='on', labelright='on')
plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey",\
         transform=plt.gcf().transFigure)
plt.legend(loc='upper left', fontsize=14)

plt.savefig('spectra.png', bbox_inches='tight', dpi=300)

#********************************************************************