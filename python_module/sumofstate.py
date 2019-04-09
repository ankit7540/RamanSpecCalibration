#!/usr/bin/python

import math   # This will import math module
import numpy as np

# Constants ------------------------------
k = np.float64(1.38064852e-23) # J/K
h = np.float64(6.626070040e-34)  # J.s
c = np.float64(2.99792458e+10)   # cm/s

# ----------------------------------------

# Laser properties------------------------
omega = 18789.9500    #cm-1
# ----------------------------------------

# Load data on the energy levels and the polarizability anisotropy

eJH2v0 = np.loadtxt("./energy_levels_and_gamma/H2eV0.txt")
eJH2v1 = np.loadtxt("./energy_levels_and_gamma/H2eV1.txt")
eJHDv0 = np.loadtxt("./energy_levels_and_gamma/HDeV0.txt")
eJHDv1 = np.loadtxt("./energy_levels_and_gamma/HDeV1.txt")
eJD2v0 = np.loadtxt("./energy_levels_and_gamma/D2eV0.txt")
eJD2v1 = np.loadtxt("./energy_levels_and_gamma/D2eV1.txt")

ME_H2_532 = np.loadtxt("./energy_levels_and_gamma/ME_gamma_532.199323_H2.txt")
ME_HD_532 = np.loadtxt("./energy_levels_and_gamma/ME_gamma_532.199323_HD.txt")
ME_D2_532 = np.loadtxt("./energy_levels_and_gamma/ME_gamma_532.199323_D2.txt")

print(eJH2v0)
print(eJH2v1)
print(eJD2v0)
print(eJD2v1)
print(eJHDv0)
print(eJHDv1)


print(ME_H2_532)
print(ME_HD_532)
print(ME_D2_532)
# ----------------------------------------

# set of functions for computing the Sum of states of molecular hydrogen and
#   its isotopologues at given temperature.
# Data on energy levels is needed for the specific molecule


def sumofstate_H2(T):
    """calculate the sum of state for H2 molecule at T """

    Q = np.float64(0.0)

    #--- nuclear spin statistics ------------
    ge = 1 	# hydrogen
    go = 3
    # ---------------------------------------

    jmaxv0 = len(eJH2v0)
    jmaxv1 = len(eJH2v1)

#    print(jmaxv0, jmaxv1)

    # contribution from the ground vibrational state
    for i in range(0, jmaxv0):
        E = eJH2v0[i]
        energy = (-1*E*h*c)
        factor = (2*i+1)*math.exp(energy/(k*T))
        if i % 2 == 0:
            factor = factor*ge
        else:
            factor = factor*go
        Q = Q+factor

    # contribution from the first vibrational state
    for i in range(0, jmaxv1):
        E = eJH2v1[i]
        energy = (-1*E*h*c)
        factor = (2*i+1)*math.exp(energy/(k*T))
        if i % 2 == 0:
            factor = factor*ge
        else :
            factor = factor*go
        Q = Q+factor

#   return the output
    return Q

#********************************************************************

# compute the temperature dependent sum of state for D2 which includes contributions
# from the ground and first vibrational state of electronic ground state.

def sumofstate_D2(T):
    """calculate the sum of state for D2 molecule at T """

    Q = np.float64(0.0)

    #--- nuclear spin statistics ------------
    ge = 6        # deuterium
    go = 3
    # ---------------------------------------

    jmaxv0 = len(eJD2v0)
    jmaxv1 = len(eJD2v1)
#    print(jmaxv0, jmaxv1)

    # contribution from the ground vibrational state
    for i in range(0, jmaxv0):
        E = eJD2v0[i]
        energy = (-1*E*h*c)
        factor = (2*i+1)*math.exp(energy/(k*T))
        if i % 2 == 0:
            factor = factor*ge
        else :
            factor = factor*go
        Q = Q+factor

    # contribution from the first vibrational state
    for i in range(0, jmaxv1):
        E = eJD2v1[i]
        energy = (-1*E*h*c)
        factor = (2*i+1)*math.exp(energy/(k*T))
        if i % 2 == 0:
            factor = factor*ge
        else :
            factor = factor*go
        Q = Q+factor

    # return the output
    return Q

#********************************************************************

# compute the temperature dependent sum of state for HD which includes contributions
# from the ground and first vibrational state of electronic ground state.

def sumofstate_HD(T):
    """calculate the sum of state for HD molecule at T """

    Q = np.float64(0.0)

    #--- nuclear spin statistics ------------
    ge = 1        # hydrogen deuteride
    go = 1
    # ---------------------------------------

    jmaxv0 = len(eJHDv0)
    jmaxv1 = len(eJHDv1)
#    print(jmaxv0, jmaxv1)

    # contribution from the ground vibrational state
    for i in range(0, jmaxv0):
        E = eJHDv0[i]
        energy = (-1*E*h*c)
        factor = (2*i+1)*math.exp(energy/(k*T))
        if i % 2 == 0:
            factor = factor*ge
        else :
            factor = factor*go
        Q = Q+factor

    # contribution from the first vibrational state

    for i in range(0,jmaxv1):
        E = eJHDv1[i]
        energy = (-1*E*h*c)
        factor = (2*i+1)*math.exp(energy/(k*T))
        if i % 2 == 0:
            factor = factor*ge
        else :
            factor = factor*go
        Q = Q+factor

    # return the output
    return Q

#********************************************************************

def spectra_H2(T, Js, Jas):
    """Compute in intensities and position for rotational Raman bands of H2 """

    #--- nuclear spin statistics ------------
    ge = 1        # hydrogen
    go = 3
    # ---------------------------------------

    sos=sumofstate_H2(T)
    print(sos)

    # define arrays for intensity, position and J_{i}
    specH2 = np.zeros(shape=(Js+Jas, 4))
    posnH2 = np.zeros(shape=(Js+Jas, 1))
    spectraH2 = np.zeros(shape=(Js+Jas, 1))

    print("Stokes line")
    # Stokes lines
    for i in range(0, Js+1):
        E = eJH2v0[i]
        energy = (-1*E*h*c)
        popn = (2*i+1)*math.exp(energy/(k*T))
        bj = (3*(i+1)*(i+2))/(2*(2*i+1)*(2*i+3))
        position = (eJH2v0[i+2]-eJH2v0[i])
        gamma = ME_H2_532[i+1][2]
        print(gamma, ME_H2_532[i+1][0], ME_H2_532[i+1][1] )
        if i % 2 == 0:
            factor = popn*ge*bj*omega*((omega-position)**3) *(gamma**2)/sos
        else :
            factor = popn*go*bj* omega*((omega-position)**3) *(gamma**2)/sos

        specH2[(Jas-1)+i][0] = i
        specH2[(Jas-1)+i][1] = position
        specH2[(Jas-1)+i][2] = factor
        specH2[(Jas-1)+i][3] = (omega-position)

        posnH2[(Jas-1)+i] = position
        spectraH2[(Jas-1)+i] = factor

    # Anti-Stokes lines
    i = Jas
    print("AS")
    for i in range(Jas, 1, -1):
        E = eJH2v0[i]
        energy = (-1*E*h*c)
        gamma = ME_H2_532[i-1][2]
        print(gamma)
        #print(i,E, i-1,gamma)
        popn = (2*i+1)*math.exp(energy/(k*T))
        bj = (3*(i)*(i-1))/(2*(2*i-1)*(2*i+1))
        position = -1*( eJH2v0[i] - eJH2v0[i-2])
        if i % 2 == 0:
            factor = popn*ge*bj*omega*((omega-position)**3) *(gamma**2)/sos
        else :
            factor = popn*go*bj*omega*((omega-position)**3) *(gamma**2)/sos

        specH2[Jas-i][0] = i
        specH2[Jas-i][1] = position
        specH2[Jas-i][2] = factor
        specH2[Jas-i][3] = (omega-position)
        posnH2[Jas-i] = position
        spectraH2[Jas-i] = factor

#   print(sum)

#   return the output

    normalize1d(spectraH2)
    print("Normalized spectra H2 :\n", spectraH2)
    np.savetxt('posnH2.txt',(posnH2), fmt='%+5.12f') #save the band position.
    np.savetxt('spectraH2.txt',(spectraH2), fmt='%+5.12f') #save the band intensity

#    return factor

#********************************************************************

def normalize1d(array_name):
    maxV = np.max(array_name)
    size = len(array_name)
    for i in range(0, size,1):
        array_name[i] = array_name[i]/maxV

#********************************************************************
#********************************************************************

def spectra_D2(T, Js, Jas):
    """Compute in intensities and position for rotational Raman bands of D2 """

    #---MOLECULE PROPERTIES------------
    ge = 2        # deuterium
    go = 1
    # ---------------------------------

    sos=sumofstate_D2(T)

    #print(sos)

    # define arrays for intensity, position and J_{i}
    specD2 = np.zeros(shape=(Js+Jas, 4))
    posnD2 = np.zeros(shape=(Js+Jas, 1))
    spectraD2 = np.zeros(shape=(Js+Jas, 1))
#    print(jmaxv0, jmaxv1)

    # Stokes lines

    for i in range(0, Js+1):
        E = eJD2v0[i]
        energy = (-1*E*h*c)
        popn = (2*i+1)*math.exp(energy/(k*T))
        bj = (3*(i+1)*(i+2))/(2*(2*i+1)*(2*i+3))
        position = (eJD2v0[i+2]-eJD2v0[i])
        gamma = ME_D2_532[i+1][2]
#        print(i,E,gamma)

        if i % 2 == 0:
            factor = popn*ge*bj*omega*((omega-position)**3) *(gamma**2)/sos
        else :
            factor = popn*go*bj* omega*((omega-position)**3) *(gamma**2)/sos

#        print (round(popn,4),"\t",round(bj,3),"\t",factor,"\t", E,"\t", math.exp(energy/(k*T)) )
        specD2[(Jas-1)+i][0] = i
        specD2[(Jas-1)+i][1] = position
        specD2[(Jas-1)+i][2] = factor
        specD2[(Jas-1)+i][3] = omega-position
        posnD2[(Jas-1)+i] = position
        spectraD2[(Jas-1)+i] = factor


    # Anti-Stokes lines
    i = Jas
    for i in range(Jas,1,-1):
        E = eJD2v0[i]
        energy = (-1*E*h*c)
        gamma = ME_D2_532[i-1][2]
#        print(i,E, i-1,gamma)
        popn = (2*i+1)*math.exp(energy/(k*T))
        bj = (3*(i)*(i-1))/(2*(2*i-1)*(2*i+1))
        position = -1*(eJD2v0[i]-eJD2v0[i-2])
        #print("Position,i,E,q1,bj,gamma:",position,i,E,math.exp(energy/(k*T)),bj,gamma)

        if i % 2 ==0:
            factor = popn*ge*bj*omega*((omega-position)**3) *(gamma**2)/sos
        else :
            factor = popn*go*bj*omega*((omega-position)**3) *(gamma**2)/sos
#        print (round(popn,4),"\t",round(bj,3),"\t",factor,"\t", E,"\t", math.exp(energy/(k*T)) )
        specD2[Jas-i][0] = i
        specD2[Jas-i][1] = position
        specD2[Jas-i][2] = factor
        specD2[Jas-i][3] = omega-position
        posnD2[Jas-i] = position
        spectraD2[Jas-i] = factor

#   return the output
#    print(specD2)
    normalize1d(spectraD2)
    print("Normalized spectra D2 :\n",spectraD2)

    np.savetxt('posnD2.txt',(posnD2), fmt='%+5.12f') #save the band position.
    np.savetxt('spectraD2.txt',(spectraD2), fmt='%+5.12f') #save the band intensity

#    return factor

#********************************************************************
#********************************************************************

def spectra_HD(T, Js, Jas):
    """Compute in intensities and position for rotational Raman bands of HD """

    #---MOLECULE PROPERTIES------------
    ge=1        # deuterium hydride (HD)
    go=1
    # ---------------------------------

    sos=sumofstate_HD(T)
    #print(sos)

    # define arrays for intensity, position and J_{i}
    specHD = np.zeros(shape=(Js+Jas, 4))
    posnHD = np.zeros(shape=(Js+Jas, 1))
    spectraHD = np.zeros(shape=(Js+Jas, 1))

#    print(jmaxv0, jmaxv1)
    # Stokes lines
    for i in range(0,Js+1):
        E = eJHDv0[i]
        energy = (-1*E*h*c)
        popn = (2*i+1)*math.exp(energy/(k*T))
        bj = (3*(i+1)*(i+2))/(2*(2*i+1)*(2*i+3))
        position = (eJHDv0[i+2]-eJHDv0[i])

        gamma = ME_HD_532[i+1][2]
#        print("Position,i,E,q1,bj,gamma:",position,i,E,math.exp(energy/(k*T)),bj,gamma)
        if i % 2 ==0:
            factor = popn*ge*bj*omega*((omega-position)**3) *(gamma**2)/sos
        else :
            factor = popn*go*bj* omega*((omega-position)**3) *(gamma**2)/sos
        #print (popn,bj,factor, E, math.exp(energy/(k*T)) )
        specHD[(Jas-1)+i][0] = i
        specHD[(Jas-1)+i][1] = position
        specHD[(Jas-1)+i][2] = factor
        specHD[(Jas-1)+i][3] = omega-position
        posnHD[(Jas-1)+i] = position
        spectraHD[(Jas-1)+i] = factor
#        print("Index on result:",(Jas-1+i))


    # Anti-Stokes lines
    i = Jas
    for i in range(Jas, 1, -1):
        E = eJHDv0[i]
        energy = (-1*E*h*c)
        gamma = ME_HD_532[i-1][2]
        popn = (2*i+1)*math.exp(energy/(k*T))
        bj = (3*(i)*(i-1))/(2*(2*i-1)*(2*i+1))
        position = -1*(eJHDv0[i]-eJHDv0[i-2])
#        print("Position,i,E,q1,bj,gamma:",position,i,E,math.exp(energy/(k*T)),bj,gamma)
        if i % 2 == 0:
            factor = popn*ge*bj*omega*((omega-position)**3) *(gamma**2)/sos
        else :
            factor = popn*go*bj*omega*((omega-position)**3) *(gamma**2)/sos
        # print (popn,bj,factor, E, math.exp(energy/(k*T)) )
        specHD[Jas-i][0] = i
        specHD[Jas-i][1] = position
        specHD[Jas-i][2] = factor
        specHD[Jas-i][3] = omega-position
        posnHD[Jas-i] = position
        spectraHD[Jas-i] = factor
#        print("Index on result:",(Jas-i))


#   return the output
#    print("unnormalized spectra:",spectraHD,specHD)
    normalize1d(spectraHD)
    print("Normalized spectra HD: \n", spectraHD)

    np.savetxt('posnHD.txt', (posnHD), fmt='%+5.12f') #save the band position.
    np.savetxt('spectraHD.txt', (spectraHD), fmt='%+5.12f') #save the band intensity

#    return factor

#********************************************************************

print(sumofstate_H2(333))

print(sumofstate_D2(333))

print(sumofstate_HD(333))

spectra_H2(578, 6, 6)
spectra_D2(578, 6, 8)
spectra_HD(578, 7, 7)


H2posn=np.loadtxt("posnH2.txt")
D2posn=np.loadtxt("posnD2.txt")
HDposn=np.loadtxt("posnHD.txt")

H2spectra=np.loadtxt("spectraH2.txt")
D2spectra=np.loadtxt("spectraD2.txt")
HDspectra=np.loadtxt("spectraHD.txt")


#print(D2spectra,HDspectra,H2spectra)

#x3=[0]
#y3=[0]
#plotter2d( H2posn, H2spectra, D2posn, D2spectra, HDposn, HDspectra,'bo-','ro-','go-',1,400)
