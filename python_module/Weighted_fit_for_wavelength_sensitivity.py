# -*- coding: utf-8 -*-

"""Module describing the weighted non-linear optimization scheme used to
determine the wavelength sensitivity of the spectrometer using a  polynomial
as a model function"""

# Load necessary modules
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

#********************************************************************

# Load data files

dataH2 = np.loadtxt("./dataH2.txt")
dataHD = np.loadtxt("./dataHD.txt")
dataD2 = np.loadtxt("./dataD2.txt")
dataO2 = np.loadtxt("./dataO2.txt")
dataO2_p = np.loadtxt("./dataO2_p.txt")

xaxis = np.loadtxt("./Ramanshift_axis.txt")

#print(len(dataH2), len(dataHD), len(dataD2), len(dataO2))

#********************************************************************

# Constants ------------------------------
scale1 = 1e4
scale2 = 1e7
scale3 = 1e9
# ----------------------------------------

def residual_quad(param):
    '''Function which computes the difference between the ratio of the 'curve'
    to the true given as 'ydata'  with that from  the function describing the
    perturbation. Here it is modeled as a quadratic function of
    the type, ( 1+ c1*x + c2*x**2 ) '''

    c1 = param[0]
    c2 = param[1]

	# - H2 -
    ratio_H2 = dataH2[:, 1]/dataH2[:, 2]
    RHS_H2 = (1.0 + c1/scale1 * dataH2[:, 3] + c2/scale2 * (dataH2[:, 3]**2))/ (1.0 +\
             c1/scale1 * dataH2[:, 4] + c2/scale2 * (dataH2[:, 4]**2))

    resd_H2 = dataH2[:, 5] * ((ratio_H2 - RHS_H2)**2)
    #print(resd_H2)
	# ------

	# - HD -
    ratio_HD = dataHD[:, 1]/dataHD[:, 2]
    RHS_HD = (1.0 + c1/scale1 * dataHD[:, 3] + c2/scale2 * (dataHD[:, 3]**2))/ (1.0 +\
             c1/scale1 * dataHD[:, 4] + c2/scale2 * (dataHD[:, 4]**2))

    resd_HD = dataHD[:, 5] * ((ratio_HD - RHS_HD)**2)
	# ------

	# - D2 -
    ratio_D2 = dataD2[:, 1]/dataD2[:, 2]
    RHS_D2 = (1.0 + c1/scale1 * dataD2[:, 3] + c2/scale2 * (dataD2[:, 3]**2))/ (1.0 +\
             c1/scale1 * dataD2[:, 4] + c2/scale2 * (dataD2[:, 4]**2))

    resd_D2 = dataD2[:, 5] * ((ratio_D2 - RHS_D2)**2)
	# ------

	# - O2 high frequency -
    ratio_O2 = dataO2[:, 1]/dataO2[:, 2]
    RHS_O2 = (1.0 + c1/scale1 * dataO2[:, 3] + c2/scale2 * (dataO2[:, 3]**2))/ (1.0 +\
             c1/scale1 * dataO2[:, 4] + c2/scale2 * (dataO2[:, 4]**2))

    resd_O2 = dataO2[:, 5] * ((ratio_O2 - RHS_O2)**2)
	# ------

    # - O2 pure rotation -
    ratio_O2p = dataO2_p[:, 1]/dataO2_p[:, 2]
    RHS_O2p = (1.0 + c1/scale1 * dataO2_p[:, 3] + c2/scale2 * (dataO2_p[:, 3]**2))/ (1.0 +\
             c1/scale1 * dataO2_p[:, 4] + c2/scale2 * (dataO2_p[:, 4]**2))

    resd_O2p = (dataO2_p[:, 5]*2.5) * ((ratio_O2p - RHS_O2p)**2)
	# ------


    resd = np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) + np.sum(resd_O2) \
    + np.sum(resd_O2p)
    return np.sum(resd)

#********************************************************************

def residual_cubic(param):
    '''Function which computes the difference between the ratio of the 'curve'
    to the true given as 'ydata'  with that from  the function describing the
    perturbation. Here it is modeled as a quadratic function of
    the type, ( 1+ c1*x + c2*x**2  + c3*x**3 ) '''

    c1 = param[0]
    c2 = param[1]
    c3 = param[2]

	# - H2 -
    ratio_H2 = dataH2[:, 1]/dataH2[:, 2]
    RHS_H2 = (1.0 + c1/scale1 * dataH2[:, 3] + c2/scale2 *(dataH2[:, 3]**2) +\
              c3/scale3 * (dataH2[:, 3]**3))/ ( 1.0 + c1/scale1 * dataH2[:, 4] +\
                          c2/scale2 * (dataH2[:, 4]**2) + c3/scale3 * ( dataH2[:, 4]**3))

    resd_H2 = dataH2[:, 5] * ((ratio_H2 - RHS_H2)**2)
    #print(resd_H2)
	# ------

	# - HD -
    ratio_HD = dataHD[:, 1]/dataHD[:, 2]
    RHS_HD = (1.0 + c1/scale1 * dataHD[:, 3] + c2/scale2 *(dataHD[:, 3]**2)+\
              c3/scale3 *(dataHD[:, 3]**3))/ ( 1.0 + c1/scale1 * dataHD[:, 4] +\
                         c2/scale2 * (dataHD[:, 4]**2)+ c3/scale3 * ( dataHD[:, 4]**3))

    resd_HD = dataHD[:, 5] * ((ratio_HD - RHS_HD)**2)
	# ------

	# - D2 -
    ratio_D2 = dataD2[:, 1]/dataD2[:, 2]
    RHS_D2 = (1.0 + c1/scale1 * dataD2[:, 3] + c2/scale2 * (dataD2[:, 3]**2)+\
              c3/scale3 * (dataD2[:, 3]**3))/ ( 1.0 + c1/scale1 * dataD2[:, 4] +\
                          c2/scale2 * (dataD2[:, 4]**2)+ c3/scale3 * ( dataD2[:, 4]**3))

    resd_D2 = dataD2[:, 5] * ((ratio_D2 - RHS_D2)**2)
	# ------

	# - O2 high frequency -
    ratio_O2 = dataO2[:, 1]/dataO2[:, 2]
    RHS_O2 = (1.0 + c1/scale1 * dataO2[:, 3] + c2/scale2 * (dataO2[:, 3]**2)+\
              c3/scale3 * (dataO2[:, 3]**3))/ ( 1.0 + c1/scale1 * dataO2[:, 4] +\
                          c2/scale2 *(dataO2[:, 4]**2)+ c3/scale3 *( dataO2[:, 4]**3))

    resd_O2 = dataO2[:, 5] * ((ratio_O2 - RHS_O2)**2)
	# ------

    # - O2 pure rotation -
    ratio_O2p = dataO2_p[:, 1]/dataO2_p[:, 2]
    RHS_O2p = (1.0 + c1/scale1 * dataO2_p[:, 3] + c2/scale2 * (dataO2_p[:, 3]**2)+\
               c3/scale3 * (dataO2_p[:, 3]**3))/ ( 1.0 + c1/scale1 * dataO2_p[:, 4] +\
                           c2/scale2 * (dataO2_p[:, 4]**2)+ c3/scale3 * ( dataO2_p[:, 4]**3))

    resd_O2p = (dataO2_p[:, 5]) * ((ratio_O2p - RHS_O2p)**2)
	# ------

    resd = np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) + np.sum(resd_O2) \
    + np.sum(resd_O2p)
    return np.sum(resd)

#***************************************************************

# Define the residual function

#***************************************************************

def run_fit_quad(init_k1, init_k2):
    # Intial guess

    param_init = np.array([ init_k1, init_k2 ])
    print("**********************************************************")
    #print("Testing the residual function with data")
    print("Initial coef :  k1={0}, k2={1} output = {2}".format(init_k1, init_k2,\
          (residual_quad(param_init))))


    print("\nOptimization run     \n")
    res = opt.minimize(residual_quad, param_init, method='Nelder-Mead', \
                              options={'xatol': 1e-9, 'fatol': 1e-9})

    print(res)
    optk1 = res.x[0]
    optk2 = res.x[1]
    print("\nOptimized result : k1={0}, k2={1} \n".format(round(optk1, 6), round(optk2, 6)))


    correction_curve_quad = 1+(optk1/scale1)*xaxis +(optk2/scale2)*(xaxis**2)  # generate the correction curve

    np.savetxt("correction_quad.txt", correction_curve_quad, fmt='%2.8f',\
               header='corrn_curve_quad', comments='')

    print("**********************************************************")

#***************************************************************
#***************************************************************

def run_fit_cubic(k1, k2, k3):
    # Intial guess
    init_k1 = k1
    init_k2 = k2
    init_k3 = k3

    param_init = np.array([init_k1, init_k2, init_k3 ])
    print("**********************************************************")
    print("Initial coef :  k1={0}, k2={1}, k3={2}, output = {3}".format(init_k1, init_k2,\
          init_k3, (residual_cubic(param_init))))

    print("\nOptimization run     \n")
    res = opt.minimize(residual_cubic, param_init, method='Nelder-Mead', \
                              options={'xatol': 1e-9, 'fatol': 1e-9})

    print(res)
    optk1 = res.x[0]
    optk2 = res.x[1]
    optk3 = res.x[2]
    print("\nOptimized result : k1={0}, k2={1}, k3={2}  ".format(round(optk1, 6),\
          round(optk2, 6), round(optk3, 6)))


    correction_curve_cub = 1+(optk1/scale1)*xaxis +(optk2/scale2)*(xaxis**2)+\
    (optk3/scale3)* (xaxis**3)  # generate the correction curve

    np.savetxt("correction_cubic.txt", correction_curve_cub, fmt='%2.8f',\
               header='corrn_curve_cubic', comments='')

    print("**********************************************************")

#***************************************************************

run_fit_quad(0.86, -3.05)
run_fit_cubic(0.76, -2.05, 0.050)

curve_quad = np.loadtxt("./correction_quad.txt", skiprows=1)
curve_cubic = np.loadtxt("./correction_cubic.txt", skiprows=1)

#********************************************************************

# Plotting the data

txt = ("*Generated from 'polDerivative.py' on the\
      \nGitHub Repository: RamanSpec_IntensityCalibration ")

# FIGURE 0 INITIALIZED

#plt.figure(0,figsize=(15, 8.5), dpi=300)


plt.figure(0)
ax0 = plt.axes()
plt.title('Fitting result', fontsize=23)
plt.plot(xaxis,  curve_quad, 'r', linewidth=3, label='quad_fit')
plt.plot(xaxis,  curve_cubic, 'b', linewidth=3, label='cubic_fit')
#plt.plot(xdata, curve, 'b', linewidth=2.0)
#plt.plot(xdata, correction_curve, 'k--', linewidth=1.5)

#plt.plot(xdata, noise, 'k', linewidth=1.0)

plt.xlabel('Ramanshift / $cm^{-1}$', fontsize=25)
plt.ylabel('Relative sensitivity', fontsize=25)
plt.grid(True)
ax0.tick_params(axis='both', labelsize =24)
plt.xlim((1755, -1120))
ax0.set_ylim([0, 1.25]) # change this if the ydata function or perturbation is changed

ax0.minorticks_on()
ax0.tick_params(which='minor', right='on')
ax0.tick_params(axis='y', labelleft='on', labelright='on')
plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey",\
         transform=plt.gcf().transFigure)
#plt.legend(loc='upper left', fontsize=24)

#plt.savefig('test_NLO.png', dpi=300)
#********************************************************************

