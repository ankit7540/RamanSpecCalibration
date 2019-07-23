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
# Modify this as user case
dataH2 = np.loadtxt("./dataH2.txt")
dataHD = np.loadtxt("./dataHD.txt")
dataD2 = np.loadtxt("./dataD2.txt")
dataO2 = np.loadtxt("./dataO2.txt")
dataO2_p = np.loadtxt("./dataO2_p.txt")
xaxis = np.loadtxt("./Ramanshift_axis_g_para2_ED.txt")

#********************************************************************

# Constants ------------------------------
scale1 = 1e4
scale2 = 1e7
scale3 = 1e9
# ----------------------------------------

#*******************************************************************
# Define the residual function
#*******************************************************************

def residual_linear(param):
    '''Function which computes the residual (as sum of squares) comparing the
    ratio of expt to theoretical intensity ratio to the sensitivity  profile
    modelled as  a line, ( 1+ c1*x ) '''

    c1 = param[0]

	# - H2 -
    ratio_H2 = dataH2[:, 1]/dataH2[:, 2]
    RHS_H2 = (1.0 + c1/scale1 * dataH2[:, 3] )/ (1.0 +\
             c1/scale1 * dataH2[:, 4] )

    resd_H2 = dataH2[:, 5] * ((ratio_H2 - RHS_H2)**2)
    #print(resd_H2)
	# ------

	# - HD -
    ratio_HD = dataHD[:, 1]/dataHD[:, 2]
    RHS_HD = (1.0 + c1/scale1 * dataHD[:, 3] )/ (1.0 +\
             c1/scale1 * dataHD[:, 4] )

    resd_HD = dataHD[:, 5] * ((ratio_HD - RHS_HD)**2)
	# ------

	# - D2 -
    ratio_D2 = dataD2[:, 1]/dataD2[:, 2]
    RHS_D2 = (1.0 + c1/scale1 * dataD2[:, 3] )/ (1.0 +\
             c1/scale1 * dataD2[:, 4] )

    resd_D2 = dataD2[:, 5] * ((ratio_D2 - RHS_D2)**2)
	# ------

	# - O2 high frequency -
    ratio_O2 = dataO2[:, 1]/dataO2[:, 2]
    RHS_O2 = (1.0 + c1/scale1 * dataO2[:, 3] )/ (1.0 +\
             c1/scale1 * dataO2[:, 4] )

    resd_O2 = dataO2[:, 5] * ((ratio_O2 - RHS_O2)**2)
	# ------

    # - O2 pure rotation -
    ratio_O2p = dataO2_p[:, 1]/dataO2_p[:, 2]
    RHS_O2p = (1.0 + c1/scale1 * dataO2_p[:, 3] )/ (1.0 +\
             c1/scale1 * dataO2_p[:, 4] )

    resd_O2p = (dataO2_p[:, 5]*2.5) * ((ratio_O2p - RHS_O2p)**2)
	# ------


    return np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) + \
    np.sum(resd_O2)  + np.sum(resd_O2p)

#***************************************************************

def residual_quadratic(param):
    '''Function which computes the residual (as sum of squares) comparing the
    ratio of expt to theoretical intensity ratio to the sensitivity  profile
    modelled as  a quadratic polynomial, ( 1 + c1*x + c2*x**2 ) '''

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


    return np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) +\
     np.sum(resd_O2)  + np.sum(resd_O2p)

#********************************************************************

def residual_cubic(param):
    '''Function which computes the residual (as sum of squares) comparing the
    ratio of expt to theoretical intensity ratio to the sensitivity  profile
    modelled as  a quadratic polynomial, ( 1 + c1*x + c2*x**2 + c3*x**3 ) '''

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

    return np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) +\
     np.sum(resd_O2)  + np.sum(resd_O2p)

#***************************************************************
#***************************************************************
# Fit functions
#***************************************************************
#***************************************************************

def run_fit_linear (init_k1 ):
    '''Function performing the actual fit using the residual_linear function
    defined earlier '''

    # init_k1 : Intial guess

    param_init = np.array([ init_k1  ])
    print("**********************************************************")
    #print("Testing the residual function with data")
    print("Initial coef :  k1={0}  output = {1}".format(init_k1, \
          (residual_linear(param_init))))


    print("\nOptimization run     \n")
    res = opt.minimize(residual_linear, param_init, method='Nelder-Mead', \
                              options={'xatol': 1e-9, 'fatol': 1e-9})

    print(res)
    optk1 = res.x[0]
    print("\nOptimized result : k1={0} \n".format(round(optk1, 6) ) )

    correction_curve_line= 1+(optk1/scale1)*xaxis     # generate the correction curve

    np.savetxt("correction_linear.txt", correction_curve_line, fmt='%2.8f',\
               header='corrn_curve_linear', comments='')

    print("**********************************************************")

#***************************************************************

def run_fit_quadratic(init_k1, init_k2):
    '''Function performing the actual fit using the residual_quadratic function
    defined earlier '''

    #  init_k1  =   initial  guess  for k1
    #  init_k2  =   initial  guess  for k2

    param_init = np.array([ init_k1, init_k2 ])
    print("**********************************************************")
    #print("Testing the residual function with data")
    print("Initial coef :  k1={0}, k2={1} output = {2}".format(init_k1, init_k2,\
          (residual_quadratic(param_init))))

    print("\nOptimization run     \n")
    res = opt.minimize(residual_quadratic, param_init, method='Nelder-Mead', \
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

def run_fit_cubic(init_k1, init_k2, init_k3):
    '''Function performing the actual fit using the residual_cubic function
    defined earlier '''

    #  init_k1  =   initial  guess  for k1
    #  init_k2  =   initial  guess  for k2
    #  init_k3  =   initial  guess  for k3

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

# Run the fit with initial guesses
run_fit_linear(0.86 )
run_fit_quadratic(0.86, -3.05)
run_fit_cubic(0.76, -2.05, 0.050)

# Load the saved correction curves for  plotting
correction_line = np.loadtxt("./correction_linear.txt", skiprows=1)
correction_quad = np.loadtxt("./correction_quad.txt", skiprows=1)
correction_cubic = np.loadtxt("./correction_cubic.txt", skiprows=1)

#********************************************************************

# Plotting the data

txt = ("*Generated from 'polDerivative.py' on the\
      \nGitHub Repository: RamanSpec_IntensityCalibration ")

# FIGURE 0 INITIALIZED

plt.figure(0)
ax0 = plt.axes()
plt.title('Fitting result', fontsize=22)

plt.plot(xaxis,  correction_line, 'r', linewidth=3, label='line_fit')
plt.plot(xaxis,  correction_quad, 'g', linewidth=3, label='quad_fit')
plt.plot(xaxis,  correction_cubic, 'b', linewidth=3, label='cubic_fit')

plt.xlabel('Ramanshift / $cm^{-1}$', fontsize=20)
plt.ylabel('Relative sensitivity', fontsize=20)
plt.grid(True)
ax0.tick_params(axis='both', labelsize =20)
plt.xlim((1755, -1120)) #  change this if the xlimit is not correct
ax0.set_ylim([0, 1.25]) # change this if the ylimit is not enough

ax0.minorticks_on()
ax0.tick_params(which='minor', right='on')
ax0.tick_params(axis='y', labelleft='on', labelright='on')
plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey",\
         transform=plt.gcf().transFigure)
plt.legend(loc='lower left', fontsize=16)

#  For saving the plot
#plt.savefig('fit_output.png', dpi=300)
#********************************************************************
