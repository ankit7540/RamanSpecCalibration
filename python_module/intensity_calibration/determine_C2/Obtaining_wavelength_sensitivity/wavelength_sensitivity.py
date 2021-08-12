# -*- coding: utf-8 -*-

"""Module describing the weighted non-linear optimization scheme used to
determine the wavelength sensitivity of the spectrometer using a  polynomial
as a model function. This file is identical in terms of code to the one present
 in the one above directory, except for detailed description."""

# Load necessary modules
import numpy as np
import scipy.optimize as opt
import logging
import matplotlib.pyplot as plt
from datetime import datetime

# Set logging ------------------------------------------
# log file will have name as 'logfile_Y_m_d__H_M'
timestamp=datetime.now().strftime('%Y-%m-%d__%H_%M')
log_filename='logfile'+timestamp;
fileh = logging.FileHandler(log_filename, 'w+')
formatter = logging.Formatter('%(message)s')
fileh.setFormatter(formatter)

log = logging.getLogger()  # root logger
for hdlr in log.handlers[:]:  # remove all old handlers
    log.removeHandler(hdlr)
log.addHandler(fileh)      # set the new handler
# ------------------------------------------------------

# Logging starts here
log.debug( datetime.now().strftime('%Y-%m-%d %H:%M:%S') )
logger= logging.getLogger( __file__ )
log.info(logger)
logging.getLogger().setLevel(logging.INFO)
log.warning('\n',)
log.error("------------ Run log ------------\n")

#********************************************************************
# Load data files
# Modify this as user case
# by default no header should be present
# if header present, add 'skiprows=n' as on line 494

# For structure of data in the input file refer to the readme in the
#  directory containing this file
dataH2 = np.loadtxt("./dataH2.txt")
dataHD = np.loadtxt("./dataHD.txt")
dataD2 = np.loadtxt("./dataD2.txt")
dataO2 = np.loadtxt("./DataO2.txt")
dataO2_p = np.loadtxt("./DataO2_pR.txt")
xaxis = np.loadtxt("./Wavenumber_axis.txt")

# Go to around line 482 for examples of
#  using the above data for fitting

#********************************************************************

# Constants ------------------------------
# these are used for scaling the polynomial coefs
scale1 = 1e4
scale2 = 1e7
scale3 = 1e9
scale4 = 1e11
# ----------------------------------------

# these are used for scaling the weights for O2
# edit as needed
# best to keep weights of all gases within the same magnitude
scale_O2_S1O1 = 2.0
scale_O2_pureRotn= 0.1
# ----------------------------------------

# common log of scaling factors
log.info("Scaling factor for coefs: \n\t for c1 = %3.3e\n\t for c2 = %3.3e\n\t for c3 = %3.3e", scale1, scale2, scale3)

log.info("\nScaling factor for weights: \n\tfor O2, s1-o1 branch = %3.3f\n\tfor O2, s1-o1 branch = %3.3f", scale_O2_S1O1, scale_O2_pureRotn )
log.error("------------ Run log ------------\n")

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

    # residual = weight * (difference of LHS and RHS, squared)
    resd_H2 = dataH2[:, 5] * ((ratio_H2 - RHS_H2)**2)

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

    # residual = weight * scaling factor * (difference of LHS and RHS, squared)
    resd_O2 = ( dataO2[:, 5] * scale_O2_S1O1 ) * ((ratio_O2 - RHS_O2)**2)
	# ------

    # - O2 pure rotation -
    ratio_O2p = dataO2_p[:, 1]/dataO2_p[:, 2]
    RHS_O2p = (1.0 + c1/scale1 * dataO2_p[:, 3] )/ (1.0 +\
             c1/scale1 * dataO2_p[:, 4] )

    resd_O2p = (dataO2_p[:, 5] * scale_O2_pureRotn ) * ((ratio_O2p - RHS_O2p)**2)
	# ------


    return np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) + np.sum(resd_O2)  + np.sum(resd_O2p)

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

    resd_O2 = (dataO2[:, 5] * scale_O2_S1O1 ) * ((ratio_O2 - RHS_O2)**2)
	# ------

    # - O2 pure rotation -
    ratio_O2p = dataO2_p[:, 1]/dataO2_p[:, 2]
    RHS_O2p = (1.0 + c1/scale1 * dataO2_p[:, 3] + c2/scale2 * (dataO2_p[:, 3]**2))/ (1.0 +\
             c1/scale1 * dataO2_p[:, 4] + c2/scale2 * (dataO2_p[:, 4]**2))

    resd_O2p = (dataO2_p[:, 5]* scale_O2_pureRotn ) * ((ratio_O2p - RHS_O2p)**2)
	# ------

    return np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) + np.sum(resd_O2)  + np.sum(resd_O2p)

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

    resd_O2 =( dataO2[:, 5]  * scale_O2_S1O1 ) * ((ratio_O2 - RHS_O2)**2)
	# ------

    # - O2 pure rotation -
    ratio_O2p = dataO2_p[:, 1]/dataO2_p[:, 2]
    RHS_O2p = (1.0 + c1/scale1 * dataO2_p[:, 3] + c2/scale2 * (dataO2_p[:, 3]**2)+\
               c3/scale3 * (dataO2_p[:, 3]**3))/ ( 1.0 + c1/scale1 * dataO2_p[:, 4] +\
                           c2/scale2 * (dataO2_p[:, 4]**2)+ c3/scale3 * ( dataO2_p[:, 4]**3))

    resd_O2p = (dataO2_p[:, 5] * scale_O2_pureRotn  ) * ((ratio_O2p - RHS_O2p)**2)
	# ------

    return np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) + np.sum(resd_O2)  + np.sum(resd_O2p)


#***************************************************************
#********************************************************************

def residual_quartic(param):
    '''Function which computes the residual (as sum of squares) comparing the
    ratio of expt to theoretical intensity ratio to the sensitivity  profile
    modelled as  a quadratic polynomial, ( 1 + c1*x + c2*x**2 + c3*x**3 + c4*x**4 ) '''

    c1 = param[0]
    c2 = param[1]
    c3 = param[2]
    c4 = param[3]

	# - H2 -
    ratio_H2 = dataH2[:, 1]/dataH2[:, 2]
    RHS_H2 = (1.0 + c1/scale1 * dataH2[:, 3] + c2/scale2 *(dataH2[:, 3]**2) +\
              c3/scale3 * (dataH2[:, 3]**3)+c4/scale4 * (dataH2[:, 3]**4))/ ( 1.0 + c1/scale1 * dataH2[:, 4] +\
                          c2/scale2 * (dataH2[:, 4]**2) + c3/scale3 * ( dataH2[:, 4]**3)+\
                          c4/scale4 * ( dataH2[:, 4]**4)    )

    resd_H2 = dataH2[:, 5] * ((ratio_H2 - RHS_H2)**2)
	# ------

	# - HD -
    ratio_HD = dataHD[:, 1]/dataHD[:, 2]
    RHS_HD = (1.0 + c1/scale1 * dataHD[:, 3] + c2/scale2 *(dataHD[:, 3]**2)+\
              c3/scale3 *(dataHD[:, 3]**3)+c4/scale4 *(dataHD[:, 3]**4))/ ( 1.0 + c1/scale1 * dataHD[:, 4] +\
                         c2/scale2 * (dataHD[:, 4]**2)+ c3/scale3 * ( dataHD[:, 4]**3) + \
                             + c4/scale4 * ( dataHD[:, 4]**4) )

    resd_HD = dataHD[:, 5] * ((ratio_HD - RHS_HD)**2)
	# ------

	# - D2 -
    ratio_D2 = dataD2[:, 1]/dataD2[:, 2]
    RHS_D2 = (1.0 + c1/scale1 * dataD2[:, 3] + c2/scale2 * (dataD2[:, 3]**2)+\
              c3/scale3 * (dataD2[:, 3]**3) + c4/scale4 * (dataD2[:, 3]**4) )/ ( 1.0 + c1/scale1 * dataD2[:, 4] +\
                          c2/scale2 * (dataD2[:, 4]**2)+ c3/scale3 * ( dataD2[:, 4]**3) +\
                              + c4/scale4 * ( dataD2[:, 4]**4) )

    resd_D2 = dataD2[:, 5] * ((ratio_D2 - RHS_D2)**2)
	# ------

	# - O2 high frequency -
    ratio_O2 = dataO2[:, 1]/dataO2[:, 2]
    RHS_O2 = (1.0 + c1/scale1 * dataO2[:, 3] + c2/scale2 * (dataO2[:, 3]**2)+\
              c3/scale3 * (dataO2[:, 3]**3) + c4/scale4 * (dataO2[:, 3]**4))/ ( 1.0 + c1/scale1 * dataO2[:, 4] +\
                          c2/scale2 *(dataO2[:, 4]**2)+ c3/scale3 *( dataO2[:, 4]**3) + \
                          c4/scale4 *( dataO2[:, 4]**4) )

    resd_O2 =( dataO2[:, 5]  * scale_O2_S1O1 ) * ((ratio_O2 - RHS_O2)**2)
	# ------

    # - O2 pure rotation -
    ratio_O2p = dataO2_p[:, 1]/dataO2_p[:, 2]
    RHS_O2p = (1.0 + c1/scale1 * dataO2_p[:, 3] + c2/scale2 * (dataO2_p[:, 3]**2)+\
               c3/scale3 * (dataO2_p[:, 3]**3) + c4/scale4 * (dataO2_p[:, 3]**4) )/ ( 1.0 + c1/scale1 * dataO2_p[:, 4] +\
                           c2/scale2 * (dataO2_p[:, 4]**2)+ c3/scale3 * ( dataO2_p[:, 4]**3) + \
                               c4/scale4 * ( dataO2_p[:, 4]**4) )

    resd_O2p = (dataO2_p[:, 5] * scale_O2_pureRotn  ) * ((ratio_O2p - RHS_O2p)**2)
	# ------

    return np.sum(resd_H2) + np.sum(resd_HD) + np.sum(resd_D2) + np.sum(resd_O2)  + np.sum(resd_O2p)

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


    print("\nOptimization run : Linear    \n")
    res = opt.minimize(residual_linear, param_init, method='Nelder-Mead', \
                              options={'xatol': 1e-9, 'fatol': 1e-9})

    print(res)
    optk1 = res.x[0]
    print("\nOptimized result : k1={0} \n".format(round(optk1, 6) ) )

    correction_curve_line= 1+(optk1/scale1)*xaxis     # generate the correction curve

    np.savetxt("correction_linear.txt", correction_curve_line, fmt='%2.8f',\
               header='corrn_curve_linear', comments='')

    print("**********************************************************")

    # save log -----------
    log.info('\n *******  Optimization run : Linear  *******')
    log.info('\n\t Initial : c1 = %4.8f\n', init_k1 )
    log.info('\n\t %s\n', res )
    log.info('\n Optimized result : c1 = %4.8f \n', optk1 )
    log.info(' *******************************************')
    # --------------------


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

    print("\nOptimization run : Quadratic    \n")
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

    # save log -----------
    log.info('\n *******  Optimization run : Quadratic  *******')
    log.info('\n\t Initial : c1 = %4.8f, c2 = %4.8f\n', init_k1, init_k2 )
    log.info('\n\t %s\n', res )
    log.info('\n Optimized result : c1 = %4.8f, c2 = %4.8f\n', optk1, optk2 )
    log.info(' *******************************************')
    # --------------------

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

    print("\nOptimization run : Cubic   \n")
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

    # save log -----------
    log.info('\n *******  Optimization run : Cubic  *******')
    log.info('\n\t Initial : c1 = %4.8f, c2 = %4.8f, c3 = %4.8f\n', init_k1, init_k2, init_k3 )
    log.info('\n\t %s\n', res )
    log.info('\n Optimized result : c1 = %4.8f, c2 = %4.8f, c3 = %4.8f\n', optk1, optk2, optk3 )
    log.info(' *******************************************')
    # --------------------

#***************************************************************
#***************************************************************

def run_fit_quartic(init_k1, init_k2, init_k3, init_k4):
    '''Function performing the actual fit using the residual_quartic function
    defined earlier '''

    #  init_k1  =   initial  guess  for k1
    #  init_k2  =   initial  guess  for k2
    #  init_k3  =   initial  guess  for k3
    #  init_k4  =   initial  guess  for k4

    param_init = np.array([init_k1, init_k2, init_k3, init_k4 ])
    print("**********************************************************")
    print("Initial coef :  k1={0}, k2={1}, k3={2}, output = {3}, , output = {4}".format(init_k1, init_k2,\
          init_k3, init_k4, (residual_cubic(param_init))))

    print("\nOptimization run : Quartic    \n")
    res = opt.minimize(residual_cubic, param_init, method='Nelder-Mead', \
                              options={'xatol': 1e-9, 'fatol': 1e-9})

    print(res)
    optk1 = res.x[0]
    optk2 = res.x[1]
    optk3 = res.x[2]
    optk4 = res.x[3]
    print("\nOptimized result : k1={0}, k2={1}, k3={2}, k4={3}  ".format(round(optk1, 6),\
          round(optk2, 6), round(optk3, 6), round(optk4, 6)))

    correction_curve_cub = 1+(optk1/scale1)*xaxis +(optk2/scale2)*(xaxis**2)+\
    (optk3/scale3)* (xaxis**3) + (optk4/scale4)* (xaxis**4)  # generate the correction curve

    np.savetxt("correction_Quartic.txt", correction_curve_cub, fmt='%2.8f',\
               header='corrn_curve_Quartic', comments='')

    print("**********************************************************")

    # save log -----------
    log.info('\n *******  Optimization run : Quartic  *******')
    log.info('\n\t Initial : c1 = %4.8f, c2 = %4.8f, c3 = %4.8f, c4 = %4.8f\n', init_k1, init_k2, init_k3, init_k4 )
    log.info('\n\t %s\n', res )
    log.info('\n Optimized result : c1 = %4.8f, c2 = %4.8f, c3 = %4.8f, c4 = %4.8f\n', optk1, optk2, optk3, optk4 )
    log.info(' *******************************************')
    # --------------------

#***************************************************************

# Run the fit with initial guesses
# change these as per your case
run_fit_linear(0.796 )
run_fit_linear(-0.196 )

run_fit_quadratic(0.886, +1.05)
run_fit_quadratic(0.886, -1.05)

run_fit_cubic(0.76, -2.05, 0.050)
run_fit_cubic(0.76, 1.205, 0.010)

run_fit_quartic(-0.176, -0.05, 0.075, 0.0001 )
run_fit_quartic(-0.076, +0.05, -0.075, 0.0001 )

# Final run with better guesses from previous runs
# edit the following guess values accordingly
run_fit_linear( -1.63 )
run_fit_quadratic(-1.06, -0.95)
run_fit_cubic(-1.04, -1.21, 0.0130)
run_fit_quartic(-0.006, +0.05, -0.075, 0.0001 )

# Log file generated contains summary of all fits
# the log file will have date and time stamp

# Load the saved correction curves for  plotting
# outputs from last run will be loaded
correction_line = np.loadtxt("./correction_linear.txt", skiprows=1)
correction_quad = np.loadtxt("./correction_quad.txt", skiprows=1)
correction_cubic = np.loadtxt("./correction_cubic.txt", skiprows=1)
correction_quartic = np.loadtxt("./correction_quartic.txt", skiprows=1)

#********************************************************************

# Plotting the data

txt = ("*Generated from 'wavelength_sensitivity.py' on the\
      \nGitHub Repository: RamanSpecCalibration ")

# FIGURE 0 INITIALIZED

plt.figure(0)
ax0 = plt.axes()
plt.title('Fitting result', fontsize=22)

plt.plot(xaxis,  correction_line, 'r', linewidth=3, label='line_fit')
plt.plot(xaxis,  correction_quad, 'g', linewidth=4.2, label='quad_fit')
plt.plot(xaxis,  correction_cubic, 'b--', linewidth=2.65, label='cubic_fit')
plt.plot(xaxis,  correction_quartic, 'k--', linewidth=2.65, label='quartic_fit')

plt.xlabel('Wavenumber / $cm^{-1}$', fontsize=20)
plt.ylabel('Relative sensitivity', fontsize=20)
plt.grid(True)

# change following as needed
ax0.tick_params(axis='both', labelsize =20)
plt.xlim((1755, -1120)) #  change this if the xlimit is not correct
ax0.set_ylim([0, 1.30]) # change this if the ylimit is not enough

ax0.minorticks_on()
ax0.tick_params(which='minor', right='on')
ax0.tick_params(axis='y', labelleft='on', labelright='on')
plt.text(0.05, 0.0095, txt, fontsize=6, color="dimgrey",\
         transform=plt.gcf().transFigure)
plt.legend(loc='lower left', fontsize=16)


#  For saving the plot
#plt.savefig('fit_output.png', dpi=120)
#********************************************************************

# Log file generated contains summary of all fits
# the log file will have date and time stamp
