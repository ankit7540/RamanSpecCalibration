# -*- coding: utf-8 -*-

"""Module describing the weighted non-linear optimization scheme used to
determine the wavelength sensitivity of the spectrometer using a  polynomial
as a model function"""

# Load necessary modules
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

#********************************************************************

def residual(param):
    '''Function which computes the difference between the ratio of the 'curve'
    to the true given as 'ydata'  with that from  the function describing the
    perturbation. Here it is modeled as a quadratic function of
    the type, ( 1+ c1*x + c2*x**2 ) '''

    c1 = param[0]
    c2 = param[1]


    rhs_n = 1.0 + c1/100 * position_pos + c2/100 * (position_pos**2)
    rhs_d = 1.0 + c1/100 * position_neg + c2/100 * (position_neg**2)

    resd = weight * (ratio_et - (rhs_n/rhs_d))**2
    return np.sum(resd)

#********************************************************************

print("************************************************\n")

# Generate artificial data for fitting

xdata = np.linspace(-10.0, 10.0, num=200) #xaxis
ydata = 1 + (-0.0075)*xdata +  (0.005) *(0.25-xdata)**2 #ydata

# Generate noise vector
noise = np.zeros(200)
for i in range(200):
    val = (ydata[i]**2)/65
    noise[i] = np.random.normal(0, val)

# True coefficients for the perturbation function
k1_true = 1.0855
k2_true = -1.6

# Generate a perturbation which affects the ydata
pert = 1.0 + (k1_true/100) *xdata + (k2_true/100) *xdata**2  # quadratic

curve = np.multiply(pert, ydata)    # affected ydata
curve = curve + noise  # with noise

#***************************************************************

# sampling from the curve at defined data points

position_pos = np.array([xdata[105], xdata[115], xdata[125], xdata[145], xdata[150]\
                         , xdata[175], xdata[192]])

position_neg = np.array([xdata[95], xdata[85], xdata[75], xdata[65], xdata[50] ,\
                         xdata[40], xdata[20]])

# ------------------------

curve_pos = np.array([curve[105], curve[115], curve[125], curve[145], curve[150] ,\
                      curve[175], curve[192]])

curve_neg = np.array([curve[95], curve[85], curve[75], curve[65], curve[50] ,\
                      curve[40], curve[20]])
# ------------------------

y_pos = np.array([ydata[105], ydata[115], ydata[125], ydata[145], ydata[150] ,\
                  ydata[175], ydata[192]])

y_neg = np.array([ydata[95], ydata[85], ydata[75], ydata[65], ydata[50] ,\
                  ydata[40], ydata[20]])
# ------------------------

err_pos = np.array([noise[105], noise[115], noise[125], noise[145], noise[150] ,\
                    noise[175], noise[192]])

err_neg = np.array([noise[95], noise[85], noise[75], noise[65], noise[50] ,\
                    noise[40], noise[20]])

# ------------------------

sigma = (y_pos / y_neg)  * ((err_pos/y_pos)**2)+((err_neg/y_neg)**2)
weight = 1/(sigma**2)

#print(weight)
weight = weight/1e5 # scale the weight

ratio_et = np.zeros(position_pos.shape)      # ratio of yval , curve / true
ratio_true = np.zeros(position_pos.shape)    # ratio of yval from ydata at the mentioned data points
ratio_curve = np.zeros(position_pos.shape)   # ratio of yval from curve at the mentioned data points

#Generate the ratio of the y-values

for i in range(len(position_pos)):
    ratio_curve[i] = (curve_pos[i] / curve_neg[i])
    ratio_true[i] = (y_pos [i] / y_neg [i] )
    ratio_et[i] = ratio_curve[i] / ratio_true[i]

print("Number of pairs of data points used :", len(ratio_et))

#***************************************************************

# Intial guess
init_k1 = 1.01
init_k2 = -2.00

param_init = np.array([ init_k1, init_k2 ])
print("Testing the residual function with data")
print("Initial coef :  k1={0}, k2={1} output = {2}".format(init_k1, init_k2,\
      (residual(param_init))))

print("\n***************************************************************")

print("\nOptimization run     \n")
res = opt.minimize(residual, param_init, method='Nelder-Mead', \
                              options={'xatol': 1e-9, 'fatol': 1e-9})

print(res)
k1 = res.x[0]
k2 = res.x[1]
print("\nOptimized result : k1={0}, k2={1} ".format(round(k1, 6), round(k2, 6)))

correction_curve = 1+(k1/100)* xdata +(k2/100)* (xdata**2)  # generate the correction curve

corrected = curve / correction_curve # perform correction to obtain the corrected array

print("***************************************************************")



#********************************************************************
# Plotting the data

txt = ("*Generated from 'polDerivative.py' on the\
      \nGitHub Repository: RamanSpec_IntensityCalibration ")

# FIGURE 0 INITIALIZED
plt.figure(0)
ax0 = plt.axes()
plt.title('Fitting results when using the ratio of \ny-values at sparsely sampled pair of data points')
plt.plot(xdata, ydata, 'r', linewidth=3)
plt.plot(xdata, curve, 'b', linewidth=2.0)
plt.plot(xdata, correction_curve, 'k--', linewidth=1.5)

plt.plot(xdata, noise, 'k', linewidth=1.0)

plt.xlabel('x-axis')
plt.ylabel('y-value')
plt.grid(True)

ax0.set_ylim([-1.5, 3.9]) # change this if the ydata function or perturbation is changed

ax0.minorticks_on()
ax0.tick_params(which='minor', right='on')
ax0.tick_params(axis='y', labelleft='on', labelright='on')
plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey", \
         transform=plt.gcf().transFigure)
plt.legend(('true', 'perturbed', 'correction_curve (from optimization)',\
            'corrected'), loc='upper left')

#********************************************************************