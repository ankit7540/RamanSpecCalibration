# -*- coding: utf-8 -*-

"""Module describing the concept of the non-linear optimization scheme used to
    determine the wavelength sensitivity of the spectrometer using a  polynomial
    as a model function"""

# Load necessary modules
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

print("***************************************************************\n")

#-------------------------------------------------------------------
#********************************************************************
# Generate noise free x and y data

xdata = np.linspace(-10.0, 10.0, num=100)       # xaxis
ydata = 1+ (0.04) *xdata                        # true ydata

# Generate a perturbation which affects the ydata
print("True coefs for the perturbation, k1 : {0}, k2 :{1}".format(0.014, -0.005))
pert = 1+ (0.014) *xdata + (-0.005) *xdata**2   # perturbation

curve = np.multiply(pert, ydata)                # affected ydata

#********************************************************************
# Asymmetric data sampling from the perturbed y-data (curve) at the
#   following data points

# true ratio (taken at 2,3,4,6,7,9 on positive side to -1.5,-2.5,-5.5,
#             -7.5, and -9.5 on the negative side respectively)

position_pos = np.array([2, 3, 4, 6, 7, 9])
position_neg = np.array([-1.5, -2.5, -3.5, -5.5, -7.5, -9.5])

#Generate numpy arrays to keep the interpolated data

ratio_true = np.zeros(position_pos.shape)    # ratio of yval from ydata at the mentioned data points
ratio_curve = np.zeros(position_pos.shape)   # ratio of yval from curve at the mentioned data points
ratio_et = np.zeros(position_pos.shape)      # ratio of yval , curve / true

for i in range(len(position_pos)):

    ratio_true[i] = (np.interp(position_pos[i], xdata, ydata)) /  \
    (np.interp(position_neg[i], xdata, ydata))

    ratio_curve[i] = (np.interp(position_pos[i], xdata, curve)) / \
    (np.interp(position_neg[i], xdata, curve))
    ratio_et[i] = ratio_curve[i] / ratio_true[i]


print("Number of pairs of data points used :", len(ratio_et))


#********************************************************************

def residual(param):
    '''Function which computes the difference between the ratio of the 'curve'
    to the true given as 'ydata'  with that from  the function describing the
    perturbation. Here it is modeled as a quadratic function of
    the type, ( 1+ c1*x + c2*x**2 ) '''

    c1 = param[0]
    c2 = param[1]

    rhs_n = 1 + c1 * position_pos + c2 * (position_pos**2)
    rhs_d = 1 + c1 * position_neg + c2 * (position_neg**2)

    #print(rhs_n)
    #print(rhs_d)

    resd = (ratio_et - (rhs_n/rhs_d))**2
    return np.sum(resd)

#********************************************************************

print("***************************************************************")

param_init = np.array([-0.0125, -0.001])
print("Testing the residual function with data")
print("Initial coef :  k1={0}, k2={1} output = {2}".format(-0.014, -0.001,\
      (residual(param_init))))

print("***************************************************************")

print("\nOptimization run     \n")
res = opt.minimize(residual, param_init, method='Nelder-Mead', \
                              options={'xatol': 1e-9, 'fatol': 1e-9})

print(res)
k1 = res.x[0]
k2 = res.x[1]
print("\nOptimized result : k1={0}, k2={1} ".format(round(k1, 6), round(k2, 6)))


correction_curve = (1+k1* xdata +k2* (xdata**2))  # generate the correction curve

corrected = curve / correction_curve # perform correction to obtain the corrected array

print("***************************************************************")

#********************************************************************
# Plotting the data

txt = ("*Generated from 'polDerivative.py' on the\
      \nGitHub Repository: RamanSpec_IntensityCalibration ")

# FIGURE 0 INITIALIZED
plt.figure(0)
ax0 = plt.axes()
plt.title('Fitting results when using the ratio of \ny-values at sparsely \
          sampled pair of data points')
plt.plot(xdata, ydata, 'r', linewidth=3)
plt.plot(xdata, curve, 'b', linewidth=2.0)
plt.plot(xdata, correction_curve, 'k--', linewidth=1.5)
plt.plot(xdata, corrected, 'Y', linewidth=1.5)
ax0.set_ylim([0.0, 1.725])

plt.xlabel('x-axis')
plt.ylabel('y-value')
plt.grid(True)

ax0.minorticks_on()
ax0.tick_params(which='minor', right='on')
ax0.tick_params(axis='y', labelleft='on', labelright='on')
plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey", \
         transform=plt.gcf().transFigure)
plt.legend(('true', 'perturbed', 'correction_curve (from optimization)',\
            'corrected'), loc='upper left')

#********************************************************************

xaxd = np.array([1, 2, 3, 4, 5, 6])

interp_fit = np.zeros(position_pos.shape)

for i in range(len(position_pos)):
    interp_fit[i] = (np.interp(position_pos[i], xdata, corrected)) /\
    (np.interp(position_neg[i], xdata, corrected))
#********************************************************************

# FIGURE 1 INITIALIZED

plt.figure(1)
ax0 = plt.axes()
plt.title('Ratio of y-values sampled  at the xaxis\n position in positive \
          and negative region')
plt.legend(loc='upper left')

plt.plot(xaxd, ratio_true, 'ro-', linewidth=4.5)
plt.plot(xaxd, ratio_curve, 'bo', linewidth=2.0)
plt.plot(xaxd, interp_fit, 'go-', linewidth=2.25)
plt.legend(loc='upper left')

plt.xlabel('index for the sampled data points')
plt.ylabel('value of the ratio')
plt.grid(True)
ax0.minorticks_on()
ax0.tick_params(axis='x', which='minor', bottom=False)

ax0.tick_params(which='minor', right='on')
ax0.tick_params(axis='y', labelleft='on', labelright='on')
plt.legend(('ratio_true',\
            'ratio_curve (perutbed)',\
            'ratio_corrected'), loc='upper left')

plt.text(0.05, 0.00001, txt, fontsize=5, color="dimgrey",\
         transform=plt.gcf().transFigure)


#********************************************************************
