# Purpose : Load the matrix containing polarizability and wavefunctions, interpolate the
# polarizability if needed, and compute the respective matrix elements.

# Load necessary modules
import sys
import numpy
from scipy import interpolate
from scipy import integrate
from scipy.optimize import leastsq

#-------------------------------------------------------------------
#********************************************************************


def residual(vars, x, data, eps_data):
    amp = vars[0]
    phaseshift = vars[1]
    freq = vars[2]
    decay = vars[3]

    model = amp * sin(x*freq + phaseshift) * exp(-x*x*decay)

    return (data-model) / eps_data

def residual_line(vars, x, data, error):
    intercept = vars[0]
    slope = vars[1]
    
    model = intercept + slope*x
    
    return (data-model) / error


def residual_quadratic(vars):
    intercept = vars[0]
    slope = vars[1]
    
    model = amp * sin(x*freq + phaseshift) * exp(-x*x*decay)
    return (data-model) / eps_data


vars = [3, 4]


out = leastsq(residual_line, vars, args=(x, data, error))