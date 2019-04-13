# Purpose : Load the matrix containing polarizability and wavefunctions, interpolate the
# polarizability if needed, and compute the respective matrix elements.

# Load necessary modules
import sys
import numpy
from scipy import interpolate
from scipy import integrate
#-------------------------------------------------------------------
#********************************************************************


def residual(vars, x, data, eps_data):
    amp = vars[0]
    phaseshift = vars[1]
    freq = vars[2]
    decay = vars[3]

    model = amp * sin(x*freq + phaseshift) * exp(-x*x*decay)

    return (data-model) / eps_data

from scipy.optimize import leastsq

vars = [10.0, 0.2, 3.0, 0.007]


out = leastsq(residual, vars, args=(x, data, eps_data))