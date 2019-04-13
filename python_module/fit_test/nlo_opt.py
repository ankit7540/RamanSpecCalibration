# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:27:16 2019

@author: Ankit Raj
"""

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

def load2D_xy(filename):
    col0 = numpy.loadtxt(filename)[:, 0]
    col1 = numpy.loadtxt(filename)[:, 1]
    return col0, col1

xd1, yd1 = load2D_xy("p1.txt")
#print(xd1, yd1)

#-------------------------------------------------------------------
#********************************************************************

xd2, yd2 = load2D_xy("p2.txt")
#print(xd2, yd2)

#-------------------------------------------------------------------
#********************************************************************

def residual(par):
    k1 = par[0]
    k2 = par[1]
    
    def poly2(k1,k2, x):
        return 1+k1*x+k2*x**2        
    
    
    
    #for i in range(len(xd1)):
    ##    print(xd1[i])
    #return(poly2(k1,k2,2))

va = numpy.array([1.0, 1.0])    
print(residual(va))    