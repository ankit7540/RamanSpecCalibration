# Purpose : Load the matrix containing polarizability and wavefunctions, interpolate the
# polarizability if needed, and compute the respective matrix elements.

# Load necessary modules
import sys
import numpy as np
import time
from scipy import interpolate
from scipy import integrate
#-------------------------------------------------------------------
#********************************************************************
# python function to compute the first derivative at the first point
# and the last point of an array. For this computation, the first 4
# points are used for the derivative for the first data point. Similarly
# last 4 (x,y) points are used for the derivative at the last point.

def fd_ends(x,y):
    ''' Parameter:
        x       =       xaxis of the array
        y       =       y data of the array
                (x and y arrays must have atleast 4 elements)
        Returns = first derivative at the first point and the last point
    '''
    if (len(x) < 4 or len(y) < 4):
        print("Error : x and y arrays must have 4 elements")

    subx=np.zeros(4)
    suby=np.zeros(4)

    for i in range(0,4):
        subx[i] = x[i]
        suby[i] = y[i]
#        print('posn',i)


    fd01= ((subx[1]*subx[3]+subx[2]*subx[3]+subx[1]*subx[2]-       \
         2*subx[0]*subx[1]-2*subx[0]*subx[3]-2*subx[0]*subx[2]    \
         +3*subx[0]**2)/(-subx[3]+subx[0])/(-subx[1]+subx[0])      \
         /(-subx[2]+subx[0])*suby[0]-(-subx[2]+subx[0])*(-subx[3] \
         +subx[0])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]  \
         -subx[2])*suby[1]+(-subx[1]+subx[0])*(-subx[3]+subx[0])  \
         /(subx[2]-subx[3])/(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]   \
         -(-subx[1]+subx[0])*(-subx[2]+subx[0])/(subx[2]-subx[3])/(subx[1] \
         -subx[3])/(-subx[3]+subx[0])*suby[3] )

    fd1 = ((subx[1]*subx[3]+subx[2]*subx[3]+subx[1]*subx[2]-2*subx[0]*subx[1] \
          -2*subx[0]*subx[3]-2*subx[0]*subx[2]+3*subx[0]**2)/(-subx[3]+subx[0]) \
          /(-subx[1]+subx[0])/(-subx[2]+subx[0])*suby[0]-(-subx[2]+subx[0])    \
          *(-subx[3]+subx[0])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]-subx[2])\
          *suby[1]+(-subx[1]+subx[0])*(-subx[3]+subx[0])/(subx[2]-subx[3])          \
          /(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]-(-subx[1]+subx[0])          \
          *(-subx[2]+subx[0])/(subx[2]-subx[3])/(subx[1]-subx[3])/(-subx[3]+subx[0])\
          *suby[3] )

    for i in range(0,4):
        subx[i] = x[int(i-4)]
        suby[i] = y[int(i-4)]
#        print (i, int(i-4))


    fdn = ( (subx[1]-subx[3])*(subx[2]-subx[3])/(-subx[3]+subx[0])/(-subx[1] \
           +subx[0])/(-subx[2]+subx[0])*suby[0]-(-subx[3]+subx[0])*(subx[2]   \
           -subx[3])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]-subx[2])   \
           *suby[1]+(-subx[3]+subx[0])*(subx[1]-subx[3])/(subx[2]-subx[3])    \
           /(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]-(-2*subx[0]*subx[3]  \
           -2*subx[1]*subx[3]-2*subx[2]*subx[3]+subx[0]*subx[1]+subx[0]       \
           *subx[2]+subx[1]*subx[2]+3*subx[3]**2)/(subx[2]-subx[3])/(subx[1]   \
           -subx[3])/(-subx[3]+subx[0])*suby[3]  )

    return(fd1,fdn)

#********************************************************************

# Spline function taken from Numerical Recipes in FORTRAN, page 109 and 110
# Numerical recipes in FORTRAN, Second Edition, Press, Teukolsky, Vetterling, Flannery
# Cambridge University Press, 1992
def spline (x, y, yp1, ypn):

    '''Parameters:
        x       =       1D vector of x-values in increasing order.
        y       =       1D vector of y-values
        n       =       number of elements in xVector (and yVector)
        yp1     =       first derivative of the interpolating function at the first segment
        ypn     =       first derivative of the interpolating function at the last segment

    '''
    nx=len(x)
    ny=len(y)

    if (nx == ny):
        n=nx
    else :
        print("Error : x and y data have different lengths in spline.")
        quit()

    u=np.zeros(n)
    y2=np.zeros(n) # this is the output
    p = 0.0
    sig = 0.0

    if yp1 > 1e30 :     # lower boundar condition 'natural'
        y2[0] = 0.0
        u[0]  = 0.0
    else :              # specified first derivative
        y2[0] = -0.5
        u[0]=(3/(x[1]-x[0])) * ( ((y[1]-y[0])/ (x[1]-x[0])) - yp1)

    for i in range(1,n-1) :     #       Decomposition loop of tridiagonal algorithm. y2 and u are temporary.
        sig = ( (x[i]-x[i-1])/(x[i+1]-x[i-1]) )
        p = (sig * y2[i-1]) + 2.0
        y2[i] = (sig-1.0)/p
        u[i] = ((6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1]) / (x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p)
        # print("first loop:",i)

    if ypn > 1e30 :     # upper boundary condition 'natural'
        qn=0.0
        un=0.0
    else :              # specified first derivative
        qn=0.5
        un = (3.0 /(x[n-1]-x[n-2]))*(ypn-((y[n-1]-y[n-2])/(x[n-1]-x[n-2])) )

    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0)

    for k in range(n-2,-1,-1):           # from second last point to the second point
        y2[k]=y2[k]*y2[k+1]+u[k]        # backsubstitution loop of tridiagonal algorithm
        # print("loop 2 :",k)

    return(y2)


#************************************************************************
# Spline interpolation function taken from Numerical Recipes in FORTRAN, page 109 and 110
# Numerical recipes in FORTRAN, Second Edition, Press, Teukolsky, Vetterling, Flannery
# Cambridge University Press, 1992

def splint(xa, ya, y2a, x):
    ''' Parameters :
        xa      =       original x-axis 1D vector
        ya      =       original y-axis 1D vector
        y2a     =       output of the spline function
        x       =       new x axis, scalar
    '''

    nxa=len(xa)
    nya=len(ya)
    ny2a=len(y2a)

    if (nxa != nya or nxa != ny2a or nya != ny2a):
        print("Error : xa or ya or y2a have incorrect dimension(s).")
        quit()

    n = nxa

    klo = int(0)
    khi = int(n-1)
    k = int(0)
    h = 0.0
    element=0.0

    while ((khi-klo) > 1) :
        k = int((khi+klo)/2)
        element=xa[k]
#        print(element,xa[k],k,x)
        if ( element > x) :
            khi = k
        else:
            klo = k

    h = xa[khi] - xa[klo]

    if h == 0 :
        print("Error : Bad xa input in splint")
        quit()

    a = (xa[khi]-x)/h
    b = (x-xa[klo])/h
    y = a*ya[klo]+b*ya[khi] + (( (a**3)-a)*y2a[klo]+( (b**3)-b)*y2a[khi])*(h**2)/6.0

    return(y) #returns the interpolated value

#************************************************************************
#************************************************************************
#************************************************************************

# Loading polarizability data and checks  ::::


# Load the polarizability data ( alpha_xx and alpha_zz)---
alpha_xx=np.loadtxt("matrix_xxf.txt")
alpha_zz=np.loadtxt("matrix_zzf.txt")
omega=np.loadtxt("freq.txt")
distance=np.loadtxt("distance.txt")

# check size of the arrays -------------------------------
#print("Dimensions of isotropy matrix :",alpha_xx.shape)
#print("Dimensions of anisotropy matrix :",alpha_zz.shape)
if not(alpha_xx.shape ==  alpha_zz.shape  or len(omega)==alpha_xx.shape[0] ):
    print("Dimension check on polarizability data matrices or wavelength file failed.")
    quit()
else :
    print("Polarizability data dimension checked.")
    print("\n")
    omega_nm=(1e7/(omega*219474.6313702000))
    omega_A=(omega_nm*10)
    print("Give  rovibME.compute  command with parameters:")
    print("\trovibME.compute(molecule, bra_v, bra_J, ket_v, ket_J, lambda, unit of lambda, operator)")
    print('\t for example:  rovibME.compute("H2",0,2,0,4,488,"n","mp")  ')
    print("\t\t")

    print('\t\tmolecule = for H2 enter "H2", for D2 enter "D2", for HD enter "HD" ')
    print("\t\tbra_v    = vibrational state, v=[0,4]")
    print("\t\tbra_J    = rotataional state, J=[0,10]")
    print("\t\tket_v    = vibrational state, v=[0,4]")
    print("\t\tket_J    = rotataional state, J=[0,10]")
    print("\t\tlambda   = wavelength in Hartree, nm or Angstrom")
    print('\t\tunit of lambda =  for  Hartree           use "H" or "h"  ')
    print('\t\t\t          for  nanometers        use "n" or "nm" ')
    print('\t\t\t          for  Angstrom          use "a" or "A"  ')
    print("\t\tAvailable wavelength range: {0} - {1} Hartree;\n   \t\t\t\t\t    {2} - {3} nm; \n   \t\t\t\t\t    {4} - {5} Angstrom".format(round(omega[0],4),round(omega[-1],4),round(omega_nm[0],4),round(omega_nm[-1],4),round(omega_A[0],4),round(omega_A[-1],4)   ))
    print('\t\toperator	= alpha(perpendicular) = alpha_xx given by "xx" or "x" ')
    print('\t\t\t          alpha(parallel) = alpha_xx given by "zz" or "z" ')
    print('\t\t\t          isotropy or mean polarizability given by "iso" or "mp" or "mean" ')
    print('\t\t\t          anisotropy or polarizability difference or gamma given by "aniso" or "g"  or "diff" ')
    print('\t\t\t          for all the above use "all"   or  "All"or "ALL" ')

    print("...ready.")
#********************************************************************
# the actual function for computation of the rovibrational matrix element.
# vl, Jl, vr , Jr are numbers
# mol, wavelength unit and operator are string, hence need quotes. 

def compute(mol, vl, Jl, vr, Jr, wavelength, wavelength_unit, operator):
    '''#  parameters:
    # mol  =    molecule (for H2 enter "H2", for D2 enter "D2", for HD enter "HD")
    # vl   =    vibrational state for the bra, vl = [0,4]
    # Jl   =    rotational state for the bra,  Jl = [0,10]
    # vr   =    vibrational state for the ket, vr = [0,4]
    # Jr   =    rotational state for the ket,  Jr = [0,10]
    # wavelength =  wavelength ( can be Hartree, nanometers or Angstrom)
    # wavelength_unit = specify unit using the specifier
                                ( for  Hartree           use "H" or "h"  )
                                ( for  nanometers        use "n" or "nm"  )
                                ( for  Angstrom          use "a" or "A"  )

    # operator   = property namely alpha_xx, alpha_zz, mean polarizability
                                   (isotropy)[\bar{alpha}], anisotropy[\gamma]
                                   Specify operator using the specifier.
                                 ( for  alpha_xx = alpha(?)         use "x"     or  "xx"  )
                                 ( for  alpha_zz          use "z"     or  "zz"  )
                                 ( for  isotropy          use "iso"   or  "mp" or "mean" )
                                 ( for  anisotropy        use "aniso" or  "g"  or "diff" )
                                 ( for  all the above     use "all"   or  "All"or "ALL"  )

    This function runs on both Python 2.7x and 3.x
    '''

    # set a dictionary for output array
    d={'output':[0]}
    #----------------------------------------------------------------
    # interpolation function defined here and is used later.
    def interpolate2D_common(input2D,originalx,finalx):
        inputSize=input2D.shape
        tempx=np.zeros((len(originalx),1))
        tempx[:,0]=originalx
        col=np.zeros((len(originalx),1))
        outputArray=np.zeros((1,inputSize[1]))
        for i in range(0,inputSize[1]):
            col[:,0]=input2D[:,i]
            der=(0.0,0.0)	# derivatives at first and last ends
            der=fd_ends(tempx,col)
            # print(i,der[0],der[1])
            secarray=spline(tempx,col,der[0],der[1])
            interp=splint(tempx,col,secarray,finalx)
            outputArray[:,i]=interp
        d['output']=outputArray.T
    #----------------------------------------------------------------

    # Load the polarizability data ( alpha_xx and alpha_zz)
    alpha_xx=np.loadtxt("./data/matrix_xxf.txt")
    alpha_zz=np.loadtxt("./data/matrix_zzf.txt")
    omega=np.loadtxt("./data/freq.txt")
    dist=np.loadtxt("./data/distance.txt")
    distance= np.asarray(dist)
    omega_nm=(1e7/(omega*219474.6313702000)) # convert the original freq to nm

    # compute the isotropy(mean polarizability) and anisotropy (gamma)
    isotropy=np.absolute(2*(np.array(alpha_xx))+np.array(alpha_zz))/3
    anisotropy=np.absolute(np.array(alpha_zz)-np.array(alpha_xx))

    # step 1: load the required wavefunctions ------------------------
    Wfn1="./wavefunctions/{0}v{1}J{2}_norm.txt".format(mol,vl,Jl)
    Wfn2="./wavefunctions/{0}v{1}J{2}_norm.txt".format(mol,vr,Jr)
    r_wave="./wavefunctions/r_wave.txt"
    #print(Wfn1,Wfn2)
    if vl < 0 or vr < 0 or vl > 4 or vr > 4 :
        print("v value out of range. vl and vr = [0,4].")
        quit()

    if Jl < 0  or Jr < 0 or Jl > 10 or Jr > 10 : 
        print("J value out of range. Jl and Jr =[0,10].")
        quit()

    if not (mol == "H2"  or mol == "HD" or mol == "D2" ): 
        print("Incorrect molecule chosen. For H2 enter H2, for D2 enter D2, for HD enter HD. Use quotes." )
        quit()

    # Proceed to load wavefunctions.
    psi1=np.loadtxt(Wfn1)
    psi2=np.loadtxt(Wfn2)
    rwave=np.loadtxt(r_wave)
    #print(len(psi1),len(psi2),len(rwave))
    #----------------------------------------------------------------

    wv=float(wavelength) # entered wavelength is a float.

    if (wavelength_unit == "h" or wavelength_unit == "H"):
        omegaFinal=(1e7/(wv*219474.6313702000))
    elif (wavelength_unit == "n" or wavelength_unit == "nm"):
        omegaFinal=wv
    elif (wavelength_unit == "a" or wavelength_unit == "A"):
        omegaFinal=(wv/10)
    elif not (wavelength_unit == "h" or wavelength_unit == "H" or wavelength_unit == "n" or wavelength_unit == "nm"
              or wavelength_unit == "a" or wavelength_unit == "A" ):
        print("Default unit of nm will be used.")
        omegaFinal=wv

    #print("Wavelength in nanometer : {0}, Hartree : {1}".format(round(omegaFinal,6) , round((1e7/(omegaFinal*219474.63137020)),6) ) )

    #print(omega_nm[0],omega_nm[-1])
    if omegaFinal < omega_nm[0] or omegaFinal > omega_nm[-1]:
        sys.exit("Requested wavelength is out of range.")
        #print("Requested wavelength is out of range.")
        #quit()
    #--------------------------------------------------------------
    n=0
    if (operator == "x" or operator == "xx" ):
        param=alpha_xx
        name=["alpha_xx"]
        n=1
    elif (operator == "z" or operator == "zz" ):
        param=alpha_zz
        name=["alpha_zz"]
        n=1
    elif (operator == "mean" or operator == "mp" or operator == "iso" ):
        param=isotropy
        name=["isotropy"]
        n=1
    elif (operator == "diff" or operator == "g" or operator == "aniso" ):
        param=anisotropy
        name=["anisotropy"]
        n=1
    elif (operator == "all" or operator == "All" or operator == "ALL") :
        list = [alpha_xx,alpha_zz,isotropy,anisotropy]
        name=["alpha_xx","alpha_zz","isotropy","anisotropy"]
        n=4
    else :
        print("Operator not correctly specified.")
        quit()

    for i in range(n):    # evaluation of  interpolation and integral

        # interpolate to the asked wavelength ----------------
        if not(n == 1):
            param=list[i]
            #print(param.shape)
            interpolate2D_common(param,omega_nm,omegaFinal)
        else :
            interpolate2D_common(param,omega_nm,omegaFinal)

        parameter=np.zeros((len(distance),1))
        temp=d['output'] 
        parameter[:,0]=temp[:,0] 
        #-----------------------------------------------------
        # OPTIONAL  export the interpolated data to txt
        # and plot image 
        #np.savetxt("intY.txt",parameter,fmt='%1.9e')
        #np.savetxt("distance.txt",distance,fmt='%1.9e')

#        x2=[0]
#        y2=[0]
#        x3=[0]
#        y3=[0]
#        plotter2d( distance,parameter, x2, y2,x3,y3,'r-','w-','w-',name[i],1,350)
        #-----------------------------------------------------

        # step 2: generate interpolated parameter for same xaxis as psi
        # and interpolate the parameter to same dimension as psi
        der2=(0.0,0.0)
        der2=fd_ends(distance,parameter)
#        print('derivatives:',der2[0],der2[1])
        secarray2=spline(distance,parameter,der2[0],der2[1])
        parameter_interp=np.zeros(len(rwave))
        for j in range(0,len(rwave)):
            parameter_interp[j]=splint(distance,parameter,secarray2,rwave[j])

        # step 3: compute the pointwise products
        p1=np.multiply(psi1,psi2)
        p2=np.multiply(p1,rwave)
        p3=np.multiply(p2,rwave)
        product=np.multiply(p3,parameter_interp)

        # function defining the integrand which uses the spline coef array to give interpolated values
        def integrand(xpoint):
            result=splint(rwave,product,secarray2,xpoint)
            return result

        # step 4: gen cubic spline coefs for product
        der2=fd_ends(rwave,product)
        secarray2=spline(rwave,product,der2[0],der2[1])

        # step 5: compute the integral using adaptive Quadrature
        result=integrate.quadrature(integrand,0.2,4.48,tol=1.0e-6,vec_func=False,maxiter=1000)
        rounderr=round(result[1],7)
        print("{0} < v={1}J={2} | {3} | v={4}J={5} >  =  {6} a.u. (err : {7}) ". format(mol,vl,Jl,name[i],vr,Jr,abs(round(result[0],7)  ) , rounderr ) )    
        
#********************************************************************

