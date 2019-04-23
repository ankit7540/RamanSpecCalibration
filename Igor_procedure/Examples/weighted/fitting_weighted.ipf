
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//******************************************************************************************************
//******************************************************************************************************

// Procedure   which generates sample x-y data, a quadratic functions as perturbation and the 
//	the perturbed function decribed as 'curve'

// Using ratio of the yvalues from  pair of x-positions the fitting is done to obtain the polynomial
//	coefficients.

// Regarding the choice of the sampling points, 
//	 * In this example, the x-points are taken in positive and negative region but these 
//	can be randomly chosen.
//	 * If more data points are sampled then the accuracy increases.
//******************************************************************************************************
//******************************************************************************************************

function generate_noisy_data( )

	make /d /o /n=100 xdata = -10+p/5							// xaxis
	make /d /o /n=100 ydata = 0+ 0.04 *xdata					// ydata, true
	make /d /o /n = 100  pert = 1+ 0.014* xdata  - 0.005*(xdata^2) 	// perturbation, true coefs : k1=0.014, k2=-0.005
	make /d /o /n = 100 noise =0
	variable i
	for ( i=0 ; i<100 ; i =i+1)
		noise[i] = enoise( (ydata[i]^2)/50 )
	endfor	
		
	
	make /d /o /n = 100 curve = (pert * ydata) + noise						// perturbed data
	
	make /d /o /n=6 posn_pos, posn_neg							
	posn_pos={ xdata[56] , xdata[64], xdata[75], xdata[80], xdata[85], xdata[90], xdata[95]}	// x-axis positions, on the positive side
	posn_neg={ xdata[47] , xdata[41], xdata[35], xdata[30], xdata[25], xdata[15], xdata[2]} // x-axis positions, on the negative side
																// these can be randomly chosen

	variable nRows = dimsize(posn_neg, 0)
	
	// waves to keep the interpolation and the ratio of yvalues	
	make  /o/d /n=(nRows) interp_pos		// interpolated y val, positive x
	make  /o/d /n=(nRows) interp_neg		// interpolated y val, negative x

	make  /o/d /n=(nRows) cY_pos		// y values in curve
	make  /o/d /n=(nRows) cY_neg		// y values in curve

	make  /o/d /n=(nRows) ratio_true			// ratio of y values, in the true
	make  /o/d /n=(nRows) ratio_curve		// ratio of  y values, in the perturbed
	make  /o/d /n=(nRows) ratio_et			//  ratio of y values, perturbed / true
	make  /o/d /n=(nRows) sigma			//  sigma_ratio
	make  /o/d /n=(nRows) weight
	
	interp_pos={ ydata[56] , ydata[64], ydata[75], ydata[80], ydata[85], ydata[90], ydata[95]}  // yValues, true
	interp_neg={ ydata[47] , ydata[41], ydata[35], ydata[30], ydata[25], ydata[15], ydata[2]}    // yValues, true 	
	
	ratio_true = 	interp_pos /  interp_neg	
			
	interp_pos={ curve[56] , curve[64], curve[75], curve[80], curve[85], curve[90], curve[95]}  // yValues, perturbed
	interp_neg={ curve[47] , curve[41], curve[35], curve[30], curve[25], curve[15], curve[2]}    //	 yValues, perturbed
	
	cY_pos = interp_pos		// yValues, perturbed, for later use
	cY_neg = interp_neg		// yValues, perturbed, for later use
	
	ratio_curve  = 	interp_pos /  interp_neg		
	ratio_et =  ratio_curve / ratio_true
	
	interp_pos={ noise[56] , noise[64], noise[75], noise[80], noise[85], noise[90], noise[95] } // has error
	interp_neg={ noise[47] , noise[41], noise[35], noise[30], noise[25], noise[15], noise[2] }   // has error	
	
	sigma = (ratio_curve) * ( ( interp_pos / cY_pos )^2 + ( interp_neg / cY_neg )^2 )
	weight = 1/(sigma^2)
	weight = weight / 1e9

		
	printf "True y values = ydata\r"
	printf "Perturbed wave = curved\r"
	printf "Perturbation =  quadratic polynomial : 1 + (0.014) x + (-0.005)  x^2\r" 
	
end	

//******************************************************************************************************
//******************************************************************************************************

//	This is the function describing the residual when two parameters (for polynomial) are used as 
//		variables

function  residual_weighted (w,  k1, k2)
	wave w			// optional wave
	variable  k1		// fit parameter
	variable k2		// fit  parameter
	
	make  /o /FREE /d /n=(7) rhs_n	// numerator, RHS
	make  /o /FREE /d /n=(7) rhs_d	//  denominator, RHS
	make  /o /FREE /d /n=(7) resd	// residual, wave
	
	wave posn_pos
	wave posn_neg
	wave ratio_et
	wave weight
	
	rhs_n = 1+k1/100*posn_pos + k2/100*(posn_pos^2)	// scaling k1 and  k2 reduces numerical error
	rhs_d = 1+k1/100*posn_neg + k2/100*(posn_neg^2)	// scaling k1 and  k2 reduces numerical error
	
	resd = weight * ( ratio_et - (rhs_n / rhs_d))^2
	
	return wavesum(resd) //  return the  sum of the  residual wave
end	

//******************************************************************************************************
//******************************************************************************************************

// WaveSum(w)
// Returns the sum of the entire wave, just like Igor’s sum function.

Function WaveSum(w)
	Wave w
	Variable i, n=numpnts(w), total=0
	
	for(i=0;i<n;i+=1)
		total += w[i]
	endfor

return total
End

//******************************************************************************************************
//******************************************************************************************************

// This  function runs the non-linear optimization process using the 'Optimize' function in  IgorPro

function run_Opt  ( temp_wave  , initial_k1,  initial_k2 ) 
	wave temp_wave	// optional wave passed to the residual function
	variable  initial_k1	// guess for k1
	variable  initial_k2	// guess for k2

	
	Optimize /A=0  /I=3000 /M={0 , 0} /T={1e-9, 1e-11}    /X={ initial_k1 , initial_k2 }  residual_weighted , temp_wave  
	//			     iterations  method  tolerance                  initial coefs                    function    optional
	
	wave  result = W_Extremum
	wave xdata
	wave curve

	// k1 was scaled by  10, k2 was scaled by 100
	printf "Optimized parameters : %5.6f, %5.6f\r", result[0] / 100, result[1] / 100
	
	make /o /d /n=100 corn_curve = 1+ ( result[0]/100) * xdata + (result[1]/100)*xdata^2 
	
	make  /o  /d /n=100 corrected  = curve / corn_curve		//  generate the  corrected wave, compare to the true ydata
end

//******************************************************************************************************
//******************************************************************************************************

