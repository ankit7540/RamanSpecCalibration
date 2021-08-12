
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

function generate_data( )

	make /d /o /n=100 xdata = -10+p/5							// xaxis
	make /d /o /n=100 ydata = 1+ 0.04 *xdata					// ydata, true
	make /d /o /n = 100  pert = 1+ 0.014* xdata  - 0.005*(xdata^2) 	// perturbation
	make /d /o /n = 100 curve = pert * ydata						// perturbed data
	
	make /d /o /n=6 posn_pos, posn_neg							
	posn_pos={2, 3, 4,  5, 6, 7, 9.25}								// x-axis positions, on the positive side
	posn_neg={-1.5, -2.5,  -3.5, -4.5, -5.5, -6.5, -8.85}				// x-axis positions, on the negative side
																// these can be randomly chosen

	variable nRows = dimsize(posn_neg, 0)
	
	// waves to keep the interpolation and the ratio of yvalues	
	make  /o/d /n=(nRows) interp_pos		// interpolated y val, positive x
	make  /o/d /n=(nRows) interp_neg		// interpolated y val, negative x
	make  /o/d /n=(nRows) ratio_true			// ratio of y values, in the true
	make  /o/d /n=(nRows) ratio_curve		// ratio of  y values, in the perturbed
	make  /o/d /n=(nRows) ratio_et			//  ratio of y values, perturbed / true
	
	Interpolate2 /i=3 /X = posn_pos /Y=interp_pos xdata, ydata
	Interpolate2 /i=3 /X = posn_neg /Y=interp_neg xdata, ydata

	ratio_true = 	interp_pos /  interp_neg	
			
	Interpolate2 /i=3 /X = posn_pos /Y=interp_pos xdata, curve
	Interpolate2 /i=3 /X = posn_neg /Y=interp_neg xdata, curve
				
	ratio_curve  = 	interp_pos /  interp_neg	
	ratio_et =  ratio_curve / ratio_true
	
	printf "True y values = ydata\r"
	printf "Perturbed wave = curved\r"
	printf "Perturbation =  quadratic polynomial : 1 + (0.014) x + (-0.005)  x^2\r" 
	
end	

//******************************************************************************************************
//******************************************************************************************************

//	This is the function describing the residual when two parameters (for polynomial) are used as 
//		variables

function  residual (w,  k1, k2)
	wave w			// optional wave
	variable  k1		// fit parameter
	variable k2		// fit  parameter
	
	make  /o /FREE /d /n=(7) rhs_n	// numerator, RHS
	make  /o /FREE /d /n=(7) rhs_d	//  denominator, RHS
	make  /o /FREE /d /n=(7) resd	// residual, wave
	
	wave posn_pos
	wave posn_neg
	wave ratio_et
	
	rhs_n = 1+k1/10*posn_pos + k2/100*(posn_pos^2)	// scaling k1 and  k2 reduces numerical error
	rhs_d = 1+k1/10*posn_neg + k2/100*(posn_neg^2)	// scaling k1 and  k2 reduces numerical error
	
	resd = ( ratio_et - (rhs_n / rhs_d))^2
	
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

function run_Opt ( temp_wave  , initial_k1,  initial_k2 ) 
	wave temp_wave	// optional wave passed to the residual function
	variable  initial_k1	// guess for k1
	variable  initial_k2	// guess for k2

	
	Optimize /A=0  /I=3000 /M={0 , 0} /T={1e-9, 1e-11}  /R={1e-2, -1e-2}  /X={ initial_k1 , initial_k2 }  residual , temp_wave  
	//			     iterations  method  tolerance          order of params         initial coefs                    function    optional
	
	wave  result = W_Extremum
	wave xdata
	wave curve

	// k1 was scaled by  10, k2 was scaled by 100
	printf "Optimized parameters : %5.6f, %5.6f\r", result[0] / 10, result[1] / 100
	
	make /o /d /n=100 corn_curve = 1+ ( result[0]/10) * xdata + (result[1]/100)*xdata^2 
	
	make  /o  /d /n=100 corrected  = curve / corn_curve		//  generate the  corrected wave, compare to the true ydata
end

//******************************************************************************************************
//******************************************************************************************************

