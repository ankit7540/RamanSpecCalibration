
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//*******************************************************************************************************************
// Define constants here : 
constant scale = 1e3
//*******************************************************************************************************************

//  Structure of the data wave :

//  Data waves  are  required  as 2D waves with columns  containing the data as following.

//	|   J   |  expt  ratio |  theoretical ratio |   J+2  transition wavenumber |  J-2  transition wavenumber | weight  | 

// For example,   (header  given for  clarity)

// | J | Expt Ratio | Theortcl  Ratio | J+2, freq | J-2, freq | Weight |
// |---|     ------------|-----------------      |-----------   |-----------   |--------    |
// | 2 | 2.4872     | 2.3827               | 414.65    | -179.07   | 11.67   |
// | 3 | 1.7910     | 1.6652               | 529.90    | -297.53   | 9.48    |
// | 4 | 1.4674     | 1.3889               | 642.81    | -414.65   | 10.72  |
// | 5 | 1.3904     | 1.2331               | 752.92    | -529.90   | 4.49    |
// | 6 | 1.2955     | 1.1287               | 859.83    | -642.81   | 3.26    |



//*******************************************************************************************************************

// Error function describing the difference between the experimental values and
//     theoretical data  utilising  a quadratic polynomial as  the wavelength
//     dependent sensitivity

//	All gases included

function Residual_Linear  (w, k1 )
	
	wave w		// optional  wave, not used here
	variable k1	// coef of  the quadratic function
	variable k2	// // coef of  the quadratic function
	
	// Linear function used  : 1+ k1*x  
	
	 // ---------------------- waves containing  data   for  fit------------------------------
	 //  Modify the location of  the  data waves	 
	wave dataH2   = root:dataH2
	wave dataD2  = root:dataD2
	wave dataHD  = root:dataHD
	wave dataO2_pR  = root:DataO2_p
	wave dataO2_o1s1		= root:DataO2
	// -----------------------------------------------------------------------------------------------		

	
	// define number of data points for each gas
	variable nH2	=	dimsize (dataH2 ,  0)
	variable nHD	=	dimsize (dataHD ,  0)
	variable nD2	=	dimsize (dataD2 ,  0)
	variable nO2_pR	=	dimsize (dataO2_pR ,  0)
	variable nO2_o1s1	=	dimsize (dataO2_o1s1 ,  0)
					
			
	// temporary waves for data 			
	make /o /d /FREE   /n=(nH2)  ratioH2, RHS_H2, resdH2
	make /o /d /FREE   /n=(nHD)  ratioHD, RHS_HD, resdHD
	make /o /d /FREE   /n=(nD2)  ratioD2, RHS_D2, resdD2
	make /o /d /FREE   /n=(nO2_pR)  ratioO2_pR, RHS_O2pR, resdO2_pR	
	make /o /d /FREE   /n=(nO2_o1s1)  ratioO2_o1s1, RHS_O2_o1s1, resdO2_o1s1

	// H2
	ratioH2[]  = dataH2[p][1] / dataH2[p][2]
	RHS_H2 []  =  ( 1+ k1* (dataH2[p][3] /  scale )  ) / ( 1+ k1* (dataH2[p][4] /  scale )  )
	resdH2[] =   dataH2[p][5]  * (ratioH2[p] - RHS_H2 [p] )^2

	// HD
	ratioHD[]  = dataHD[p][1] / dataHD[p][2]
	RHS_HD []  =  ( 1+ k1* (dataHD[p][3] /  scale )  ) / ( 1+ k1* (dataHD[p][4] /  scale )  )
	resdHD[] =   dataHD[p][5]  * (ratioHD[p] - RHS_HD [p] )^2
	
	//D2
	ratioD2[]  = dataD2[p][1] / dataD2[p][2]
	RHS_D2 []  =  ( 1+ k1* (dataD2[p][3] /  scale )  ) / ( 1+ k1* (dataD2[p][4] /  scale )  )
	resdD2[] =   dataD2[p][5]  * (ratioD2[p] - RHS_D2 [p] )^2
	
	//O2_pureRotation
	ratioO2_pR[]  = dataO2_pR[p][1]/dataO2_pR[p][2]
	RHS_O2pR []  =  ( 1+ k1* (dataO2_pR[p][3] /  scale )  ) / ( 1+ k1* (dataO2_pR[p][4] /  scale )  )
	resdO2_pR[] =   dataO2_pR[p][5]  * (ratioO2_pR[p] - RHS_O2pR [p] )^2
		
	
	//O2_O1S1
	ratioO2_o1s1[]  = dataO2_o1s1[p][1]/dataO2_o1s1[p][2]
	RHS_O2_o1s1 []  =  ( 1+ k1* (dataO2_o1s1[p][3] /  scale )  ) /  ( 1+ k1* (dataO2_o1s1[p][4] /  scale )  )
	resdO2_o1s1[] =   dataO2_o1s1[p][5]  * (ratioO2_o1s1[p] - RHS_O2_o1s1 [p] )^2
		
	return  WaveSum( resdH2 )  + WaveSum( resdHD ) + WaveSum( resdD2 )  + WaveSum( resdO2_pR )  + WaveSum(  resdO2_o1s1 )
	
end	


//*******************************************************************************************************************
//*******************************************************************************************************************

// Error function describing the difference between the experimental values and
//     theoretical data  utilising  a quadratic polynomial as  the wavelength
//     dependent sensitivity

//	All gases included

function Residual_Quadratic  (w, k1, k2)
	
	wave w		// optional  wave, not used here
	variable k1	// coef of  the quadratic function
	variable k2	// // coef of  the quadratic function
	
	// Quadratic function used  : 1+ k1*x + k2*x^2
	
	 // ---------------------- waves containing  data   for  fit------------------------------
	 //  Modify the location of  the  data waves
	wave dataH2   = root:dataH2
	wave dataD2  = root:dataD2
	wave dataHD  = root:dataHD
	wave dataO2_pR  = root:DataO2_p
	wave dataO2_o1s1		= root:DataO2
	// -----------------------------------------------------------------------------------------------		

	
	// define number of data points for each gas
	variable nH2	=	dimsize (dataH2 ,  0)
	variable nHD	=	dimsize (dataHD ,  0)
	variable nD2	=	dimsize (dataD2 ,  0)
	variable nO2_pR	=	dimsize (dataO2_pR ,  0)
	variable nO2_o1s1	=	dimsize (dataO2_o1s1 ,  0)
					
			
	// temporary waves for data 			
	make /o /d /FREE   /n=(nH2)  ratioH2, RHS_H2, resdH2
	make /o /d /FREE   /n=(nHD)  ratioHD, RHS_HD, resdHD
	make /o /d /FREE   /n=(nD2)  ratioD2, RHS_D2, resdD2
	make /o /d /FREE   /n=(nO2_pR)  ratioO2_pR, RHS_O2pR, resdO2_pR	
	make /o /d /FREE   /n=(nO2_o1s1)  ratioO2_o1s1, RHS_O2_o1s1, resdO2_o1s1

	// H2
	ratioH2[]  = dataH2[p][1] / dataH2[p][2]
	RHS_H2 []  =  ( 1+ k1* (dataH2[p][3] /  scale ) + k2*( (dataH2[p][3] /scale)^2)   ) / ( 1+ k1* (dataH2[p][4] /  scale ) + k2*( (dataH2[p][4] /scale)^2)   )
	resdH2[] =   dataH2[p][5]  * (ratioH2[p] - RHS_H2 [p] )^2

	// HD
	ratioHD[]  = dataHD[p][1] / dataHD[p][2]
	RHS_HD []  =  ( 1+ k1* (dataHD[p][3] /  scale ) + k2*( (dataHD[p][3] /scale)^2) ) / ( 1+ k1* (dataHD[p][4] /  scale ) + k2*( (dataHD[p][4] /scale)^2) )
	resdHD[] =   dataHD[p][5]  * (ratioHD[p] - RHS_HD [p] )^2
	
	//D2
	ratioD2[]  = dataD2[p][1] / dataD2[p][2]
	RHS_D2 []  =  ( 1+ k1* (dataD2[p][3] /  scale ) + k2*( (dataD2[p][3] /scale)^2) ) / ( 1+ k1* (dataD2[p][4] /  scale ) + k2*( (dataD2[p][4] /scale)^2) )
	resdD2[] =   dataD2[p][5]  * (ratioD2[p] - RHS_D2 [p] )^2
	
	//O2_pureRotation
	ratioO2_pR[]  = dataO2_pR[p][1]/dataO2_pR[p][2]
	RHS_O2pR []  =  ( 1+ k1* (dataO2_pR[p][3] /  scale ) + k2*( (dataO2_pR[p][3] /scale)^2) ) / ( 1+ k1* (dataO2_pR[p][4] /  scale ) + k2*( (dataO2_pR[p][4] /scale)^2) )
	resdO2_pR[] =   dataO2_pR[p][5]  * (ratioO2_pR[p] - RHS_O2pR [p] )^2
		
	
	//O2_O1S1
	ratioO2_o1s1[]  = dataO2_o1s1[p][1]/dataO2_o1s1[p][2]
	RHS_O2_o1s1 []  =  ( 1+ k1* (dataO2_o1s1[p][3] /  scale ) + k2*( (dataO2_o1s1[p][3] /scale)^2) ) /  ( 1+ k1* (dataO2_o1s1[p][4] /  scale ) + k2*( (dataO2_o1s1[p][4] /scale)^2) )
	resdO2_o1s1[] =   dataO2_o1s1[p][5]  * (ratioO2_o1s1[p] - RHS_O2_o1s1 [p] )^2
		
	return  WaveSum( resdH2 )  + WaveSum( resdHD ) + WaveSum( resdD2 )  + WaveSum( resdO2_pR )  + WaveSum(  resdO2_o1s1 )
	
end	

//*******************************************************************************************************************
//*******************************************************************************************************************

// Error function describing the difference between the experimental values and theoretical data
//	along with some function describing the response function

//	All gases included

function Residual_Cubic ( w,  k1,  k2,  k3)
	
	wave w
	variable k1
	variable k2
	variable k3
	
	// Cubic function used  : 1+ k1*x + k2*x^2 + k3*x^3
	
	 // ---------------------- waves containing  data   for  fit------------------------------
	 //  Modify the location of  the  data waves	 
	wave dataH2   = root:dataH2
	wave dataD2  = root:dataD2
	wave dataHD  = root:dataHD
	wave dataO2_pR  = root:DataO2_p
	wave dataO2_o1s1		= root:DataO2
	// -----------------------------------------------------------------------------------------------		

	
	// define number of data points for each gas
	variable nH2	=	dimsize (dataH2 ,  0)
	variable nHD	=	dimsize (dataHD ,  0)
	variable nD2	=	dimsize (dataD2 ,  0)
	variable nO2_pR	=	dimsize (dataO2_pR ,  0)
	variable nO2_o1s1	=	dimsize (dataO2_o1s1 ,  0)
					
			
	// temporary waves for data 			
	make /o /d /FREE   /n=(nH2)  ratioH2, RHS_H2, resdH2
	make /o /d /FREE   /n=(nHD)  ratioHD, RHS_HD, resdHD
	make /o /d /FREE   /n=(nD2)  ratioD2, RHS_D2, resdD2
	make /o /d /FREE   /n=(nO2_pR)  ratioO2_pR, RHS_O2pR, resdO2_pR	
	make /o /d /FREE   /n=(nO2_o1s1)  ratioO2_o1s1, RHS_O2_o1s1, resdO2_o1s1

	// H2
	ratioH2[]  = dataH2[p][1] / dataH2[p][2]
	RHS_H2 []  =  ( 1+ k1* (dataH2[p][3] /  scale ) + k2*( (dataH2[p][3] /scale)^2) + k3*( (dataH2[p][3] /scale)^3) ) / ( 1+ k1* (dataH2[p][4] /  scale ) + k2*( (dataH2[p][4] /scale)^2)  + k3*( (dataH2[p][4] /scale)^3)   )
	resdH2[] =   dataH2[p][5]  * (ratioH2[p] - RHS_H2 [p] )^2

	// HD
	ratioHD[]  = dataHD[p][1] / dataHD[p][2]
	RHS_HD []  =  ( 1+ k1* (dataHD[p][3] /  scale ) + k2*( (dataHD[p][3] /scale)^2) + k3*( (dataHD[p][3] /scale)^3) ) / ( 1+ k1* (dataHD[p][4] /  scale ) + k2*( (dataHD[p][4] /scale)^2)+ k3*( (dataHD[p][4] /scale)^3) )
	resdHD[] =   dataHD[p][5]  * (ratioHD[p] - RHS_HD [p] )^2
	
	//D2
	ratioD2[]  = dataD2[p][1] / dataD2[p][2]
	RHS_D2 []  =  ( 1+ k1* (dataD2[p][3] /  scale ) + k2*( (dataD2[p][3] /scale)^2) +k3*( (dataD2[p][3] /scale)^3) ) / ( 1+ k1* (dataD2[p][4] /  scale ) + k2*( (dataD2[p][4] /scale)^2)+ k3*( (dataD2[p][4] /scale)^3) )
	resdD2[] =   dataD2[p][5]  * (ratioD2[p] - RHS_D2 [p] )^2
	
	//O2_pureRotation
	ratioO2_pR[]  = dataO2_pR[p][1]/dataO2_pR[p][2]
	RHS_O2pR []  =  ( 1+ k1* (dataO2_pR[p][3] /  scale ) + k2*( (dataO2_pR[p][3] /scale)^2) +k3*( (dataO2_pR[p][3] /scale)^3) ) / ( 1+ k1* (dataO2_pR[p][4] /  scale ) + k2*( (dataO2_pR[p][4] /scale)^2) + k3*( (dataO2_pR[p][4] /scale)^3) )
	resdO2_pR[] =   dataO2_pR[p][5]  * (ratioO2_pR[p] - RHS_O2pR [p] )^2
		
	
	//O2_O1S1
	ratioO2_o1s1[]  = dataO2_o1s1[p][1]/dataO2_o1s1[p][2]
	RHS_O2_o1s1 []  =  ( 1+ k1* (dataO2_o1s1[p][3] /  scale ) + k2*( (dataO2_o1s1[p][3] /scale)^2)  + k3*( (dataO2_o1s1[p][3] /scale)^3) ) /  ( 1+ k1* (dataO2_o1s1[p][4] /  scale ) + k2*( (dataO2_o1s1[p][4] /scale)^2) + k3*( (dataO2_o1s1[p][4] /scale)^3) )
	resdO2_o1s1[] =   dataO2_o1s1[p][5]  * (ratioO2_o1s1[p] - RHS_O2_o1s1 [p] )^2
	
	return  WaveSum( resdH2 )  + WaveSum( resdHD ) + WaveSum( resdD2 )  + WaveSum( resdO2_pR )  + WaveSum(  resdO2_o1s1 )
	
end	


//*******************************************************************************************************************
//*******************************************************************************************************************

// This  function runs the non-linear optimization process using the 'Optimize' function in  IgorPro

//	Example :
//	make  /n=1 temp=0 // this wave is passed to the optimize function, not used in present implementation
//	run_Opt_Linear  ( temp   , 0.115   ) 

function run_Opt_Linear  ( temp_wave  , initial_k1 ) 
	wave temp_wave	// optional wave passed to the residual function
	variable  initial_k1	// guess for k1
	variable  initial_k2	// guess for k2

	
	Optimize /A=0  /I=3000 /M={0 , 0} /T={1e-9, 1e-11}    /X={ initial_k1  }  Residual_Linear , temp_wave  
	//			     iterations  method  tolerance                  initial coefs                    function                       optional
	
	wave  result = W_Extremum

	printf "Obtained  coef for the  linear function : %5.6e \r", result[0] / scale 

end

//*******************************************************************************************************************
//*******************************************************************************************************************

// This  function runs the non-linear optimization process using the 'Optimize' function in  IgorPro

//	Example :
//	make  /n=1 temp=0 // this wave is passed to the optimize function, not used in present implementation
//	run_Opt_Quadratic  ( temp   , 0.115 , 1e-2  ) 

function run_Opt_Quadratic  ( temp_wave  , initial_k1,  initial_k2 ) 
	wave temp_wave	// optional wave passed to the residual function
	variable  initial_k1	// guess for k1
	variable  initial_k2	// guess for k2

	
	Optimize /A=0  /I=3000 /M={0 , 0} /T={1e-9, 1e-11}    /X={ initial_k1 , initial_k2 }  Residual_Quadratic , temp_wave  
	//			     iterations  method  tolerance                  initial coefs                    function                       optional
	
	wave  result = W_Extremum

	printf "Obtained  coefs for the quadratic polynomial : %5.6e, %5.6e\r", result[0] / scale, result[1] / (scale^2)

end

//*******************************************************************************************************************
//*******************************************************************************************************************

// This  function runs the non-linear optimization process using the 'Optimize' function in  IgorPro

//	Example :
//	make  /n=1 temp=0 // this wave is passed to the optimize function, not used in present implementation
//	run_Opt_Cubic  ( temp   , 0.115 , 1e-2,  1e-4 ) 


function run_Opt_Cubic  ( temp_wave  , initial_k1,  initial_k2,  initial_k3 ) 
	wave temp_wave	// optional wave passed to the residual function
	variable  initial_k1	// guess for k1
	variable  initial_k2	// guess for k2
	variable  initial_k3	// guess for k3

	
	Optimize /A=0  /I=3000 /M={0 , 0} /T={1e-9, 1e-11}    /X={ initial_k1 , initial_k2,  initial_k3 }  Residual_Cubic , temp_wave  
	//			     iterations  method  tolerance                  initial coefs                    function                       optional
	
	wave  result = W_Extremum

	printf "Obtained  coefs for the cubic polynomial : %5.6e, %5.6e, %5.6e\r", result[0] / scale, result[1] / (scale^2), result[ 2 ] / (scale^3)

end

//*******************************************************************************************************************
//*******************************************************************************************************************

// WaveSum(w)
// Returns the sum of the entire wave, just like Igor’s sum function.

//	Parameters:
//		w	=	1D wave
Function  WaveSum(w)
	Wave w
	Variable i, n=numpnts(w), total=0
	
	for(i=0;i<n;i+=1)
		total += w[i]
	endfor
	
	return total
End

//*******************************************************************************************************************
//*******************************************************************************************************************
