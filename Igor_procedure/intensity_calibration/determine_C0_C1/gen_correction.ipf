

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------

// Main function
// Example :
//    gen_C0_C1 ( x_wavenumber, 632.8 , root:whitelight  ,   670 , maskWave = mask1D  )
//                  wave        laser_nm   wave             pixel   maskWave

// returns :  intensity_corr : 1D wave which corresponds to (C0/C1), to be multiplied to the Raman spectra.

function gen_C0_C1 ( ramanshift , laser_nm , wl_wave ,   norm_pnt , [ maskWave , set_mask_nan] )

	wave ramanshift		// relative, vector
	variable laser_nm		// laser wavelength, nm, scalar
	wave wl_wave			// white_light spectra, un-normalized, 1D or 2D, vector
	variable norm_pnt		// normalization index, pnt, scalar
	
	wave maskWave //  optional, supply mask for wl fitting , for eg.  maskWave = mask_1D, vector
	variable  set_mask_nan //  optional, 0 or 1. When set to 1, the masked region is set to nan in the output 
							//    intensity_corr wave. When set to 0, intensity_corr wave is not modified for the mask
							//    region.
	
	variable nPnts = dimsize(ramanshift, 0)
	variable mid = nPnts/2
	
	
	// initialization
	killwaves /Z  W_coef, W_ParamConfidenceInterval, W_sigma, fit
	
	// generate the C0 correction
	gen_C0(ramanshift ,   norm_pnt)
	wave C0
	
	// wl spectra (check averaging)
	variable nCols_wl = dimsize(wl_wave,1)
	if (nCols_wl > 1)
		print "\twl_wave is 2D. Averaging.\r"
		do_average_custom(wl_wave)
		string newWave= nameofwave(wl_wave)+"_avg"
		printf "\tAveraged wave : %s\r", newWave
		wave wl_wave = $newWave
	endif	
	
	// wl spectra (normalization) and multiply with C0
	make /d /o /n=(nPnts) wl_norm = wl_wave/ wavemax(wl_wave)
	wl_norm=wl_norm*C0
	wave wl_norm
	
	// wl spectra (subset)
	variable result  // returned value is assigned to result
	
	// Determine if the optional parameters were supplied
	if( ParamIsDefault(maskWave) ) 
		
		print "\t Mask wave not supplied. Working with full wave.\r"
		result = gen_C1(ramanshift, laser_nm ,  wl_norm ,  norm_pnt)

	elseif ( waveExists ( maskWave ))
		print "\tMask wave is supplied. Using it for fit.\r"
		result = gen_C1_mask (ramanshift, laser_nm ,  wl_norm ,   norm_pnt , maskWave)
	endif	
	
	
	//--------------------------------------------------------------
	// wave references (C0 and C1 should exist in present folder)
	wave C0
	wave C1
	//--------------------------------------------------------------
	
	if (result ==2 )
	
		if (set_mask_nan == 1)
			printf "\n\n\tset_mask_nan parameter is set to 1. Masked region in the intensity corr will be set to Nan.\r"
			
			wave C1_out = set_nan_correspond_to_mask (C1, maskWave)
			wave C1= C1_out
		else
			printf "\n\n\tset_mask_nan parameter not set to 1.\r"
		endif	
	endif	
	//--------------------------------------------------------------	
	make /d /o /n=(nPnts) intensity_corr = (C0/C1)
	
	printf "\n\tintensity_corr generated. Done. This is equal to (C0/C1)\r"
	
	// cleanup
	killwaves  /Z abs_wavenumber_subset, wl_norm_subset,  W_coef, fit, coefs
	killwaves /Z W_ParamConfidenceInterval, W_sigma

end


// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------

STATIC function gen_C0(ramanshift ,   norm_pnt)

	wave ramanshift		// relative
	variable norm_pnt		// normalization index, pnt
	
	variable nPnts = dimsize(ramanshift, 0)
	variable mid = nPnts/2
	
	make   /FREE /o /n=(nPnts-1) wavenumber_spacing =0
	wavenumber_spacing = ramanshift [p] -  ramanshift [p+1]	
	
	// normalization
	make /d /o /n=(nPnts-1) waveum_corr = ( wavenumber_spacing / wavenumber_spacing [ norm_pnt ])
	make /FREE /d /n=(nPnts-1) xax=p
	// ------------------------	
	
	// computing the value at the last point
	CurveFit /Q  poly 3,  waveum_corr /X=xax  
	wave W_coef
	variable val_plus1_pnt =   poly(W_coef, nPnts-1 ) 
	//print val_plus1_pnt
	
	InsertPoints (nPnts-1),1, waveum_corr
	waveum_corr [nPnts-1] = val_plus1_pnt
	
	duplicate /O waveum_corr, C0
	
	killwaves /Z wavenum_corr


end

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------

STATIC function gen_C1(ramanshift, laser_nm ,  wl_input   norm_pnt)

	wave ramanshift		// relative
	variable laser_nm		// laser wavelength (nm)
	wave wl_input		// normalized wl wave with C0 multiplied
	variable norm_pnt		// normalization index, pnt
	
	variable nPnts = dimsize(ramanshift, 0)
	variable mid = nPnts/2
	

	// cleanup from prev run
	killwindow /Z genC0C1	
	
	// make abs_wavenumber, for xaxis
	
	make /o /d /n=(nPnts) abs_wavenumber = ((1e7/laser_nm)-ramanshift)*100
	make /o /d /n=2 coefs
	
	// defining initial coefs
	coefs[0] = 0.856e-18
	coefs[1] = 2759 // 0.856e-18
	
		
	Display /K=1 wl_input vs abs_wavenumber
	FuncFit photons_per_unit_wavenum_abs coefs wl_input /X=abs_wavenumber 
	
	printf "\tFit (not using mask), Obtained fit coef:%5.6e, %5.6e",coefs[0],coefs[1]
	
	make /o /d /n=(nPnts) fit = photons_per_unit_wavenum_abs(coefs,abs_wavenumber)
	appendtograph fit vs abs_wavenumber
	
	// customize graph
	ModifyGraph  grid=1,mirror=1,axThick=1.2,lblPosMode(bottom)=2;DelayUpdate
	ModifyGraph  lblMargin(bottom)=2,gridHair=1,manTick(left)={0,0.1,0,1},manMinor(left)={3,2},manTick(bottom)={0,0.02,6,2}
	ModifyGraph    manMinor(bottom)={3,2},gridRGB(left)=(13107,13107,13107),gridRGB(bottom)=(17476,17476,17476); 
	Label  left "Intensity";DelayUpdate
	Label  bottom "Wavenumber / m\\S-1"
	ModifyGraph  fSize=25
	ModifyGraph  rgb(fit)=(0,0,65535)
	ModifyGraph  lsize=2
	SetAxis   left 0,*
	DoUpdate;		
	
	// generate C1
	make /o /d /n=(nPnts) C1 = wl_input / fit

	// remove fit and residual  wave
	string fname="fit_"+nameofwave(wl_input)
	killwaves /Z  $fname	

	fname="Res_"+nameofwave(wl_input)
	killwaves /Z  $fname	

	
	if (waveexists(C1) && coefs[1]< 6000 )   // 6000 is the temperature of the lamp in K, too large value may be unreasonable.
		return 1
	endif		


end
	 
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------	

STATIC function gen_C1_mask (ramanshift, laser_nm ,  wl_input  , norm_pnt , maskWave )

	wave ramanshift		// relative
	variable laser_nm		// laser wavelength (nm)
	wave wl_input		// normalized wl wave with C0 multiplied
	variable norm_pnt		// normalization index, pnt
	wave maskWave
	
	variable nPnts = dimsize(ramanshift, 0)
	variable mid = nPnts/2
	
	// cleanup from prev run
	killwindow /Z genC0C1
	
	// make abs_wavenumber, for xaxis
	
	make /o /d /n=(nPnts) abs_wavenumber = ((1e7/laser_nm)-ramanshift)*100
	make /o /d /n=2 coefs
	
	// defining initial coefs
	coefs[0] = 0.856e-18
	coefs[1] = 2659 // 0.856e-18 
	
	
	Display /N=genC0C1 /K=1 wl_input vs abs_wavenumber

	// perform fitting, with mask applied
	FuncFit photons_per_unit_wavenum_abs coefs wl_input /X=abs_wavenumber /M=maskWave
	
	printf "\tMasked wl fit, obtained coefs : %5.6e, %5.6e",coefs[0],coefs[1]
	
	// generate the fit curve
	make /o /d /n=(nPnts) fit = photons_per_unit_wavenum_abs(coefs,abs_wavenumber)

	// generate C1
	make /o /d /n=(nPnts) C1 = wl_input / fit
	
	//customize the graph
	appendtograph /W=genC0C1 fit vs abs_wavenumber
	
	ModifyGraph /W=genC0C1 grid=1,mirror=1,axThick=1.2,lblPosMode(bottom)=2;DelayUpdate
	ModifyGraph /W=genC0C1 lblMargin(bottom)=2,gridHair=1,manTick(left)={0,0.1,0,1},manMinor(left)={3,2},manTick(bottom)={0,0.02,6,2}
	ModifyGraph  /W=genC0C1  manMinor(bottom)={3,2},gridRGB(left)=(13107,13107,13107),gridRGB(bottom)=(17476,17476,17476); 
	Label /W=genC0C1 left "Intensity";DelayUpdate
	Label /W=genC0C1 bottom "Wavenumber / m\\S-1"
	ModifyGraph /W=genC0C1 fSize=25
	ModifyGraph /W=genC0C1 rgb(fit)=(0,0,65535)
	ModifyGraph /W=genC0C1 lsize=2
	SetAxis /W=genC0C1  left 0,*
	DoUpdate;
	
	// generate C1
	make /o /d /n=(nPnts) C1 = wl_input / fit
	
	// remove fit wave
	string fname="fit_"+nameofwave(wl_input)
	killwaves /Z  $fname

	fname="Res_"+nameofwave(wl_input)
	killwaves /Z  $fname		
	
	if (waveexists(C1) && coefs[1]< 6000 ) // 6000 is the temperature of the lamp in K, too large value may be unreasonable.
		return 2
	endif	
	

end
	 

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------	 
// Custom fit function 

Function photons_per_unit_wavenum_abs(w0,w) : FitFunc
	Wave w0
	Variable w

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(w) = a0*599584916*w^2/(exp(0.1438776877e-1*w/T)-1)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ w
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w0[0] = a0
	//CurveFitDialog/ w0[1] = T

	return w0[0]*599584916*w^2/(exp(0.1438776877e-1*w/w0[1])-1)
End

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------

//  for averaging 2D input wave across columns	
STATIC function do_average_custom(input)
	wave input
	
	variable nRows
	variable nCols
	variable i
	
	nRows=dimsize(input, 0)
	nCols=dimsize(input, 1)	

	
	string name=nameofwave (input)
	
	string new_name
	sprintf new_name, "%s_avg",name
	printf "\t%s\t%g\t%g\r" new_name, nRows, nCols
	
	make /o/d /n=(nRows) $new_name
	wave output=$new_name
	
	make /o/d /n=(nRows) temp=0
	
	for (i=0 ; i<nCols ;  i=i+1)
		temp=temp+input [p][(i)]
	endfor
	
	output=temp / nCols
	
	killwaves /Z temp
	
end

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------	

// set points in input 1D wavto nan, for points where mask has values == 0
// returns wave reference
STATIC function /WAVE set_nan_correspond_to_mask (input, mask)
	wave input
	wave mask
	
	variable nRows=dimsize (input, 0)
	variable i

	for (i=0 ; i<nRows ; i=i+1)
		if (mask[i]==0)
			input[i] = nan		
		endif

	endfor	
	return input
end	
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
