
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//***************************************************************************************
//***************************************************************************************


//  Fitting of  the reference wavenumber vs  the  band  position in pixel
//	including the error in both y and  x

function  wavenumber_fit ( ref_wave, ref_error , pixel_posn , pixel_error ,  nPoints)
	wave ref_wave		//	reference wavenumber
	wave ref_error 		//	error in the  reference wavenumber
	wave pixel_posn 	//	band  position in  pixel
	wave pixel_error 	//	error in  the band  position in pixel
	variable  nPoints	//	Number of points on  the CCD

	variable  nRows

	nRows  = dimsize (pixel_posn, 0)

	make /d  /o /n=(nRows) residual_y, residual_x

	// Change the initial coefs as required
	make /d /o /n=4 initCoef={  -1045  ,  1.5  , -1e-4  , 1e-008  }

	CurveFit /ODR=2   poly 4 ,  kwCWave=initCoef  , ref_wave /X=pixel_posn  /W=ref_error   /XW=pixel_error   /I=1 /D /R=residual_y  /XR=residual_x

  // fit  coefs  are  updated in the initCoef wave
	make /d /o /n=(nPoints)   fit  = poly( initCoef , p )

	make  /d /o /n=(nPoints)   pixel_x = p


	variable yr_min = abs(  wavemin (residual_y) )
	variable yr_max= abs (wavemax (residual_y))
	variable xr_min = abs(wavemin (residual_x))
	variable xr_max= abs(wavemax (residual_x))

	// print  xr_min, xr_max,  yr_min, yr_max

	// limits for the residuals in x and y --------
	variable yResd_limit
	variable xResd_limit

	if (yr_min < yr_max)
		yResd_limit = yr_max
	else
		yResd_limit = yr_min
	endif

	if (xr_min < xr_max)
		xResd_limit = xr_max
	else
		xResd_limit = xr_min
	endif
	// -----------------------------------------------------
	// print  xResd_limit, yResd_limit


	display ref_wave vs  pixel_posn		// Plot the fit and residuals
	AppendToGraph/L=l1 residual_y vs pixel_posn;
	AppendToGraph/R=r1 residual_x vs pixel_posn;


	ModifyGraph axisEnab(left)={0,0.78}
	ModifyGraph axisEnab(l1)={0.82,1}
	ModifyGraph  axisEnab(r1)={0.82,1}
	ModifyGraph freePos(l1)=0
	ModifyGraph freePos(r1)=0

	string  wname=nameofwave(ref_wave)
	ModifyGraph mode=3,marker=19,msize=3
	ModifyGraph rgb($wname)=(65535,0,0)
	ModifyGraph useMrkStrokeRGB($wname)=1
	ModifyGraph rgb(residual_y)=(16385,28398,65535)

	SetAxis l1 (-1*yResd_limit),  yResd_limit
	SetAxis r1 (-1*xResd_limit),  xResd_limit

	ModifyGraph mirror(left)=1
	ModifyGraph mirror(bottom)=1
	ModifyGraph axThick(left)=1.2
	ModifyGraph axThick(bottom)=1.2
	ModifyGraph manTick(left)={0,500,0,0}
	ModifyGraph manMinor(left)={3,2}
	ModifyGraph manTick(bottom)={0,200,0,0}
	ModifyGraph manMinor(bottom)={3,2}
	ModifyGraph manTick(l1)={0,20,-3,0}
	ModifyGraph manMinor(l1)={0,50}
	ModifyGraph manTick(r1)={0,0.5,0,1}
	ModifyGraph manMinor(r1)={3,2}

	Label left "Wavenumber  / cm\\S-1";DelayUpdate
	Label bottom "Pixel index";DelayUpdate
	Label l1 "Residual y";DelayUpdate
	Label r1 "Residual x"
	ModifyGraph lblPosMode(l1)=2,lblPosMode(r1)=2
	ModifyGraph lblMargin(l1)=12,lblMargin(r1)=12
	ModifyGraph axThick=1.2

	ModifyGraph grid(left)=1,grid(bottom)=1,grid(l1)=1
	ModifyGraph  gridRGB(left)=(4369,4369,4369)
	ModifyGraph gridRGB(bottom)=(4369,4369,4369)
	ModifyGraph gridRGB(l1)=(4369,4369,4369)


	AppendToGraph fit vs pixel_x
	ModifyGraph lsize(fit)=2,rgb(fit)=(0,0,0)
	ModifyGraph lsize(fit)=1.5

	ErrorBars freq XY,wave=( pixel_error , pixel_error ),wave=( ref_error, ref_error )
	ErrorBars/T=2/L=2/X=10/Y=10 freq XY,wave=(  pixel_error , pixel_error   ),wave=(  ref_error, ref_error  )

	ModifyGraph useMrkStrokeRGB(residual_y)=1,useMrkStrokeRGB(residual_x)=1
	ModifyGraph rgb(residual_x)=(65535,43688,32768),rgb(fit)=(4369,4369,4369)
	ModifyGraph rgb(freq)=(65535,0,0)
	ModifyGraph axRGB(l1)=(1,16019,65535),axRGB(r1)=(65535,21845,0)
	ModifyGraph tlblRGB(l1)=(1,16019,65535),tlblRGB(r1)=(65535,21845,0)
	ModifyGraph alblRGB(l1)=(1,16019,65535),alblRGB(r1)=(65535,21845,0)
	Legend/C/N=text1/A=RB
	ModifyGraph manTick(r1)={0,0.1,0,1},manMinor(r1)={0,50}
	ModifyGraph zero(l1)=1,axRGB(l1)=(0,0,0)

	duplicate  /o  fit,  Wavenumber_axis

	printf  "Wavenumber_axis wave  generated.\r\n"
end
//***************************************************************************************
//***************************************************************************************
