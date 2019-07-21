
// Functions for computing the pure rotational Raman intensities of H2, HD and D2.

// set params for the rotational intensity analysis.----------------------------

	function rotational_analysis()

			string cdf=getdatafolder(1)
			
			string preset_folder="root:Params:RotnRamanIns:"
			//	check if exists
			if (datafolderexists(preset_folder)==1)
				// print "folder exists"
				setdatafolder root:Params:RotnRamanIns:
			else
				newdatafolder /O/S root:Params
				newdatafolder /O/S RotnRamanIns
				newdatafolder /O  opt
			endif

			// make waves and set variables.
			variable /G temp= 298.0
			variable /G wavenumber=18789.9380
		
		 print "Folder (root:Params:RotnRamanIns)  generated. Please load the data files in the following folder. \r"
		 print "root:Params:RotnRamanIns  \r"
		setdatafolder $cdf
			
	end

//**********************************************************************************************************
//	Define constants

	constant k_PlanckConstant = 1.38064852e-23 // J/K
	constant k_h = 6.626070040e-34 // J.s
	constant k_c = 2.99792458e+10 //cm/s

//**********************************************************************************************************
//**********************************************************************************************************

// compute the sum of states for H2 at a given temperature for upto Jf state.

function SumOfstates_H2 (T)
	variable T
	// 	Argument
	//	T	: temperature, Kelvin

	variable  ge, go
	variable q1, su
	variable i1, Ji, E1
	variable jMaxv0
	variable jMaxv1

	// --- Nuclear spin statictics ------------------------------------------------------
	ge=1 	// molecular hydrogen
	go=3
	//---------------------------------------------------------------------------------------

	wave ejv0 = root:Params:RotnRamanIns:H2eV0
	wave ejv1 = root:Params:RotnRamanIns:H2eV1

	jMaxv0 = dimsize(ejv0, 0)
	jMaxv1 = dimsize(ejv1, 0)

	make /o /d /n = (jMaxv0 + jMaxv1) factor
	wave fac = factor 

	// ------------------------------------------------------------

	// ground vibrational state
	for (i1=0; i1 < (jMaxv0) ; i1=i1+1)
		Ji=i1
		E1= ejv0[Ji]	
		q1=(exp(( -1* E1* k_h * k_c  )/(T* k_PlanckConstant )) )
		fac[(i1)]=(q1*(2*Ji+1))	// rotational degeneracy term 
		
		if (mod(Ji,2)==0)		//	even
			fac[i1]=ge*fac[i1]
		elseif (mod(Ji,2)==1)	//	odd
			fac[i1]=go*fac[i1]
		endif	
			
		su=fac[i1]+su	
	endfor

	// first vibrational state
	for (i1=0; i1<(jMaxv1) ; i1=i1+1)

		Ji = i1
		E1= ejv1[Ji]		//	v=1
		q1 = (exp(( -1* E1* k_h * k_c  )/(T* k_PlanckConstant )) )
		fac[(i1)] = (q1)*(2*(Ji)+1 )	// rotational degeneracy term 
		
		if (mod(Ji,2)==0)			//	even
			fac[jMaxv0+i1] = ge*fac[i1]
			su = fac[jMaxv0+i1] +su
		else						//	odd
			fac[jMaxv0+i1] = go*fac[i1]
			su = fac[jMaxv0+i1] +su
		endif
	endfor

	// ------------------------------------------------------------
	killwaves /Z fac
	return su

	end
//**********************************************************************************************************
//**********************************************************************************************************

// compute the sum of states for HD at a given temperature for upto Jf state.

function SumOfstates_HD (T)
	variable T
	// 	Argument
	//	T	: temperature, Kelvin

	variable  ge, go
	variable q1, su
	variable i1, Ji, E1


	// --- Nuclear spin statictics ------------------------------------------------------
	ge=1 // hydrigen deuteride
	go=1
	//---------------------------------------------------------------------------------------

	// waves containing the energy for each rotational state in the ground vib. state
	wave ejv0 = root:Params:RotnRamanIns:HDeV0
	wave ejv1 = root:Params:RotnRamanIns:HDeV1

	variable jMaxv0 = dimsize(ejv0, 0)
	variable jMaxv1 = dimsize(ejv1, 0)

	make /o /d /n=(jMaxv0 + jMaxv1) factor
	wave fac=factor 
	// ------------------------------------------------------------

	// ground vibrational state
	for (i1=0; i1 < (jMaxv0) ; i1=i1+1)
		Ji=i1
		E1= ejv0[Ji]	
		q1=(exp(( -1* E1* k_h * k_c  )/(T* k_PlanckConstant )) )
		fac[(i1)]=(q1*(2*Ji+1))	// rotational degeneracy term 
		
		if (mod(Ji,2)==0)		//	even
			fac[i1]=ge*fac[i1]
		elseif (mod(Ji,2)==1)	//	odd
			fac[i1]=go*fac[i1]
		endif	
			
		su=fac[i1]+su	
	endfor

	// first vibrational state
	for (i1=0; i1<(jMaxv1) ; i1=i1+1)

		Ji = i1
		E1= ejv1[Ji]		//	v=1
		q1 = (exp(( -1* E1* k_h * k_c  )/(T* k_PlanckConstant )) )
		fac[(i1)] = (q1)*(2*(Ji)+1 )	// rotational degeneracy term 
		
		if (mod(Ji,2)==0)			//	even
			fac[jMaxv0+i1] = ge*fac[i1]
			su = fac[jMaxv0+i1] +su
		else						//	odd
			fac[jMaxv0+i1] = go*fac[i1]
			su = fac[jMaxv0+i1] +su
		endif
	endfor

	// ------------------------------------------------------------
	killwaves /Z fac
	return su

	end

//**********************************************************************************************************
//**********************************************************************************************************

// compute the sum of states for D2 at a given temperature for upto Jf state.

function SumOfstates_D2 (T)
	variable T
	// 	Argument
	//	T	: temperature, Kelvin

	variable  ge,go
	variable q1,su
	variable i1,Ji, E1

	// --- Nuclear spin statictics ------------------------------------------------------
	ge=6 // molecular deuterium
	go=3
	//---------------------------------------------------------------------------------------

	// waves containing the energy for each rotational state in the ground vib. state
	wave ejv0 = root:Params:RotnRamanIns:D2eV0
	wave ejv1 = root:Params:RotnRamanIns:D2eV1

	variable jMaxv0 = dimsize(ejv0,0)
	variable jMaxv1 = dimsize(ejv1,0)

	make /o /d /n=(jMaxv0 + jMaxv1) factor
	wave fac=factor 
	// ------------------------------------------------------------

	// ground vibrational state
	for (i1=0; i1 < (jMaxv0) ; i1=i1+1)
		Ji=i1
		E1= ejv0[Ji]	
		q1=(exp(( -1* E1* k_h * k_c  )/(T* k_PlanckConstant )) )
		fac[(i1)]=(q1*(2*Ji+1))	// rotational degeneracy term 
		
		if (mod(Ji,2)==0)		//	even
			fac[i1]=ge*fac[i1]
		elseif (mod(Ji,2)==1)	//	odd
			fac[i1]=go*fac[i1]
		endif	
			
		su=fac[i1]+su	
	endfor

	// first vibrational state
	for (i1=0; i1<(jMaxv1) ; i1=i1+1)

		Ji = i1
		E1= ejv1[Ji]		//	v=1
		q1 = (exp(( -1* E1* k_h * k_c  )/(T* k_PlanckConstant )) )
		fac[(i1)] = (q1)*(2*(Ji)+1 )	// rotational degeneracy term 
		
		if (mod(Ji,2)==0)			//	even
			fac[jMaxv0+i1] = ge*fac[i1]
			su = fac[jMaxv0+i1] +su
		else						//	odd
			fac[jMaxv0+i1] = go*fac[i1]
			su = fac[jMaxv0+i1] +su
		endif
	endfor

	// ------------------------------------------------------------
	killwaves /Z fac
	return su

	end
//**********************************************************************************************************
//**********************************************************************************************************

// Procedure to calculate the intensity of rotational Raman bands of hydrogen at a given temperature.

function spectra_H2 ( T , nf_s, nf_as)   // for HYDROGEN
	variable T , nf_s, nf_as
	
	//	Arguments
	//	T		: temperature in Kelvin.
	//	nf_s	: the last J state for stokes sides.
	//	nf_as	: the last J state for anti-stokes sides.

	NVAR  omega=root:Params:RotnRamanIns:wavenumber
	variable omega_sc = omega / 1e4
	variable ge, go
	variable i1, Ji, sos 
	variable E, energy, popn, bj, posn, anisotropyME, factor

	// -------        set external waves here        -------------
	// matrix elements of gamma for H2.
	wave gw = root:Params:RotnRamanIns:'ME_gamma_532.199323_H2'
	wave ej  = root:Params:RotnRamanIns:H2eV0
	//----------------------------------------------------------------

	// --- Nuclear spin statictics -------------------------------
	ge=1	//	molecular hydrogen 
	go=3
	//----------------------------------------------------------------

	make /o /d /n = ((nf_s+nf_as), 4)	specH2
	make /o /d /n = (nf_s+nf_as)		posnH2		
	make /o /d /n = (nf_s+nf_as)		spectraH2

	//	sum of states for the temperature
	sos = SumOfstates_H2(T)

	// Stokes bands
	for (i1=0; i1<( nf_s +1) ; i1=i1+1)

		Ji=i1
		E = ej[Ji]
		energy = -1* E* k_h * k_c
		popn = (2*i1+1)* exp(energy/(T* k_PlanckConstant))
		bj=( 3*(Ji+1)*(Ji+2)) / (2*(2*Ji+1)*(2*Ji+3) )
		posn=(ej[(Ji+2)]- ej[(Ji)] )
		anisotropyME = gw[i1+1][2]

		if (mod(Ji,2)==0) 		//	even
			factor	=	 (popn * ge * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		elseif (mod(Ji,2)==1)	//	odd
			factor	=	 (popn * go * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		endif

			specH2 [((nf_as-1)+i1)][0]	=  i1
			specH2 [((nf_as-1)+i1)][1]	=  posn
			specH2 [((nf_as-1)+i1)][2]	=  factor		
			specH2 [((nf_as-1)+i1)][3]	=  omega - posn
			
			posnH2	[((nf_as-1)+i1)]		=  posn
			spectraH2 [((nf_as-1)+i1)]	=  factor
	endfor

	// ------------------------------------------------------------------------------
	// anti-Stokes bands 
	for (i1=(nf_as); i1>(1) ; i1=i1-1)   //

		Ji=i1
		E = ej[Ji]
		energy = -1* E* k_h * k_c	
		popn = (2*i1+1)* exp(energy/(T* k_PlanckConstant))
		bj=(3*Ji*(Ji-1))/(2*(2*Ji-1)*(2*Ji+1))
		posn = -1*(ej[(Ji)]- ej[(Ji-2)] )
		anisotropyME = gw[i1-1][2]
		
		if (mod(Ji,2)==0) 		//	even
			factor	=	 (popn * ge * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		elseif (mod(Ji,2)==1)	//	odd
			factor	=	 (popn * go * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		endif

			specH2 [(nf_as-i1)][0]	=  i1
			specH2 [(nf_as-i1)][1]	=  posn
			specH2 [(nf_as-i1)][2]	=  factor		
			specH2 [(nf_as-i1)][3]	=  omega - posn
			
			posnH2	[(nf_as-i1)]		=  posn
			spectraH2 [(nf_as-i1)]	=  factor
	endfor

	normalize_wave_to_one( spectraH2 )
	specH2 [ ][2]	= spectraH2[p]		// assign normalized intensity

	end			//	
	
//**********************************************************************************************************
//**********************************************************************************************************

// Procedure to calculate the intensity of rotational Raman spectrum of DEUTERIUM HYDRIDE (HD)

function spectra_HD(T, nf_s, nf_as)   // for hydrogen deuteride
	variable T, nf_s, nf_as

	//	Arguments
	//	T		: temperature in Kelvin.
	//	nf_s	: the last J state for stokes sides.
	//	nf_as	: the last J state for anti-stokes sides.

	NVAR  omega=root:Params:RotnRamanIns:wavenumber
	variable omega_sc = omega / 1e4
	variable ge, go
	variable i1, Ji, sos 
	variable E, energy, popn, bj, posn, anisotropyME, factor

	// -------        set external waves here        -------------
	// matrix elements of gamma for HD.
	wave gw= root:Params:RotnRamanIns:'ME_gamma_532.199323_HD'
	wave ej = root:Params:RotnRamanIns:HDeV0
	// ------------------------------------------------------------------

	// --- Nuclear spin statictics ----------------------------------------
	ge=1	// hydrogen deuteride
	go=1
	//-------------------------------------------------------------------------
	make /o /d /n = ((nf_s+nf_as), 4)	specHD
	make /o /d /n = (nf_s+nf_as)		posnHD		
	make /o /d /n = (nf_s+nf_as)		spectraHD

	//	sum of states for the temperature
	sos = SumOfstates_HD(T)

	// Stokes bands
	for (i1=0; i1<( nf_s +1) ; i1=i1+1)

		Ji=i1
		E = ej[Ji]
		energy = -1* E* k_h * k_c
		popn = (2*i1+1)* exp(energy/(T* k_PlanckConstant))
		bj=( 3*(Ji+1)*(Ji+2)) / (2*(2*Ji+1)*(2*Ji+3) )
		posn=(ej[(Ji+2)]- ej[(Ji)] )
		anisotropyME = gw[i1+1][2]

		if (mod(Ji,2)==0) //	even
			factor	=	 (popn * ge * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		elseif (mod(Ji,2)==1)	//	odd
			factor	=	 (popn * go * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		endif

			specHD [((nf_as-1)+i1)][0]	=  i1
			specHD [((nf_as-1)+i1)][1]	=  posn
			specHD [((nf_as-1)+i1)][2]	=  factor		
			specHD [((nf_as-1)+i1)][3]	=  omega - posn
			
			posnHD	[((nf_as-1)+i1)]		=  posn
			spectraHD [((nf_as-1)+i1)]	=  factor
	endfor

	// ------------------------------------------------------------------------------
	// anti-Stokes bands 
	for (i1=(nf_as); i1>(1) ; i1=i1-1)   //

		Ji=i1
		E = ej[Ji]
		energy = -1* E* k_h * k_c	
		popn = (2*i1+1)* exp(energy/(T* k_PlanckConstant))
		bj=(3*Ji*(Ji-1))/(2*(2*Ji-1)*(2*Ji+1))
		posn = -1*(ej[(Ji)]- ej[(Ji-2)] )
		anisotropyME = gw[i1-1][2]
		
		if (mod(Ji,2)==0) //	even
			factor	=	 (popn * ge * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		elseif (mod(Ji,2)==1)	//	odd
			factor	=	 (popn * go * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		endif

			specHD [(nf_as-i1)][0]	=  i1
			specHD [(nf_as-i1)][1]	=  posn
			specHD [(nf_as-i1)][2]	=  factor		
			specHD [(nf_as-i1)][3]	=  omega - posn
			
			posnHD	[(nf_as-i1)]		=  posn
			spectraHD [(nf_as-i1)]	=  factor
	endfor

	normalize_wave_to_one( spectraHD )
	specHD [ ][2]	= spectraHD[p]

	end		//	

//**********************************************************************************************************
//**********************************************************************************************************

// Procedure to calculate the intensity of rotational Raman spectrum of Deuterium

function spectra_D2(T , nf_s, nf_as)   // for molecular deuterium
	variable T, nf_s, nf_as

	//	Arguments
	//	T		: temperature in Kelvin.
	//	nf_s	: the last J state for stokes sides.
	//	nf_as	: the last J state for anti-stokes sides.

	NVAR  omega=root:Params:RotnRamanIns:wavenumber
	variable omega_sc = omega / 1e4
	variable ge, go
	variable i1, Ji, sos 
	variable E, energy, popn, bj, posn, anisotropyME, factor

	// -------        set external waves here        ----------------
	// matrix elements of gamma for D2.
	wave gw= root:Params:RotnRamanIns:'ME_gamma_532.199323_D2'
	wave ej =root:Params:RotnRamanIns:D2eV0
	// ------------------------------------------------------------------

	// --- Nuclear spin statictics ------------------------------------
	ge=6	// deuterium
	go=3
	//---------------------------------------------------------------------

	make /o /d /n = ((nf_s+nf_as), 4)	specD2
	make /o /d /n = (nf_s+nf_as)		posnD2		
	make /o /d /n = (nf_s+nf_as)		spectraD2

	//	sum of states for the temperature
	sos = SumOfstates_D2(T)

	// Stokes bands
	for (i1=0; i1<( nf_s +1) ; i1=i1+1)

		Ji=i1
		E = ej[Ji]
		energy = -1* E* k_h * k_c
		popn = (2*i1+1)* exp(energy/(T* k_PlanckConstant))
		bj=( 3*(Ji+1)*(Ji+2)) / (2*(2*Ji+1)*(2*Ji+3) )
		posn=(ej[(Ji+2)]- ej[(Ji)] )
		anisotropyME = gw[i1+1][2]

		if (mod(Ji,2)==0) //	even
			factor	=	 (popn * ge * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		elseif (mod(Ji,2)==1)	//	odd
			factor	=	 (popn * go * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		endif

			specD2 [((nf_as-1)+i1)][0]	=  i1
			specD2 [((nf_as-1)+i1)][1]	=  posn
			specD2 [((nf_as-1)+i1)][2]	=  factor		
			specD2 [((nf_as-1)+i1)][3]	=  omega - posn
			
			posnD2	[((nf_as-1)+i1)]		=  posn
			spectraD2 [((nf_as-1)+i1)]	=  factor
	endfor

	// ------------------------------------------------------------------------------
	// anti-Stokes bands 
	for (i1=(nf_as); i1>(1) ; i1=i1-1)   //

		Ji=i1
		E = ej[Ji]
		energy = -1* E* k_h * k_c	
		popn = (2*i1+1)* exp(energy/(T* k_PlanckConstant))
		bj=(3*Ji*(Ji-1))/(2*(2*Ji-1)*(2*Ji+1))
		posn = -1*(ej[(Ji)]- ej[(Ji-2)] )
		anisotropyME = gw[i1-1][2]
		
		if (mod(Ji,2)==0) //	even
			factor	=	 (popn * ge * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		elseif (mod(Ji,2)==1)	//	odd
			factor	=	 (popn * go * bj* omega_sc * (( omega_sc - posn/1e4)^3)*   anisotropyME^2) / sos
		endif

			specD2 [(nf_as-i1)][0]	=  i1
			specD2 [(nf_as-i1)][1]	=  posn
			specD2 [(nf_as-i1)][2]	=  factor		
			specD2 [(nf_as-i1)][3]	=  omega - posn
			
			posnD2	[(nf_as-i1)]		=  posn
			spectraD2 [(nf_as-i1)]	=  factor
	endfor


	normalize_wave_to_one( spectraD2 )
	specD2 [ ][2]	= spectraD2[p]		// assign normalized intensity

	end			//	

//**********************************************************************************************************
//**********************************************************************************************************
//	Function to normalize a 1D wave using the max value in that wave

function normalize_wave_to_one (wave_name)
		wave wave_name
		
		variable max_val = WaveMax( wave_name)
		wave_name = wave_name / max_val 

end

//**********************************************************************************************************
//**********************************************************************************************************