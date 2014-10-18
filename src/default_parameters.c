/* Here all default parameters are defined.
 Changes from these default values must be introduced
 from a parameter file, read by read_parameters.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

void default_parameters(void)
{ 

	strcpy(GeomName,"HomSlab");
	strcpy(OutLong,"../out/long/HomSlab.ip");
	strcpy(OutShort,"../out/short/shortfile.ip");
	strcpy(GHIFile,"../data/gamma");
	strcpy(UVzFile,"../data/HI_table/redshift.txt");
	strcpy(UVeFile,"../data/HI_table/logenergy_ryd.txt");
	strcpy(UVspecFile,"../data/HI_table/logflux.txt");
	strcpy(Set_Tolerance,"no");
	redshift		= 0;
	dX                      = 0;
	IncShield		= 0;
	IncUV			= 0;
	nout_max 		= 1e6;
	np_max			= 1e6;
	np_min 			= 100;
	Static			= 1;
	NPhotons 		= 1;
	NCells			= 1;
	NCellX   		= 1;
	NCellY   		= 1;
	NCellZ			= 1;
	XPeriodic		= 0;
	YPeriodic		= 0;
	ZPeriodic		= 0;
	VarXCrit                =0;
	xSize			= 3e18;
	ySize			= 3e18;
	zSize			= 3e18;
	RSphere			= 3e18;
	Temp			= 0.0;
	b				= 0.0;
	vmax			= 0.0;
	Vrot			= 0.0;
	x_mean			= 0.0;
	ColDens                 = 0.0;
	mean_nH			= 10; 
	mean_nH_static	= 0.0;
	xcritval		= 0.0;
	Ndot			= 1e8;
	R_inner			= 0.0;
	R_Peak                   =  1.5e18;
	R_Static		= 1.0;
	StartAtCenter 	= 1.0;	// If equal to 0 then photons are produced in random positions.
	GetVpRej		= 1.0;	// Use Rejection method to compute vp

	DefinePosAndDir = 0;
	Ix0				= 0.;
	Iy0				= 0.;
	Iz0				= 0.;
	Ith0			= 0.;
	Iph0			= 0.;


//	Velocity profile parameters
	alpha_vprof		= 1.4;
	sigma_vprof		= 150.0;
	
	
    strcpy(tabledir,"../data/tables/");
	strcpy(sxfile,"Hx_T_");
	strcpy(OutMode,"Short");	//It could be 'Short' or 'Long'
	strcpy(IncDust,"Yes");		// Include dust, 'Yes' or 'No'
	

#ifdef TEST_RII
	x_test 		= 0.0;
	nscatmax	= 100000;
	sprintf(OutLong,"../out/long/test_RII/xi_%.1f",x_test);
#endif	

}
