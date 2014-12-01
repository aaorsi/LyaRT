#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

void init_photon(photon *P, int ip, double xi0, double xi1, double xi2)
{
	double xi,xp0;
	int upcone,i;

	
	if (strcmp(GeomName ,"Biconical_Wind")==0)
	{
			
		upcone = (xi0 < 0.5) ? 1 : -1;	// This random number chooses the hemisphere
//			th0 = acos(upcone*(xi1*(cos(app_angle)-1)+1)); // theta is constrained to be drawn from 0,app_angle, or pi,pi-app_angle
		th0 = (upcone == 1) ? acos(xi1*(1-cos(app_angle))+cos(app_angle) ) : \
								  Pi - acos(xi1*(1-cos(app_angle))+ cos(app_angle));
		P[ip].th = th0;
	}
	else 
	{
		th0 = acos(2*xi1 - 1);
		P[ip].th = th0;
	}

   	ph0 = 2*Pi*xi2;
	P[ip].ph = ph0;

	if (StartAtCenter== 1)
	{
//			idc = (int) (NCells-1)/2;
//			if (gtype == 1)

		idc = 0;
		idc_old = 0;	
		rx0 = 0.;			
		ry0 = 0.;			
		rz0 = 0.;		

		P[ip].x = 0.;
		P[ip].y = 0.;
		P[ip].z = 0.;
		P[ip].idc = 0;

	}
	else
	{
		rx0 = Ix0;
		ry0 = Iy0;
		rz0 = Iz0;
	
		th0	= Ith0;
		ph0 = Iph0;
	
		r0  = sqrt(SQR(rx0) + SQR(ry0) + SQR(rz0));
				
		idc = 0;
		idc_old = 0;
		P[ip].x = Ix0;
		P[ip].y = Iy0;
		P[ip].z = Iz0;

		while (r0 <= CellArr[i].r + dr/2.)
		{
			idc = i;
			i++;
		}
		P[ip].idc = idc;

	}


//		b = sqrt(v_th^2 + v_turb^2). If defined, then vth = b

			

//      printf("Initial angles: %f %f\n",th0 * 180/Pi,ph0*180/Pi);
//      printf("Initial position: (%f %f %f)\n",(rx0),(ry0),(rz0));

// 		Initial frequency

		xp = 0; // initial frequency in the rest frame of the fluid 
		xp0=xp;		
// Flat Continuum
#ifdef XCONTINUUM
		//dX=96.0;
		xi1 = gsl_rng_uniform (r);
		
		xp = (xi1-0.5)*dX;
		xp0=xp; 
#endif		
		
#ifdef XINDIST
		
		xi1 = gsl_rng_uniform (r);
		xi2 = gsl_rng_uniform (r);
		
		xp	= x_mean + Vrot * sqrt(-2*log(xi1))*cos(2*Pi*xi2);
		xp0=xp; 
		dprintl(ip);
		dprintd(xp);
		printf("------");
		//#else
		//xp = 0; // initial frequency in the rest frame of the fluid
		//xp0=xp; 
		//printf("0 xp: %f xp0=%f \n",xp,xp0);
#endif 
		x = xp;
		xp0=xp; 
		P[ip].xp = xp;
		//printf("1 xp: %f xp0=%f \n",xp,xp0);
		xabs = fabs(x);

		inter 	= 1;

		if (StartAtCenter == 1)
		 {

			r0 				= 1;
			radius 			= 0;

			vbulk_x			= 0;
			vbulk_y			= 0;
			vbulk_z			= 0;
			vbulk0_x		= 0;
			vbulk0_y		= 0;
			vbulk0_z		= 0;
			
			px0 			= 0.;
			py0 			= 0.;
			pz0 			= 0.;

		}
		else
		{

			px0 = rx0/r0;
			py0 = ry0/r0;
			pz0 = rz0/r0;

			radius = r0;

			vbulk_x		= CellArr[idc].vbulk[0] * px0;
			vbulk_y		= CellArr[idc].vbulk[0] * py0;
			vbulk_z		= CellArr[idc].vbulk[0] * pz0;
			vbulk0_x	= vbulk_x;
			vbulk0_y	= vbulk_y;
			vbulk0_z	= vbulk_z;
			
		}	
}
