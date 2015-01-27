#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

void dust_interaction(photon *P, int ip)
{
	double xilow,xihigh,mulow,muhigh,hglow,hghigh;
	long i;

	interd = (xi1 <= Albedo) ? 1 : 0;
	if (interd == 1)	//scattering
	{

		i = 0;
		while (xi2 >= HGList[i] && i < nHG)
		{
			xilow = HGList[i];
			hglow = HGList[i + NTAB2];
			i++;
		}
		xihigh = HGList[i];
		hghigh = HGList[i+NTAB2];
		th0 = (hghigh-hglow)/(xihigh-xilow) * (xi2- xihigh) + hghigh;
		if (i == nHG)
			th0 = hglow;
   		ph0 = 2*Pi*xi3;

		P[ip].th = th0;
		P[ip].ph = ph0;
		end_syg = 0;
	}
	else 
	{
		absorbed:
		P[ip].xp += (ni*vbulk_x + nj*vbulk_y + nk*vbulk_z)/vth;
//			printf("Photon has been absorbed by dust after %ld scatters.\n",nscat);
		inter = 3;
//		(void) time(&t2);
//		op_time = (float) (t2-t1)/60.;
//				printf("Total calculation took %f minutes.\n",op_time);
		if (strcmp(OutMode,"Long")==0)
		{
			record_data_long(nscat,ip,x,rxf,ryf,rzf,radius,upar,uper1,uper2,H_x,\
			xp,inter);
			dprintl(nscat);
			printf("\n");
		}
								
//							record_data_short(nscat,ip,x,rxf,ryf,rzf,radius,H_x,inter,flag_zero,op_time);
//							record_data_short(nscat,ip,x,ph0,th0,radius,inter,op_time,flag_zero);

		XArr[ip] = P[ip].xp;
		X0Arr[ip]=xp0;
		InterArr[ip] = inter;
		NscatArr[ip] = nscat;
#ifdef GETPOSDIR
		PosArr[3*ip] 	= rxf;
		PosArr[3*ip+1]	= ryf;
		PosArr[3*ip+2]	= rzf;
		AngleArr[2*ip]	= th0;
		AngleArr[2*ip+1]= ph0;	

//						record_data_pos(nscat,ip,x,ph0,th0,rxf,ryf,rzf,inter,op_time,flag_zero);
//#else
//						record_data_short(nscat,ip,x,ph0,th0,radius,inter,op_time,flag_zero);
#endif
/*
#ifdef GETPOSDIR
						record_data_pos(nscat,ip,x,ph0,th0,rxf,ryf,rzf,inter,op_time,flag_zero);
#else
						record_data_short(nscat,ip,x,ph0,th0,radius,inter,op_time,flag_zero);
#endif
*/	
		idc = -1;
		P[ip].idc = idc;
		end_syg = 1;
	}
}						
