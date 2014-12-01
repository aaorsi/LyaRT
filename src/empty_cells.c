#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

void empty_cells(photon *P, int ip)
{
//  The following lines define the behaviour of the photons when crossing (practically) empty cells
// flag_zero counts the number of times the photon crosses an empty cell. 

	if (CellArr[idc].nH < 1e-10 && idc != -1)
	{ 		  
		
		px0 = P[ip].x/r0;
		py0 = P[ip].y/r0;
		pz0 = P[ip].z/r0;

		// This corresponds to the treatment of backscatterings.
		if (idc == 0 && idc_old == 1)
		{
//			thc = acos(P[ip].z/sqrt(SQR(P[ip].x) + SQR(P[ip].y) + SQR(P[ip].z)));
			garg = -px0*P[ip].ni  -py0*P[ip].nj  -pz0*P[ip].nk;
//			b_ = th0 - thc - Pi;
//			s_ = 2 * CellArr[0].r * cos(b_);
			s_ = 2 * CellArr[1].r * garg;

			P[ip].x += s_*P[ip].ni;
			P[ip].y += s_*P[ip].nj;
			P[ip].z += s_*P[ip].nk;

			r0   = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z) );
//			idc = (idc_old > idc && idc != 0) ? CellArr[idc].neig[0] : CellArr[idc].neig[1];
			idc = 1;

		}
		// The following in general for almost empty cells.
		else if (flag_zero >0 && idc != 0)	
		{

			px0 = P[ip].x/r0;
			py0 = P[ip].y/r0;
			pz0 = P[ip].z/r0;
			r0   = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z) );
			garg = -px0*ni - py0*nj - pz0*nk;
						
			g_   = (garg >= 1.0) ? 0. : (garg <=-1.0) ? Pi : acos(garg);
			b_   = Pi - 2*g_;
			s_   = (g_ == 0 ) ? 2 * CellArr[0].r : r0 * sin(b_)/sin(g_);

			if ( flag_zero > 1e20)	
			{ 
				printf("+++++++++++++++++++++++++++++\n");	
				printf("WARNING: flag_zero>NCells\n");
				dprintd(xp);	
				dprintd(t_0);
				dprintd(radius);
				dprintd(r0);
				dprintd(rE);
				dprintd(px0);
				dprintd(py0);
				dprintd(pz0);
				dprintd(ni);
				dprintd(nj);
				dprintd(nk);
				dprintd(g_);
				dprintd(b_);
				dprintd(s_);
				dprintl(idc);
				dprintl(flag_zero);
				dprintl(nscat);
				dprintd(dr);
				printf("EXITING...\n");
				exit(0);
		   	}  

			P[ip].x += s_*ni;
			P[ip].y += s_*nj;
			P[ip].z += s_*nk;

			idc = (idc_old > idc && idc != 0) ? CellArr[idc].neig[0] : CellArr[idc].neig[1];
							
		}
		else
		{	
		//	Empty cell at the start	
			idc = (idc_old > idc && idc != 0) ? CellArr[idc].neig[0] : CellArr[idc].neig[1];
			P[ip].x = CellArr[idc].r*P[ip].ni;
			P[ip].y = CellArr[idc].r*P[ip].nj;
			P[ip].z = CellArr[idc].r*P[ip].nk;
			r0 = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z) );
			
		}

		px0 = P[ip].x/r0;
		py0 = P[ip].y/r0;
		pz0 = P[ip].z/r0;

		vbulk0_x = vbulk_x;
		vbulk0_y = vbulk_y;
		vbulk0_z = vbulk_z;

		vbulk_x = px0 * CellArr[idc].vbulk[0];
		vbulk_y = py0 * CellArr[idc].vbulk[0];
		vbulk_z = pz0 * CellArr[idc].vbulk[0];
		
		P[ip].xp += (P[ip].ni * (vbulk0_x-vbulk_x) + P[ip].nj*(vbulk0_y-vbulk_y) + P[ip].nk*(vbulk0_z-vbulk_z) )/vth;

		flag_zero++;

					
	}
}
