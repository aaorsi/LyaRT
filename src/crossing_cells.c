#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

double crossing_cells(photon *P, int gtype, int ip)
{

	double comp,rlow,rhigh;
	double b1,b2,r0_old,r0_alt,s_old;
	double th_c,b3,b4,s_excess,s_edge;


	rx0 = P[ip].x;
	ry0 = P[ip].y;
	rz0 = P[ip].z;
	r0 = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z));
	

	ni = P[ip].ni;
	nj = P[ip].nj;
	nk = P[ip].nk;

	switch(gtype)
	{
		case 1:

//			First: easy case, identify when photon moves along the radial direction
			
			px0 = P[ip].x/r0;
			py0 = P[ip].y/r0;
			pz0 = P[ip].z/r0;

			comp = (r0 - CellArr[idc].r)/CellArr[idc].r;
	
			rlow = ( comp > -0.01) ? CellArr[idc].r : CellArr[idc-1].r;
			rhigh = rlow + dr;

			garg = -px0*ni  -py0*nj  -pz0*nk;
			g_   = (garg >= 1.0) ? 0. : ((garg <= -1.0) ? Pi : acos(garg));
			if (idc > 0)
			{
				aarg = rlow/r0;
				acrit = (aarg >= 1.0) ? Pi/2. : ( (aarg <= -1.00) ? -Pi/2. : asin(aarg) );
				comp = (fabs(garg) > 1) ? fabs(garg) - 1.0 : 1.0 - fabs(garg);
				rE   = (comp < 1e-4) ? ((g_ >= acrit) ? rhigh : rlow) : -1.0;
				s_   = (garg > 0) ? (r0 - rlow) : (rhigh - r0);		// s_, the distance traveled, is recomputed if rE = -1 here.
			}
			else
			{
//	if the photon is in idc == 0				
				rE = rhigh;
				b2 = asin(r0 * sin(g_)/ rE);
				b1 = Pi - (b2 + g_);
				s_ = rE * sin(b1)/sin(g_);
				acrit = 0.0;
			}		

//			When photon is not moving radially, identified by rE = -1:
			if (idc > 0 && rE == -1.0)
			{

				rE = (g_ >= acrit ) ? rhigh : rlow;

				b2 = asin(r0 * sin(g_)/ rE);
				b2 = (g_ >= acrit) ? b2 : Pi - b2;
				b1 = Pi - (b2 + g_);
				s_ = rE * sin(b1)/sin(g_);

			}

//	When the distance computed, s_, is smaller than the one corresponding to t_0, s, then the photon
//  crosses to another cell:

			if (fabs(s_) < s)
			{
				idc_old = idc;
				idc = (g_ >= acrit ) ? CellArr[idc].neig[1] : CellArr[idc].neig[0];

				t_used = s_sum * s_;
				t_0 -= t_used;
			}
			else 
			{		
				r0_alt = sqrt( SQR(P[ip].x+s_*P[ip].ni) + SQR(P[ip].y+s_*P[ip].nj) + SQR(P[ip].z+s_*P[ip].nk));
				s_old = s_;
				s_ = s;
				t_0 = -1;
				idc_old = idc;
			
			}

			if (t_0 <= 0)
				t_0 = -1;


			P[ip].x +=  s_*P[ip].ni;
      P[ip].y +=  s_*P[ip].nj;
	    P[ip].z +=  s_*P[ip].nk;	

			r0_old = r0;
			r0 = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z));

			if (idc == -1)
			{
				t_0 = -1;
				EscCond = 1;
				break;
			}

			if ((r0 > 1.01*CellArr[idc+1].r || r0 < 0.99*CellArr[idc].r) && idc != NCells-1)
			{	
				printf("photon position does not corresponds with cell identifier. Check!\n");
				dprintd(r0);
				dprintd(CellArr[idc].r);
				dprintd(CellArr[idc+1].r);
				dprinti(idc);
				dprinti(nscat);
				printf("stop\n");
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
		
//				if (idc != idc_old)
					P[ip].xp += (P[ip].ni * (vbulk0_x-vbulk_x) + P[ip].nj*(vbulk0_y-vbulk_y) + \
						P[ip].nk*(vbulk0_z-vbulk_z))/vth;

			if isnan(P[ip].xp)
				printf("crossing_cells.c: xp is nan\n");

		
			if (strcmp(GeomName,"Biconical_Wind")==0 && idc != 0)
			{		
				th_c = acos(pz0);
				EscCond = (pz0 > 0) ? ( (th_c > app_angle) ? 1 : 0) : \
						  ((th_c < Pi - app_angle) ? 1 : 0);

				if (EscCond == 1)
				{

/*					ObscuredCond = (pz0>0) ? ( (th0 > Pi/2.) ? 1 : 0) :\
								((th0 < Pi/2.) ? 1 : 0);
					
					if (ObscuredCond == 1)
					{
//						printf("Photon %d absorbed through perpendicular plane coming from %s side of cone:\ntheta=%f phi=%f r/RSphere=%f app=%f, nscat=%ld\n\n",\
						ip,(pz0 > 0) ? "+" : "-",th0*180./Pi,ph0*180./Pi,r0/RSphere,app_angle*180./Pi,nscat);
						idc = -1;
						goto absorbed;
					}
*/

					t_0 = -1;
					idc = -1;
/*
					b3 = th_c - app_angle;
					b4 = Pi - (b3 + b2);
					s_excess = r0 * sin(b3)/sin(b4);
					s_edge   = s_ - s_excess;

					P[ip].x = rx0 + s_edge*P[ip].ni;
					P[ip].y = ry0 + s_edge*P[ip].nj;
					P[ip].z = rz0 + s_edge*P[ip].nk;

					r0 = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z));
//								printf("corrected position: %f,%f,%f\n",rx0,ry0,rz0); 
*/

				}
			}
				
/*			if (r0 < CellArr[idc].r || r0 > CellArr[idc].r+dr)
			{
				printf("ERROR: Something wrong happened, r0 is outside the limits of cell %d\n",idc);
				dprintd(CellArr[idc].r);		
				dprintd(CellArr[idc].r+dr);
			}
//							dprintd(r0);
							r0 = (r0 < CellArr[idc].r) ? CellArr[idc].r : CellArr[idc].r+dr;	
							rx0 = px0 * r0;
							ry0 = py0 * r0;
							rz0 = pz0 * r0;
			}
*/
			break;
						
		case 0:	

			t_0 = -1;
			vbulk_x = 0;
			vbulk_y = 0;
			vbulk_z = 0;

			P[ip].x +=  s*P[ip].ni;
        	P[ip].y +=  s*P[ip].nj;
	    	P[ip].z +=  s*P[ip].nk;	

			r0_old = r0;
			r0 = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z));

			s_ = s;

			if (fabs(P[ip].z) >= zSize/2.)
			{
				EscCond =  1;	
				idc		= -1;	
				break;
			}
						
			break;

		case 2:
			t_0 = -1;
			vbulk_x = 0;
			vbulk_y = 0;
			vbulk_z = 0;

			P[ip].x +=  s*P[ip].ni;
      P[ip].y +=  s*P[ip].nj;
	    P[ip].z +=  s*P[ip].nk;	

			r0_old = r0;
			r0 = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z));

			s_ = s;

			if (r0 >= RSphere)
			{
				EscCond 	=  1;	
				idc		= -1;	
				break;
			}
			break;
	}

	return(s_);

}

