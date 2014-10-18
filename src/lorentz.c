// This function will compute the new direction of the photon (theta, phi) 
// expressed in the other frame (observer or fluid).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"


void lorentz(int dir)
{
	// dir is +1 when the transformation goes from the moving frame back to the lab frame, and -1 if the opposite

	double vnorm,vbx,vby,vbz,gf,npar,alpha_f,nper_bx,nper_by,nper_bz,npern,\
	npar_f,nper_f,norm_f;	

	vnorm = sqrt(SQR(vbulk_x) + SQR(vbulk_y) + SQR(vbulk_z));

	if(vbulk_x == 0. && vbulk_y == 0. && vbulk_z == 0.)
	{
//		printf("zero velocity\n");
//		dprinti(nscat);
		return;
	}
	
	vbx = vbulk_x/vnorm;
	vby = vbulk_y/vnorm;
	vbz = vbulk_z/vnorm;
			
	gf = 1./sqrt(1 - SQR(H*radius/c));
	
	npar = ni*vbx + nj*vby + nk*vbz;
	alpha_f = acos(npar);

	if (alpha_f == 0)
	{
//		printf("alpha_f = %f\n",alpha_f);
		return;
	}
	
	nper_bx = (ni - cos(alpha_f)*vbx)/sin(alpha_f);
	nper_by = (nj - cos(alpha_f)*vby)/sin(alpha_f);
	nper_bz = (nk - cos(alpha_f)*vbz)/sin(alpha_f);
	
	npern = sqrt( SQR(nper_bx) + SQR(nper_by) + SQR(nper_bz));
	
	nper_bx = nper_bx/npern;
	nper_by = nper_by/npern;
	nper_bz = nper_bz/npern;
		
	npar_f = (c * cos(alpha_f) + dir*H*radius)/(1 + dir*cos(alpha_f)*H*radius/c);
	nper_f = (c * sin(alpha_f))/(gf*(1 + dir*cos(alpha_f)*H*radius/c));
							
	ni = (1./c) * ( npar_f * vbx + nper_f * nper_bx);
	nj = (1./c) * ( npar_f * vby + nper_f * nper_by);
	nk = (1./c) * ( npar_f * vbz + nper_f * nper_bz);
	
	norm_f = sqrt(SQR(ni) + SQR(nj) + SQR(nk));
	
	ni = ni/norm_f;
	nj = nj/norm_f;
	nk = nk/norm_f;
	
	th0 = acos(nk);
	ph0 = atan2(nj,ni);
	ph0 = (ph0 < 0.) ? (ph0 + 2*Pi) : ph0;

	return;	
}
