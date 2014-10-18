/* This code computes the velocity component of the atom parallel to the 
photon's direction by using the rejection method (see Zheng & Miralda-Escude 2002
and Laursen et al. 2009 for references.)

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

float vp_rejection(double x,double a,double xi1,double xi2,double xi3)
{
	double u0,theta_0,theta,p,u,RF,xminus,xcw;
	double loga, Piotwo,expsqu0;
	xminus = 1;

	loga 	= log10(a);
	Piotwo	= Pi/2.; 

	xcw = 1.59 - 0.60*loga - 0.03*SQR(loga);

	if (x < 0.)
	{
		xminus = x;
		x = fabs(x);
	}
	
	if ( x >=0 && x < 0.2)
		u0 = 0;
	else if (x >= 0.2 && x < xcw)
		u0 = x - 0.01*pow(a,1./6)*exp(1.2*x);
	else if (x >= xcw)
		u0 = 3.5;	// In Laursen's paper, u0 = 4.5 here. 

	expsqu0 = exp(-SQR(u0));

	theta_0 = atan((u0 - x)/a);
	p = (theta_0 + Piotwo) * 1./((1. - expsqu0)*theta_0 + (1. + expsqu0)*Piotwo ) ;
	
	if (xi1 <= p)
		theta = xi2 * (theta_0+Piotwo) - Piotwo;
	else
		theta = xi2 * (Piotwo-theta_0) + theta_0;


	u = a * tan(theta) + x;


	if ( u <= u0)
		RF = exp(-SQR(u));
	else
		RF = exp(-SQR(u))/expsqu0;


	if (xi3 > RF)
	{
		u = -999.0;
		return(u);
	}
	else
	{
	if (xminus <0)
		u = -u;
	return(u);
	}
}
	
