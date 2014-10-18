/* This corrects the final position of a photon escaping from the 
   Biconical outflow to lie exactly at the edge of the cone. */

#include <math.h>
#include "allvars.h"
#include "proto.h"
void find_edgecone(float px0,float py0, float pz0, float r0, float th0, float ph0, double pos[])
{
	double thp, php, re, alpha, phe, app;

	thp = acos(pz0/r0);		// Polar angle of photon position prior to correction
	php = atan2(py0,px0);	// azimuthal angle of photon position

	app = (pz0 > 0) ? app_angle : Pi - app_angle;

	// re is the distance to the edge where the photon actually escaped from the cone
	re = r0 * sin(Pi/2.-thp) / (cos(app)- (sin(thp-app)*cos(Pi-th0)/sin(th0-thp))); 
	// r0 is the distance traveled from re to r0
	r0 = re * sin(thp-app)/sin(th0-thp);
	
	alpha = asin( (r0/re) * sin(Pi-th0)/sin(app) * sin(ph0-php));
	// phe os the azimuthal angle of the photon at the edge. (the polar one is app)
	phe = php-alpha;	   

	// updating positions:
	pos[0] = re * sin(app)*cos(phe);
	pos[1] = re * sin(app)*sin(phe);
	pos[2] = re * cos(app);

	
}

