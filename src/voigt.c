#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

double voigt(double * HList, double xp)
{
	long i;
	double x_0,x_1,H0,H1,H_x,axp;
	double q,z,sqrz,pi_z;
	double sqrx;

	axp = fabs(xp);

	if (axp > MAXX)
		H_x = a_par/(sqrtPi*SQR(xp));		// Power law approximation when x > MAXX, numerical calculation from a lookup table otherwise
	else
	{

#ifdef HAPPROX

		sqrx = SQR(xp);
		z 	 = (sqrx - hc1)/(sqrx + hc2);
		sqrz = SQR(z);
		pi_z = hc3 * (sqrz*sqrz) - hc4 * sqrz*z + hc5*sqrz + hc6 * z;
		q = (z > 0 ) ? (1 + 21/sqrx) * (a_par*pi_z)/(Pi*(sqrx+1)) : 0;
		H_x = q * sqrtPi + exp(-sqrx);

#else

		i = 0;
		while (axp >= HList[i] && i < nHList)
		{
			x_0 = HList[i];
			H0 = HList[i + NTAB2];
			i++;
		}
		x_1 = HList[i];
		H1 = HList[i+NTAB2];
		H_x = (H1 - H0)/(x_1 - x_0) * (axp - x_1) + H1;
		if (i == nHList)
			H_x = H0;

#endif

	}
	
	return(H_x);
}
	
