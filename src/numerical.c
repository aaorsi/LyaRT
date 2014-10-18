/* Different numerical functions used in the code */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

double int_tabulated_trap_rule(double *xarr, double *yarr, long n)
{
// Simple tabulated integral (not even a trapezoidal rule). 
// xarr and yarr must have the same # of elements

    double sum;
    double h,maxy,miny;
	long i;

    sum = 0.;
	for (i=0;i<n-1;i++)
	{
       h = xarr[i+1] - xarr[i];
	   maxy = MAX(yarr[i+1],yarr[i]);
	   miny = MIN(yarr[i+1],yarr[i]);


	   sum += h * (miny + 0.5*(maxy - miny));

//     sum += h * 0.5 * (yarr[i+1] + yarr[i]);
    }
	
	return(sum);
}

