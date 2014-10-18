/* This function receives a random number between 0<xi<1
and returns a value for polar angle off the initial photon's direction
when scattered by an hydrogen atom. */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

void get_dipolar(double *dipolar_list)
{
	FILE *fd;
	long i;
	char dipfile[NCHARMAX];
	char buf[NCHARMAX];
	strcpy(dipfile,tabledir);
	strcat(dipfile,"dipolar/dipolar");
	if (fd = fopen(dipfile,"r"))
	{
		dprints(dipfile);
		fgets(buf,20,fd);
		i = 0;
		while (!feof(fd))
		{
			fscanf(fd,"%lf %lf",&dipolar_list[i],&dipolar_list[i + NTAB1]);

			i++;
/*
			if (xi_list <= ran)
			{		
				xi0 = xi_list;	
				mu0 = mu_list;
			}
			else
			{	
				xi1 = xi_list;
				mu1 = mu_list;
				break;
			}					
*/
		}
	}
	else
	{
		printf("Cannot open file %s\n Exiting..",dipfile);
		exit(0);
	}

	fclose(fd);
	ndip = i-1;	
//	mu = (mu1 - mu0)/(xi1 - xi0) * (ran - xi1) + mu1;

//	return mu;
}



		

