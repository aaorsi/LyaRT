/* Read a look-up table with tabulated values of
 sigma_x as a function of x and T */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"



void get_H(char *SigmaFile, double *HList)
{
	FILE *fd;	
	int i;
	char buf[NCHARMAX];
	double Hi;
	double H1;


	if(fd = fopen(SigmaFile,"r"))
	{
		fgets(buf,50,fd);
		i=0;

 	    while (!feof(fd))
        {
            fscanf(fd,"%lf %lf",&HList[i],&HList[i + NTAB2]);
			i++;
		}
    }
    else
    {
        printf("Cannot open file %s\n Exiting..",SigmaFile);
        exit(0);
    }
	fclose(fd);
	nHList = i-1;
}	
			
	
