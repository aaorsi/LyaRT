/* To output data containing the position and frequency of the photon.
This routine will create one file for each photon, so is not recommended
when the number of photons is large */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"


void record_data_long(long nscat, long ip, float x, float rxf, float ryf, float rzf, float radius, float upar,\
float uper1, float uper2, float H_x,float xp, int inter)
{

	FILE *fout;
	char OutFile[NCHARMAX];
	char sp[10];	
	char itype[20];
	
	sprintf(sp,"%ld",ip);
	strcpy(OutFile,OutLong);
	strcat(OutFile,sp);

	switch(inter)
	{	
		case 1:
		strcpy(itype,"Hydrogen");
		break;
		case 2:
		strcpy(itype,"Dust");
		break;
		case 3:
		strcpy(itype,"Absorbed_by_dust");
		break;
		case 4:
		strcpy(itype,"Photon_escaped");	
		break;
	}


	if (fout = fopen(OutFile,"a"))
	{
		if (nscat == 0)
		{
			fprintf(fout,"# Photon %ld\n",ip);
			fprintf(fout,"# NScatt	x[freq]		H(x,a)	(rx)	(ry)	(rz)	radius	upar	uper1	uper2 xp Interaction \n");
		}
		
		fprintf(fout," %ld\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %s\n",nscat,x,H_x,(rxf),(ryf),(rzf),\
		radius,upar,uper1,uper2,xp,itype);	
		fclose(fout);
	}
	else
	{
		printf("File %s could not be opened. fout = %d\n Exiting\n\n",OutFile,fout);
		fclose(fout);
		exit(0);
	}
}
		
