/* To output data containing the position and frequency of the photon.
This routine will creatie one single file for all photons.
For a more detailed output, see RECORD_DATA_LONG*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"


//void record_data_short(long nscat, long ip, float x, float rxf, float ryf, float rzf, float radius, \
//float H_x,int inter,float op_time, long flag_zero)

void record_data_short(long nscat,long ip,float x,float xp0,float ph0, float th0, float radius, int inter,float op_time,long flag_zero)
{

	FILE *fout;
	char OutFile[NCHARMAX];
	char sp[10];	
	char itype[20];
	
	sprintf(sp,"%ld",ip);
	strcpy(OutFile,OutShort);

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
//		printf("Writing in file %s\n",OutFile);
//		dprints(OutShort);
		if (ip == 0)
		{
			fprintf(fout, \
//		"# NScatt	x[freq]		H(x,a)		log(rx)		log(ry)		log(rz)	log(radius)		Time[min]	Interaction \n");
		"# NScatt	x[freq]		Phi[rad]	Theta[rad]	log10(R)	Time[min]	NBackScatt x0[freq]	Interacion \n");
		}
		
//		fprintf(fout," %ld\t %f\t %f\t %f\t %f\t %ld\t %s\n",nscat,x,H_x,log10(radius),op_time,flag_zero,itype);	
//		fprintf(fout," %ld\t %f\t %f\t %f\t %f\t %f\t %ld\t %s\n",nscat,x,ph0,th0,log10(radius),op_time,flag_zero,itype);
		fprintf(fout," %ld\t %f\t %s\t %f\n",nscat,x,itype,xp0);
		fclose(fout);
	}
	else
	{
		printf("File %s could not be opened. fout = %d\n Exiting\n\n",OutFile,fout);
		fclose(fout);
		exit(0);
	}
}
		
