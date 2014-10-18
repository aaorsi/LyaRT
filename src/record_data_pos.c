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

//void record_data_pos(long nscat,long ip,float x,float ph0, float th0, float rx, float ry, float rz, int inter,float op_time,long flag_zero)
void record_data_pos(long Np)
{

	FILE *fout;
	char OutFile[NCHARMAX];
	long i;
//	sprintf(sp,"%ld",ip);
	strcpy(OutFile,OutShort);
	
	printf("Output file : %s\n",OutFile);
	fflush(stdout);

	Np = (Np < NPhotons) ? Np : NPhotons;

	if (fout = fopen(OutFile,(char *) "w"))
	{

		fwrite(&NPhotons,sizeof(long),1,fout);
		fwrite(&Np,sizeof(long),1,fout);

		printf("Np = %ld\n",Np);
/*		for (i = 0; i<NPhotons; i++)
		{
			fwrite(&XArr[i],sizeof(float),NPhotons,fout);
			fwrite(&NscatArr[i],sizeof
*/		
		fwrite(InterArr,sizeof(int),Np,fout);
		fwrite(XArr,sizeof(float),Np,fout);
		fwrite(NscatArr,sizeof(int),Np,fout);
// Los siguientes dos sirven principalmente para crear mapas de brillos superficial, puedes comentarlos.
		fwrite(PosArr,sizeof(float),Np*3,fout);
		fwrite(AngleArr,sizeof(float),Np*2,fout);
		fwrite(X0Arr,sizeof(float),Np,fout);
/*
		fwrite(&nscat,sizeof(long),1,fout);
		fwrite(&x,sizeof(float),1,fout);
		fwrite(&ph0,sizeof(float),1,fout);
		fwrite(&th0,sizeof(float),1,fout);
		fwrite(&rx,sizeof(float),1,fout);
		fwrite(&ry,sizeof(float),1,fout);
		fwrite(&rz,sizeof(float),1,fout);
//		fprintf(fout," %ld\t %f\t %f\t %f\t %f\t %f\t %f\t %s\n",nscat,x, \
//								  ph0,th0,rx,ry,rz,itype );
*/

		fclose(fout);
	}
	else
	{
		printf("File %s could not be opened. fout = %d\n Exiting\n\n",OutFile,fout);
		fclose(fout);
		exit(0);
	}
}
		
