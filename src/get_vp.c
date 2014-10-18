/* This function receives a random number between 0<xi<1
and the frequency of the photon, 
and returns a value for the velocity of the atom parallel to 
the photon's direction. 
This is done by interpolating values 3 times. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

void get_vp(float *xlist_vp,float *vplist)
{
	FILE *flist;
	FILE *fvp;
	char str_x[9],buf[9],vpfile_low[NCHARMAX],vpaux[NCHARMAX];
	char x_str_low[9],x_str_high[9];
	long j,i;
	
    char vp_list[NCHARMAX];
	char vpfile[NCHARMAX];
	float xl,xl_aux,xi_list,u_list,xi,upar,xlow,xhigh;
	float xi0,upar0,xi1,uplow,upar1,x1,uphigh;

	strcpy(vp_list,tabledir);
	strcat(vp_list,"vp/t_1.0/vp_xlist");
		
	strcpy(vpfile,tabledir);
	strcat(vpfile,"vp/t_1.0/vp_x");


	if (flist = fopen(vp_list,"r"))
	{
	i = 0;
		while (!feof(flist))
		{
			fscanf(flist,"%s",str_x);
//			xl_aux = atof(buf);
			
			xlist_vp[i] = atof(str_x);
/*
			if (xl_aux <= x)
			{
				xlow = xl_aux;
				strcpy(x_str_low,buf);
			}
			else
			{
				xhigh = xl_aux;
				strcpy(x_str_high,buf);
				break;
			}
		}
*/
		i++;	
		}
	}
	else
	{
		printf("GET_VP.C: File %s (vp_list)  could not be opened. Exiting...\n\n\n",vp_list);
		exit(0);
	}
	nxlist = i-1;
	fclose(flist);

	strcpy(vpaux,vpfile);

	for (j = 0; j < nxlist; j++)
	{
		strcpy(vpfile,vpaux);		
		i = 0;

		sprintf(str_x,"%.3f",xlist_vp[j]);
		
		strcat(vpfile,str_x);
//		strcat(vpfile_high,x_str_high);
		if (fvp = fopen(vpfile,"r"))
		{
			fgets(buf,100,fvp);
			while(!feof(fvp))
			{
				fscanf(fvp,"%f %f",&vplist[2*(NTAB1*j + i)],&vplist[2*(NTAB1*j+i)+1]);
/*
			if (xi_list <= ran)
			{
				xi0 = xi_list;	
				upar0 = u_list;
			}
			else
			{
				xi1 = xi_list;
				upar1 = u_list;
				break;
			}
*/
			i++;
			}	
		}
		else
		{
			printf("GET_VP.C: File %s (vpfile) could not be opened, exiting...\n\n\n",vpfile);
			exit(0);
		}
	fclose(fvp);
	nvplist = i-1;
	}
/*
	uplow = ((upar1 - upar0)/(xi1 - xi0))*(x-x1) + upar1;

	if (fvp = fopen(vpfile_high,"r"))
	{
		fgets(buf,100,fvp);
		while(!feof(fvp))
		{
			fscanf(fvp,"%f %f",&xi_list,&u_list);
			if (xi_list <= ran)
			{
				xi0 = xi_list;	
				upar0 = u_list;
			}
			else
			{
				xi1 = xi_list;
				upar1 = u_list;
				break;
			}
		}	
	}
	else
	{
		printf("GET_VP.C: File %s (vpfile_high) could not be opened, exiting...\n\n\n",vpfile_high);
		exit(0);
	}

	fclose(fvp);
	uphigh = ((upar1 - upar0)/(xi1 - xi0))*(x-x1) + upar1;

	upar = ((uphigh - uplow)/(xhigh - xlow))*(x - xhigh) + uphigh ;


#if DEBUG

	dprintf(ran);
	dprintf(upar);
	dprintf(uplow);
	dprintf(uphigh);
	exit(0);

#endif

	return upar;
*/
}	
