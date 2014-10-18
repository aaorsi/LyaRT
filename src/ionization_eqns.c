/* Functions used to include the effect of the UV background in the neutral hydrogen profile 
   defined in define_geometry.c. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"


float Gamma_HI(float z)	// Photoionization rate [s-1] contribution of the UV background. No self-shielding.
{
	FILE *fd;	
	long i;
	float zlow,zhigh,glow,ghigh,g_,z_;
	float Gamma_z;
	char buf[NCHARMAX];
	zlow = 0.;
	dprints(GHIFile);
	if (fd = fopen(GHIFile,"r"))
	{
		fgets(buf,100,fd);		
		fgets(buf,100,fd);		
		i = 0;
		while(!feof(fd))
		{
			fscanf(fd,"%f %f\n",&g_,&z_);
			if(z_ < z)
			{
				zlow = z_;		
				glow = g_;	
			}
			else
			{
				zhigh = z_;
				ghigh = g_;
				break;
			}
			i++;
		}
		Gamma_z = (ghigh - glow)/(zhigh - zlow) * (z - zhigh) + ghigh;	
//		printf("Gamma_z_thin: zlow %f zhigh %f\n",zlow,zhigh);
	}
	else
	{
		printf("Cannot open file %s\n Exiting...",GHIFile);
		exit(0);
	}
	fclose(fd);
	return(Gamma_z);
}	


float beta_HI(float T) // collisional ionization rate.
{
	
	double T_HI,T5,rate;

	T_HI = 1.57807e5;
    T5 = T * 1.0e-5;

    rate = 5.85e-11 * sqrt(T) / (1.0 + sqrt(T5)) * exp(-T_HI/T);
	return(rate);
}


void set_spec(float redshift, float * uv_z, float *uv_logenergy, spec *minispec)
{
/* This function will read the UV spectrum at different redshifts
   and it will store a mini-spectrum, to be used to compute \Gamma_HI
*/
	FILE *fz,*fe,*fs;	
	float  flux_low[NEN], flux_high[NEN];
	float zlow,zhigh;
	int izlow,izhigh;
	float flow,fhigh;
	double interp_flux;
	char delim[] = "\t";
	char *res = NULL;
	char buf[1000];
	long i,i_;

	for (i = 0;i<NEN;i++)
	{
		flux_low[i] = 0.;
		flux_high[i] = 0.;
	}

/*
char str[] = "now # is the time for   all #   good men to come to the # aid of their country";
res = strtok( str, delim );
while( res != NULL ) {
    printf( "result is \"%s\"\n", res );
    res= strtok( NULL, delim );
}
	exit(0);
*/
	
	printf("Reading %s\n",UVzFile);
	if (fz = fopen(UVzFile,"r"))
	{
		i = 0;
		while(!feof(fz))
		{
			fscanf(fz,"%f\n",&uv_z[i]);
			i++;
		}
		UV_Nz = i-1;
		dprinti(UV_Nz);
	}
	else
	{
		printf("File %s cannot be open.\n Exiting...",UVzFile);
		exit(0);
	}
	fclose(fz);

	printf("Reading %s\n",UVeFile);
	if (fe = fopen(UVeFile,"r"))
	{
		i = 0;
		while(!feof(fe))
		{
			fscanf(fe,"%f\n",&uv_logenergy[i]);
			i++;
		}
		UV_Ne = i-1;
		dprinti(UV_Ne);
	}
	else
	{
		printf("File %s cannot be open.\n Exiting...",UVeFile);
		exit(0);
	}
	fclose(fe);

	i = 0;
	while (uv_z[i] < redshift)
		i++;

	zlow = uv_z[i-1];
	izlow = i-1;
	zhigh = uv_z[i];
	izhigh = i;
	
	printf(" redshift %f < %f < %f\n",zlow,redshift,zhigh);
	dprinti(i);
		

	if(uv_z[UV_Nz-1] < redshift)
	{
		printf("Warning: redshift %f outside redshift range in UV tables!\n",redshift);
		zlow = uv_z[UV_Nz-1];					
		zhigh = zlow;
		izlow = UV_Nz-1;
		izhigh = izlow;
	}

	printf("Reading UV spectrum from %s\n",UVspecFile);
	if (fs = fopen(UVspecFile,"r"))
	{
		i_ = 0;
		while(i_ < UV_Ne)
		{
			fgets(buf,1000,fs);
			res = strtok(buf,delim);
			i = 0;
			while (i < izlow)
			{
				res = strtok( NULL,delim);
//				printf("flux[i] %f i %d izlow %d\n",atof(res), i, izlow);
				i++;
			}
			flux_low[i_] = atof(res);
			res = strtok( NULL ,delim);
			flux_high[i_] = atof(res);
			i_++;
		}
	}
	else
	{	
		printf("File %s cannot be opened \n Exiting...",UVspecFile);
		exit(0);
	}	


	fclose(fs);

	printf("UV_spectrum: zlow %f zhigh %f\n",zlow,zhigh);

	i_ = 0;
	for (i = 0; i<UV_Ne; i++)
	{
		if ( uv_logenergy[i] >= minlogryd && \
			 uv_logenergy[i] <= maxlogryd)
		{
			minispec[i_].logryd = uv_logenergy[i];
			interp_flux = (flux_high[i] - flux_low[i])/(zhigh - zlow) * \
								(redshift - zhigh) + flux_high[i];
			minispec[i_].flux = pow(10,interp_flux);
			minispec[i_].sigmaHI = HI_photo_cs_analytic(pow(10,minispec[i_].logryd));
			minispec[i_].sigma_ratio = minispec[i_].sigmaHI/HI_photo_cs_analytic(1.00);			
			i_++;
		}
	}
	Nspec = i_-1;
	printf("Nspec %d UV_Ne %d\n",Nspec,UV_Ne);
}

float HI_photo_cs_verner(float Ryd)
{
/*	HI photoionization cross section. 	
	Taken from Verner et al (1996). 
	Written in general form, so it's easily adaptable to other elements if neccesary.
*/
	double Eth 		= 13.6;
    double Emax 	= 5.0e4;
    double E0 		= 4.298e-1;
    double sig0 	= 5.475e4;
    double ya 		= 3.288e1;
	double P 		= 2.963;
	double yw 		= 0;
	double y0		= 0;
	double y1 		= 0;
	double x,y;
	float sigma;

	x = (Ryd / E0) - y0;
	y = sqrt(SQR(x) + SQR(y1));
	
	sigma = sig0 *SQR(x - 1) + SQR(yw) * pow(y,0.5*P-5.5) * pow(1 + sqrt(y/ya),-P);
	return(sigma);
}


double HI_photo_cs_analytic(double Ryd)
{
/* Analytical HI photoionization cross section
	(see the Osterbrock (2006) AGN book, chapter 2, pag. 20)
*/

    double sigma;
    double A0 = 6.30e-18;
    double eps,ieps,efac,num,den;

    eps = sqrt( Ryd - 1.0 );
    ieps = 1.0/eps;
    efac = 4.0 - ( 4.0 * atan(eps) ) * ieps;
    num = exp(efac);
    den = 1.0 - exp(-2.0 * Pi * ieps);
    sigma = (Ryd == 1.00) ? A0 : (A0 * num / (den *pow(Ryd,4) ) );

	return(sigma);
}

float gamma_shield(spec *minispec, double tau0NH)
// Function to compute \Gamma including the effect of self-shielding. 
{
    double gammaHI;
	long i;
	double *xarr, *yarr;
	xarr = (double*) malloc(Nspec * sizeof(double));
	yarr = (double*) malloc(Nspec * sizeof(double));
//	double xarr[Nspec], yarr[Nspec];

	for (i=0;i<Nspec;i++)
	{
		xarr[i] = minispec[i].logryd;
		yarr[i] = exp(-tau0NH * minispec[i].sigma_ratio) *minispec[i].flux * minispec[i].sigmaHI;
//		yarr[i] = minispec[i].flux * minispec[i].sigmaHI;
//		dprintd(yarr[i]);
	}


    gammaHI = (4. * Pi * log(10.0) / Planck ) * int_tabulated_trap_rule(xarr,yarr,Nspec);
	
	free(xarr);
	free(yarr);

	return(gammaHI);

}

float get_xnHII(double Gamma_z, double beta_T, int id, cell *CellArr)
{
	double y,P,Q,R;
	double x_nHII;	
	y = 0;
	Q = beta_T * CellArr[id].nH  - Gamma_z - (beta_T + alpha_A)*CellArr[id].nH*y;
    R = -(beta_T + alpha_A)*CellArr[id].nH;
    P = Gamma_z + beta_T*CellArr[id].nH*y;
    x_nHII = (CellArr[id].nH == 0) ? 0. : (- Q - sqrt(SQR(Q) - 4*R*P))/(2*R);
	return(x_nHII);
}


