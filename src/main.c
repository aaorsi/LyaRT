/* Code to compute the radiative transfer of Ly-alpha photons
using a Monte Carlo algorithm. */


// Writing since: 12/09/09 
// Last Update: (still writing v0)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"
int main(int argc, char *argv[])
{

	char HFile[NCHARMAX],tabledir_aux[NCHARMAX],Tstr[NCHARMAX];
	int first_time=0;
	float pmin,pmax;
	double tau_ini, x_medium;
	double pos3d[3];
	double s_nH,sx_cont,s_nd;
	double sq_mu,sq_mu2,sq_b,sq_a,sq_arg;
	double dust_fac = ext_zstar/zstar,P_H;
	long ip,i,j,k;
	const gsl_rng_type * T;
	time_t t0,t1,t2,t4,tc1,tc2;

	int NScatRec = 1e4;


	long ilow,ihigh;
	double x_1,x_0,H1,H0,xl,Hl,xlow,xhigh,rvp,vpl,vp0l,vp1l,ran0l,ran1l,vp0h,vp1h,ran0h,ran1h,vplow,vphigh;
	float tot_tab,used_tab;
//
	double Vth=sqrt(2*KB*Temp/MP);
	int upcone;

	c = 299792.482;   // km s-1
	Inv_c = 1./c;	
	Kth=Vth/(c*100000.0);
	

	gsl_rng *r;
	gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
	
   if(argc < 3)
    {
        printf(" Input parameter missing. \n ");
		printf(" Usage: > LyaRt [ParFile] [SEED] \n Exiting...");
        exit(0);
    }


/*	Paramter file is read, it defines most of the variables
	found in allvars.h */

	strcpy(ParFile,argv[1]);
	Seed0 = atoi(argv[2]);
	default_parameters();
#ifndef SILENT
	dprints(OutShort);
	dprints(OutLong);
#endif

	read_parameters(ParFile);
#ifndef SILENT
	dprints(OutShort);
	dprints(OutLong);
#endif

//	The following are definitions specific to particular geometries. 

	if (strcmp(GeomName,"Wind2")==0 || strcmp(GeomName,"Wind3")==0) 
	{	
		dr = (RSphere-R_inner) / (NCells-1);
		NCellSphere = R_inner/dr + 1.0;
		NCells += NCellSphere;
	}

	if (strcmp(GeomName,"Wind4")==0)
	{	
		NCellSphere		= 100;
		dr = R_Static / (NCellSphere-1);
		NCells = (long) ((RSphere - R_inner)/dr + 1.0);
		
		NCells += NCellSphere+1;
		if (NCells > 1e4)
			printf("WARNING: NCells = %ld value too high\n",NCells);
	}		
		
	if (strcmp(GeomName,"ThinShell2")==0)
	{
		NCellSphere = 500;
		dr1	= R_inner / (NCellSphere-1);
		dr2 = (RSphere-R_inner)/(NCells-1);
		NCells+= NCellSphere;
	}

//	---

	P		= (photon*) malloc(NPhotons * sizeof(photon));
	CellArr = (cell*) malloc(NCells * sizeof(cell));
	define_geometry(GeomName,CellArr);

	gsl_rng_set(r,Seed0);	//Initializing the random generator with seed Seed0?

	//Finding range of emission probability in cells

	pmin = CellArr[0].p_lya;
	pmax = pmin;

	for (i=1;i<NCells;i++)
	{	
		if(CellArr[i].p_lya < pmin)
			pmin = CellArr[i].p_lya;
		if(CellArr[i].p_lya > pmax)
			pmax = CellArr[i].p_lya;
	}

	strcpy(tabledir_aux,tabledir);	
	strcat(tabledir_aux,sxfile);
	strcpy(HFile,tabledir_aux);	

	if (Temp >= 0)
		sprintf(Tstr,"%.2f",Temp);
	else
		sprintf(Tstr,"%.1f",Temp);

	if (Temp > 0 && Temp < 1)
		sprintf(Tstr,"%.3f",Temp);

	strcat(HFile,Tstr);

#ifndef SILENT
	printf(" Storing tabulated data... \n");
#endif

	HList = (double*) malloc(NTAB2*sizeof(double)*2);
	DipList = (double*) malloc(NTAB1*2*sizeof(double));
	HGList = (float*) malloc(NTAB2*2*sizeof(float));

#ifndef HAPPROX	
	get_H(HFile,HList);
#ifndef SILENT	
	printf("GET_H read\n");
#endif

#endif

	get_dipolar(DipList);	
#ifndef SILENT
	printf("GET_DIPOLAR read\n");
#endif	
	get_HG(HGList);	
#ifndef SILENT
	printf("GET_HG read\n");
#endif

	// Print parameters read
#ifndef SILENT
	print_parameters(CellArr);
#endif

#ifdef TEST_RII
	printf("TESTING REDISTRIBUTION FUNCTION for x0 = %f\n",x_test);
#endif

	//Loop over NPhotons
#ifndef SILENT
	printf(" Loop over %ld photons started...\n",NPhotons);
#endif

	cx = 0.;
	cy = 0.;
	cz = 0.;

	nu0 = 2.47e15;
	
	gtype = (strcmp(GeomName,"HomSlab")==0) ? 0 : 1;
	if (strcmp(GeomName,"HomSphere")==0 && NCells == 1) gtype = 2;
#ifndef SILENT
	dprinti(gtype);	
#endif	
	fflush(stdout);
	nout = 0;
	
#ifdef TCONST
	TPar = sx_const * pow(CellArr[0].T,-0.5);
#endif	

	XArr 		= (float *) malloc(NPhotons*sizeof(float));
	X0Arr           = (float *) malloc(NPhotons*sizeof(float));
	InterArr	= (int *) malloc(NPhotons*sizeof(int)); 
	NscatArr	= (int *) malloc(NPhotons*sizeof(int));
	BscatArr	= (int *) malloc(NPhotons*sizeof(int));
#ifdef GETPOSDIR
	PosArr		= (float *) malloc(NPhotons*3*sizeof(float));
	AngleArr	= (float *) malloc(NPhotons*2*sizeof(float));
//	PosArrLong	= (float *) malloc(NPhotons*NScatRec*3*sizeof(float));
#endif

	for (i = 0; i< NPhotons; i++)
	{
		XArr[i]   = 0;
		X0Arr[i]   = 0;
		InterArr[i] = -1;
		NscatArr[i] = -1;
		BscatArr[i] = -1;
#ifdef GETPOSDIR
		PosArr[3*i] 	= 0;
		PosArr[3*i+1] 	= 0;
		PosArr[3*i+2]	= 0;
		AngleArr[2*i]	= 0;
		AngleArr[2*i+1]	= 0;
/*		for (j = 0;j<NScatRec;j++)
		{
			PosArrLong[0+ 3*i + 3*NPhotons*j] = 0;
			PosArrLong[1+ 3*i + 3*NPhotons*j] = 0;
			PosArrLong[2+ 3*i + 3*NPhotons*j] = 0;
		}
*/
#endif
	}


	if (b == 0.)
	{
		a_par = 4.693e-4 * sqrt(1./CellArr[idc].T);
		vth = 12.85 * sqrt(CellArr[idc].T);
	}
	else 
	{	
		a_par = 4.693e-4 * sqrt(1./CellArr[idc].T);

//			a_par = 4.693e-4 * (12.85/b);
//			a_par = 4.7e-4 * (12.85/b);
//			vth	  = ;

		vth = 12.85 * sqrt(CellArr[idc].T);
	}

	if (VarXCrit==1)
	{ 
        x_medium = xp0 - vmax/vth;
		tau_ini	 = sx_const *  ColDens * pow(CellArr[idc].T,-0.5) * voigt(HList,x_medium);

		if (tau_ini < 1e4)
		  xcrit = 3;
		if (tau_ini >1e4 && tau_ini < 1e6)
		  xcrit = 6;
		if (tau_ini > 1e6)
		  xcrit = 10;
		//printf("1 xcrit: %f \n",xcrit);
	}

	
	(void) time(&t0);

  for (ip=0;ip < NPhotons;ip++)
	{
        nscat = 0;
        //printf("Nuevo foton\n");
		
/*		if( ip > 999 & ip % 1000 == 0)
	  	{
	   		printf("Photon %ld\n",ip);
	    	if( first_time==0 )
	     	{
				first_time=1;
				printf("Estimated time to finish: %f [min] \n",1000*op_time * (NPhotons - ip) );
	      	}
	  	}
*/

//		Initialise photon's frequency, direction and position.

		xi0 = gsl_rng_uniform (r);
		xi1 = gsl_rng_uniform (r);
		xi2 = gsl_rng_uniform (r);

			
		xi3 = gsl_rng_uniform (r);
		xi4 = gsl_rng_uniform (r);
		xi5 = gsl_rng_uniform (r);

		init_photon(P,ip, xi0, xi1, xi2);		  

		dnd = (vth*nu0)/c;
		flag_zero = 0;
		
		(void) time(&t1);
		i = 0;
		xi1 = gsl_rng_uniform (r);
		xi2 = gsl_rng_uniform (r);
		xi3 = gsl_rng_uniform (r);
		idc_old = 0;
		while(idc != -1)
		{
			t_0 = - log(gsl_rng_uniform (r));
		
			EscCond = 0;

			P[ip].ni = sin(th0)*cos(ph0);
			P[ip].nj = sin(th0)*sin(ph0);
			P[ip].nk = cos(th0);

			while(t_0 > 0.)
			{  	
#ifdef TIMELIMIT
		(void) time(&t4);
		t4 = (float) (t4-t0)/60.;
		if (t4 > MAXTIME)
		{
			printf("MAXTIME reached, forcing output and exit\n");
			fflush(stdout);
#ifdef GETPOSDIR
			record_data_pos(ip);
#ifndef SILENT
			printf("File %s written\n",OutShort);
#endif
			free(PosArr);
			free(AngleArr);
#else
			record_data_short();
#endif
			exit(0);
		}
#endif
						

				// Shortcut to compute crossing of empty cells, where no optical depth (t_0) is *used*.
				empty_cells(P,ip);			
				BscatArr[ip] = flag_zero;
				H_x = voigt(HList, P[ip].xp);
	
				if isnan(H_x) 
					printf("H_x is nan\n");

//				H_x = 0.;     // H_x disables H scattering. 

                //Siddhartha walawala change to....

                double nH_sid , r0_sid , px0_sid , py0_sid , pz0_sid , th_c_sid ;
                int SidCond;

                nH_sid = CellArr[idc].nH;

                r0_sid = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z));

                px0_sid = P[ip].x/r0_sid;
                py0_sid = P[ip].y/r0_sid;
                pz0_sid = P[ip].z/r0_sid;

                
                th_c_sid = acos( pz0_sid );

                if (strcmp(GeomName,"Anti_Wind")==0 && idc != 0)
                {    
                SidCond = (pz0_sid > 0) ? ( (th_c_sid > app_angle ) ? 1 : 0) : \ 
                          ((th_c_sid < Pi - app_angle ) ? 1 : 0);

                    if ( SidCond == 1 )
                    {
                        nH_sid = nH_sid * f_Anti_wind;
                    }
 
                }


                if (strcmp(GeomName,"Slab_plus_bicone")==0 && idc != 0)
                {    
                SidCond = (pz0_sid > 0) ? ( (th_c_sid > app_angle ) ? 1 : 0) : \ 
                          ((th_c_sid < Pi - app_angle ) ? 1 : 0); 

                    if ( SidCond == 1 ) 
                    {
                        nH_sid = HOMOGENIUOS_DENSITY_SID ;
                        CellArr[idc].vbulk[0] = V_dense_region_sid ; 
                    }
 
                    if ( SidCond == 0 ) 
                    {
                        CellArr[idc].vbulk[0] = V_transparent_region_sid ;
                    }
                    

                }
     

                //printf("%f , %f , %f\n"  , th_c_sid , nH_sid , CellArr[idc].vbulk[0]  );


	
#ifdef TCONST 
				
				//s_nH = TPar * H_x*CellArr[idc].nH;
				s_nH = TPar * H_x* nH_sid ;

#else	
				s_nH = sx_const * pow(CellArr[idc].T,-0.5)*H_x* nH_sid ;
				//s_nH = sx_const * pow(CellArr[idc].T,-0.5)*H_x*CellArr[idc].nH;
#endif

#ifdef TAUGUIDERDONI
				s_nd = Ext_Ratio * pow(CellArr[idc].z/zstar,spar) * nH_sid /NHConst;
				//s_nd = Ext_Ratio * pow(CellArr[idc].z/zstar,spar) * CellArr[idc].nH/NHConst;
#else
				s_nd =  ext_zstar * CellArr[idc].z * nH_sid ;
				//s_nd =  ext_zstar * CellArr[idc].z * CellArr[idc].nH;
                //Siddhartha walawala change to.... HERE. I changed this
#endif
	
				s_sum = s_nH + s_nd;
				s = t_0 / s_sum;

//				radius = sqrt( SQR(P[ip].x) + SQR(P[ip].y)+ SQR(P[ip].z));
				
				// The position of the photon after consuming its optical depth is computed here:
				s_ = crossing_cells(P,gtype,ip);
	            
//				rx0 = P[ip].x;
//				ry0 = P[ip].y;
//				rz0 = P[ip].z;

				ni = P[ip].ni;
				nj = P[ip].nj;	
				nk = P[ip].nk;	

			}

	        if (EscCond == 1)
   			{
			
#ifndef TEST_RII	
				idc = -1;
#else
				idc = 0;
#endif

	   			if (idc == -1)
   	   			{
escape:
//					if (idc == idc_old)	
						P[ip].xp += (P[ip].ni*vbulk_x + P[ip].nj*vbulk_y + P[ip].nk*vbulk_z)/vth;
					x = P[ip].xp;
//				x = c*(g*nup - nu0)/(vth*nu0) - g*(nup/nu0)*(ni*vbulk_x + nj*vbulk_y + nk*vbulk_z)/vth;

//				dprintd(x);
//				printf("Photon has escaped in %ld scatterings\n",nscat);
// ................................

					inter = 4;
					(void) time(&t2);
					op_time = (float) (t2-t1)/60.;
//			printf("Total calculation took %f minutes.\n",op_time);
					if (strcmp(OutMode,"Long")==0)
					{
						printf("\n");
						record_data_long(nscat,ip,x,P[ip].x,P[ip].y,P[ip].z,r0,upar,uper1,uper2,H_x,P[ip].xp,inter);
					}
					else
					{
//					record_data_short(nscat,ip,x,rxf,ryf,rzf,radius,H_x,inter,op_time,flag_zero);
						XArr[ip] = P[ip].xp;
						X0Arr[ip]=xp0;
						InterArr[ip] = inter;
						NscatArr[ip] = nscat;
					
#ifdef GETPOSDIR
						PosArr[3*ip] 	= P[ip].x;
						PosArr[3*ip+1]	= P[ip].y;
						PosArr[3*ip+2]	= P[ip].z;
						AngleArr[2*ip]	= P[ip].th;
						AngleArr[2*ip+1]= P[ip].ph;	
//						record_data_pos(nscat,ip,x,ph0,th0,rxf,ryf,rzf,inter,op_time,flag_zero);
//#else
//						record_data_short(nscat,ip,x,ph0,th0,radius,inter,op_time,flag_zero);
#endif
					
					}
					nout++;
					goto end;
					break;
				}
   			}
   			else
   			{

/*				rx0 = rxf;
				ry0 = ryf;
				rz0 = rzf;
*/			
			}

			P_H = s_nH / s_sum;
			xi1 = gsl_rng_uniform (r);
	
			if (strcmp(IncDust,"Yes")==0)
				inter = (xi1 <= P_H)? 1 : 2 ;
			else
				inter = 1;
		
			switch(inter)	// inter = 1 means interacting with hydrogen, inter = 2 means dust.
			{
				case  1:
	
#ifdef USEREJECTION
			upar = -999.0;
			i = 0;
			do
			{
				xi0 = gsl_rng_uniform (r);
				xi1 = gsl_rng_uniform (r);
				xi2 = gsl_rng_uniform (r);
				upar = vp_rejection(P[ip].xp,a_par,xi0,xi1,xi2);
				i++;
			} while(upar == -999.0);
//	if (i > 100000)
//						printf("VP_REJECTION TOOK %d tries to get upar = %f, x = %f,  nscatter = %d\n" \
						 ,i,upar,xp,nscat);
#endif
			
				xi0 = gsl_rng_uniform (r);
				xi1 = gsl_rng_uniform (r);
				xi2 = gsl_rng_uniform (r);
				xi3 = gsl_rng_uniform (r);
				xi3 = (xi3 == 0.0) ? gsl_rng_uniform (r) : xi3;
				xi4 = gsl_rng_uniform (r);
				xi5 = gsl_rng_uniform (r);
				xi6 = gsl_rng_uniform (r);
				xi7 = gsl_rng_uniform (r);
				
				scattering_hydrogen(P,ip);
				break;
				
				case 2:
//de aqui, hasta el endif --> sidd
#ifdef USEREJECTION
			upar = -999.0;
			i = 0;
			do
			{
				xi0 = gsl_rng_uniform (r);
				xi1 = gsl_rng_uniform (r);
				xi2 = gsl_rng_uniform (r);
				upar = vp_rejection(P[ip].xp,a_par,xi0,xi1,xi2);
				i++;
			} while(upar == -999.0);
//	if (i > 100000)
//						printf("VP_REJECTION TOOK %d tries to get upar = %f, x = %f,  nscatter = %d\n" \
						 ,i,upar,xp,nscat);
#endif
// todo este ifdef lo ha puesto sidd

				xi1 = gsl_rng_uniform (r);
				xi2 = gsl_rng_uniform (r);
				xi3 = gsl_rng_uniform (r);                          // sidd
				xi3 = (xi3 == 0.0) ? gsl_rng_uniform (r) : xi3;     // sidd
				xi4 = gsl_rng_uniform (r);                          // sidd
				xi5 = gsl_rng_uniform (r);                          // sidd
				xi6 = gsl_rng_uniform (r);                          // sidd
				xi7 = gsl_rng_uniform (r);                          // sidd

				dust_interaction(P,ip);
				if (end_syg ==1)
					goto end;
				break;
			}
	
//			dprinti(nscat);

/*
			if (nscat == 1000000 || nscat == 10000000 )
			{
				printf("********************************************\n");
				printf("Number of scatterings so far: %ld\n",nscat);
				printf("Photon's location %f %f %f\n",rxf,ryf,rzf);
				dprintf(x);
				dprintf(xp);
				dprintf(th0);
				dprintf(ph0);
			}
*/
				
#ifdef WRITEALL
			if (ip < 1000)
			{
				P[ip].xp += (P[ip].ni*vbulk_x + P[ip].nj*vbulk_y + P[ip].nk*vbulk_z)/vth;
				record_data_long(nscat,ip,x,P[ip].x,P[ip].y,P[ip].z,r0,upar,uper1,uper2,H_x,\
				P[ip].xp,inter);
			}
#endif
			if (gtype ==1)
			{
				if (r0 < 0.99*R_inner)
				{
					printf("scatter %d of photon %d occurs within empty zone. Something is not ok\n",nscat,ip);
					printf("rx/r_inner %f ry/r_inner %f rz/r_inner %f radius/r_inner %f\nidc %d\n",\
					P[ip].x/R_inner,P[ip].y/R_inner,P[ip].z/R_inner,r0/R_inner,idc);
				}
			}
			nscat++;

#ifdef TIMELIMIT
		(void) time(&t4);
		t4 = (float) (t4-t0)/60.;
		if (t4 > MAXTIME)
		{
			printf("MAXTIME reached, forcing output and exit\n");
			fflush(stdout);
#ifdef GETPOSDIR
			record_data_pos(ip);
#ifndef SILENT
			printf("File %s written\n",OutShort);
#endif
			free(PosArr);
			free(AngleArr);
#else
			record_data_short();
#endif
			exit(0);
		}
#endif
		}				
		

		end:
//	printf("done\n");


	

		if (strcmp(Set_Tolerance,"yes") == 0)
		{
//		if ((ip >= np_min && nout >= nout_max) || \
//			(ip >= np_max))

//			if ((ip >= np_min && nout >= nout_max) || 
			if ( 1 == 2 ) //)
//			(ip >= np_max && nout >= (int) nout_max/5) || \

			{
				printf("Max number of absorbed photons or limit on the number of photons reached,  exiting...\n");

#ifdef GETPOSDIR
#ifndef SILENT
				printf("Writing data\n");
#endif
				record_data_pos(ip);
#ifndef SILENT
				printf("File %s written\n",OutShort);
#endif
				free(PosArr);
				free(AngleArr);
#else
//				record_data_short();
#endif
				free(HList);
				free(DipList);
				free(HGList);
				free(CellArr);
				exit(0);
			}
		}		

#ifndef SILENT
	printf("%ld %f %d\n",ip,P[ip].xp,P[ip].nscat);
#endif

	}
#ifndef SILENT
	printf("Writing data \n");
#endif

#ifdef GETPOSDIR
	record_data_pos(ip);
#ifndef SILENT
	printf("File %s written\n",OutShort);
#endif

	free(PosArr);
	free(AngleArr);
#else
	record_data_short();
#endif
	free(XArr);
   	free(X0Arr);
	free(InterArr);
	free(NscatArr);

	free(HList);
	free(DipList);
	free(HGList);
	free(CellArr);
	return(0);
}





