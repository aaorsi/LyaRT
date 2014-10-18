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

	char ParFile[NCHARMAX],HFile[NCHARMAX],tabledir_aux[NCHARMAX],Tstr[NCHARMAX];
	int Seed0,inter,interd,gtype;
	int first_time=0;
	cell *CellArr;
	float pmin,pmax;
	double tau_ini, x_medium;
	double xi,xi0,xi1,xi2,xi3,xi4,xi5,xi6,xi7,xp0;
	double r0,rE,a_,g_,b_,s_,px,py,pz,s1,s2,s3,s4,s5,px0,py0,pz0;
	double pos3d[3];
	int EscCond,ObscuredCond;
	double x,t_0,t_used,t_av,H_x,tau0,s_sum,e_x,e_y,e_z,r_edge,xcrit;
	double X,xp;
	double s,op_time,sf,dnd,nup,nu0,nu_,vel;
	double vth,n_ko,nko0,nvpar,nvper,na,vpn;
	double ko_x,ko_y,ko_z;
	double rx0,ry0,rz0,rxf,ryf,rzf,radius,cx,cy,cz,acrit,aarg;
	double x0,y0,z0,th_;
	double xcw,th_aux,th_range[2],u0,u_par,g_u,accept,u_p1,u_p2;
	double s_nH,sx_cont,s_nd,g;
	double upar,uper1,uper2,utot,alpha;
	double k_p_par,k_p_per,mu,mu2,ko1,ko2,ko3,k_po_par,k_po_per,k_o_par,k_o_per,b_ang;
	double n1_i,n1_j,n1_k,n2_i,n2_j,n2_k,v_i,v_j,v_k,v_pi,v_pj,v_pk,v_po_i,v_po_j,v_po_k;
	double a1,a2,a3,vpar_i,vpar_j,vpar_k,vper_i,vper_j,vper_k,v_ko,uper,vper;
	double sq_mu,sq_mu2,sq_b,sq_a,sq_arg;
	double dust_fac = ext_zstar/zstar,P_H;
	long ip,i,j,k,idc_old,idc,nbins,nout;
	const gsl_rng_type * T;
	time_t t0,t1,t2,t4,tc1,tc2;

	long flag_zero;
	float *HGList;
	double *HList, *DipList;
	int NScatRec = 1e4;


	long ilow,ihigh;
	double x_1,x_0,H1,H0,xl,Hl,xlow,xhigh,rvp,vpl,vp0l,vp1l,ran0l,ran1l,vp0h,vp1h,ran0h,ran1h,vplow,vphigh;
	double xilow,xihigh,mulow,muhigh,hglow,hghigh,xabs,xpabs;
	float tot_tab,used_tab;
	double Inv_c,TPar,Arg1_uper,Arg2_uper,cosph0,sinph0,costh0,sinth0,cosalpha,sinalpha,Inv_sinalpha, \
			cosmu,sinmu,cosmu2,sinmu2,costh_,sinth_,Inv_kpoper,gcond,garg;
//
	double Vth=sqrt(2*KB*Temp/MP);
	double Kth,zeta,th_c;
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
	dprints(OutShort);
	dprints(OutLong);
	
	read_parameters(ParFile);
	dprints(OutShort);
	dprints(OutLong);

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

	printf(" Storing tabulated data... \n");

	HList = (double*) malloc(NTAB2*sizeof(double)*2);
	DipList = (double*) malloc(NTAB1*2*sizeof(double));
	HGList = (float*) malloc(NTAB2*2*sizeof(float));

#ifndef HAPPROX	
	get_H(HFile,HList);
	printf("GET_H read\n");
#endif

	get_dipolar(DipList);	
	printf("GET_DIPOLAR read\n");
	get_HG(HGList);	
	printf("GET_HG read\n");

	// Print parameters read
	print_parameters();

#ifdef TEST_RII
	printf("TESTING REDISTRIBUTION FUNCTION for x0 = %f\n",x_test);
#endif

	//Loop over NPhotons
	printf(" Loop over %ld photons started...\n",NPhotons);

	cx = 0.;
	cy = 0.;
	cz = 0.;

	nu0 = 2.47e15;
	
	gtype = (strcmp(GeomName,"HomSlab")==0) ? 0 : 1;
	if (strcmp(GeomName,"HomSphere")==0) gtype = 2;

	dprinti(gtype);	
	
	fflush(stdout);
	nout = 0;
	
#ifdef TCONST
	TPar = sx_const * pow(CellArr[0].T,-0.5);
#endif	

	XArr 		= (float *) malloc(NPhotons*sizeof(float));
	X0Arr           = (float *) malloc(NPhotons*sizeof(float));
	InterArr	= (int *) malloc(NPhotons*sizeof(int)); 
	NscatArr	= (int *) malloc(NPhotons*sizeof(int));
#ifdef GETPOSDIR
	PosArr		= (float *) malloc(NPhotons*3*sizeof(float));
	AngleArr	= (float *) malloc(NPhotons*2*sizeof(float));
	PosArrLong	= (float *) malloc(NPhotons*NScatRec*3*sizeof(float));
#endif

	for (i = 0; i< NPhotons; i++)
	{
		XArr[i]   = 0;
		X0Arr[i]   = 0;
		InterArr[i] = -1;
		NscatArr[i] = -1;
#ifdef GETPOSDIR
		PosArr[3*i] 	= 0;
		PosArr[3*i+1] 	= 0;
		PosArr[3*i+2]	= 0;
		AngleArr[2*i]	= 0;
		AngleArr[2*i+1]	= 0;
		for (j = 0;j<NScatRec;j++)
		{
			PosArrLong[0+ 3*i + 3*NPhotons*j] = 0;
			PosArrLong[1+ 3*i + 3*NPhotons*j] = 0;
			PosArrLong[2+ 3*i + 3*NPhotons*j] = 0;
		}
#endif
	}
	
	(void) time(&t0);
	for (ip=0;ip < NPhotons;ip++)
	{
	        nscat = 0;
		
		if( ip > 999 & ip % 1000 == 0)
		  {
		    printf("Photon %ld\n",ip);
		    if( first_time==0 )
		      {
			first_time=1;
			printf("Estimated time to finish: %f [min] \n",1000*op_time * (NPhotons - ip) );
		      }
		  }
		  

		flag_zero = 0;

		(void) time(&t1);

		i = 0;

		xi1 = gsl_rng_uniform (r);
		xi2 = gsl_rng_uniform (r);
		xi3 = gsl_rng_uniform (r);

		if (strcmp(GeomName ,"Biconical_Wind")==0)
		{
			
			xi0 = gsl_rng_uniform (r);
			upcone = (xi0 < 0.5) ? 1 : -1;	// This random number chooses the hemisphere
//			th0 = acos(upcone*(xi1*(cos(app_angle)-1)+1)); // theta is constrained to be drawn from 0,app_angle, or pi,pi-app_angle
			th0 = (upcone == 1) ? acos(xi1*(1-cos(app_angle))+cos(app_angle) ) : \
								  Pi - acos(xi1*(1-cos(app_angle))+ cos(app_angle));
		}
		else th0 = acos(2*xi1 - 1);

   		ph0 = 2*Pi*xi2;

		if (StartAtCenter== 1)
		{
//			idc = (int) (NCells-1)/2;
//			if (gtype == 1)

			idc = 0;
			idc_old = 0;	
			rx0 = 0.;			
			ry0 = 0.;			
			rz0 = 0.;		

  	    }
		else
		{
			if (DefinePosAndDir == 1)
			{
				rx0 = Ix0;
				ry0 = Iy0;
				rz0 = Iz0;
	
				th0	= Ith0;
				ph0 = Iph0;
	
				r0  = sqrt(SQR(rx0) + SQR(ry0) + SQR(rz0));
				
				idc = 0;
				idc_old = 0;
				for (i = 0;i<NCells;i++)
				{
					if (r0 <= CellArr[i].r + dr/2.)
					{
						idc = i;
						break;
					}
				}

//				dprinti(idc);
//				dprintf(r0);

			}
			else 				
			{

				if (gtype == 1)
				{
					xi1 = gsl_rng_uniform (r);
//				idc = (long) gsl_rng_uniform (r) * (NCellSphere-1);
					xi2  = (gsl_rng_uniform (r) / 3.);
					r0	= pow(3 * xi2,1./3) * R_inner;
					if (strcmp(GeomName ,"Biconical_Wind")==0)
					{
				
						xi0 = gsl_rng_uniform (r);
						upcone = (xi0 < 0.5) ? 1 : -1;	// This random number chooses the hemisphere
//						th0 = acos(upcone*(xi1*(cos(app_angle)-1)+1));

						th0 = (upcone == 1) ? acos(xi1*(1-cos(app_angle))+cos(app_angle) ) : \
								  Pi - acos(xi1*(1-cos(app_angle))+ cos(app_angle));
					}
					else th0 = acos(2*xi1 - 1);
    	  			ph0 = 2*Pi*gsl_rng_uniform (r) ;

					idc = (long) (r0 / dr);
					idc_old = 0;

					rx0	= r0 * sin(th0)*cos(ph0);
					ry0	= r0 * sin(th0)*sin(ph0);
					rz0 = r0 * cos(th0);
				
					if (idc >= NCellSphere)
					{
						printf("ERROR: Starting point outside inner radius\n Check!");	
						exit(0);
					}

/*				dprinti(idc);
				dprinti(NCellSphere);
				dprintd(r0);
				dprintd(R_inner);
*/
				}
				else
				{
					do {
						xi2 = gsl_rng_uniform (r) * (NCells-1);
						i++;
					   } while ( CellArr[ (long) xi2].p_lya > xi1);
	
		
					idc = (int) xi2;
					printf("%ld tries finding starting cell id %d for photon %ld\n",i,(int) xi2,ip);

				    rx0 = (gsl_rng_uniform (r) - 0.5)* xSize;
			   	    ry0 = (gsl_rng_uniform (r) - 0.5)* ySize;
   			 		rz0 = (gsl_rng_uniform (r) - 0.5)* zSize;
				}
			}
		}


//		b = sqrt(v_th^2 + v_turb^2). If defined, then vth = b

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
			
		dnd = (vth*nu0)/c;

//      printf("Initial angles: %f %f\n",th0 * 180/Pi,ph0*180/Pi);
//      printf("Initial position: (%f %f %f)\n",(rx0),(ry0),(rz0));

// 		Initial frequency

		xp = 0; // initial frequency in the rest frame of the fluid 
		xp0=xp;		
// Flat Continuum
#ifdef XCONTINUUM
		//dX=96.0;
		xi1 = gsl_rng_uniform (r);
		
		xp = (xi1-0.5)*dX;
		xp0=xp; 




		//printf("xp: %f xp0=%f \n",xp,xp0);
		//#endif
#endif		
		
#ifdef XINDIST
		
		xi1 = gsl_rng_uniform (r);
		xi2 = gsl_rng_uniform (r);
		
		xp	= x_mean + Vrot * sqrt(-2*log(xi1))*cos(2*Pi*xi2);
		xp0=xp; 
		dprintl(ip);
		dprintd(xp);
		printf("------");
		//#else
		//xp = 0; // initial frequency in the rest frame of the fluid
		//xp0=xp; 
		//printf("0 xp: %f xp0=%f \n",xp,xp0);
#endif 
		x = xp;
		xp0=xp; 
		//printf("1 xp: %f xp0=%f \n",xp,xp0);
		xabs = fabs(x);

		rxf 	= rx0;
		ryf 	= ry0;
		rzf 	= rz0;

		inter 	= 1;

	// AQUI:
		
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

		if (StartAtCenter == 1)
		 {

			r0 				= 1;
			radius 			= 0;

			vbulk_x			= 0;
			vbulk_y			= 0;
			vbulk_z			= 0;
			vbulk0_x		= 0;
			vbulk0_y		= 0;
			vbulk0_z		= 0;
			
			px0 			= 0.;
			py0 			= 0.;
			pz0 			= 0.;

		}
		else
		{

			px0 = rx0/r0;
			py0 = ry0/r0;
			pz0 = rz0/r0;

			radius = r0;

			vbulk_x		= CellArr[idc].vbulk[0] * px0;
			vbulk_y		= CellArr[idc].vbulk[0] * py0;
			vbulk_z		= CellArr[idc].vbulk[0] * pz0;
			vbulk0_x	= vbulk_x;
			vbulk0_y	= vbulk_y;
			vbulk0_z	= vbulk_z;
			
		}	

		while(idc != -1)
		{
			t_0 = - log(gsl_rng_uniform (r));
		
			EscCond = 0;

			ni = sin(th0)*cos(ph0);
			nj = sin(th0)*sin(ph0);
			nk = cos(th0);

			while(t_0 > 0.)
			{  	

				
			//  The following lines define the behaviour of the photons when crossing (practically) empty cells

				if (CellArr[idc].nH < 1e-20 && idc != -1)
				{  
				  
					if (flag_zero >0)	
					{
						garg = -px0*ni - py0*nj - pz0*nk;
						
						g_   = (garg >= 1.0) ? 0. : (garg <=-1.0) ? Pi : acos(garg);
						b_   = Pi - 2*g_;
						s_   = (g_ == 0 ) ? 2 * CellArr[0].r : r0 * sin(b_)/sin(g_);

						if ( flag_zero > 1e20)	
						{ 
							printf("+++++++++++++++++++++++++++++\n");	
							printf("WARNING: flag_zero>NCells\n");
							dprintd(xp);	
							dprintd(t_0);
							dprintd(radius);
							dprintd(r0);
							dprintd(rE);
							dprintd(px0);
							dprintd(py0);
							dprintd(pz0);
							dprintd(ni);
							dprintd(nj);
							dprintd(nk);
							dprintd(g_);
							dprintd(b_);
							dprintd(s_);
							dprintl(idc);
							dprintl(flag_zero);
							dprintl(nscat);
							dprintd(rx0);
							dprintd(ry0);
							dprintd(rz0);
							dprintd(dr);
							printf("EXITING...\n");
							exit(0);
					   	}  

						rx0 += s_*ni;
						ry0 += s_*nj;
						rz0 += s_*nk;

						r0   = sqrt( SQR(rx0) + SQR(ry0)+ SQR(rz0) );
						idc = (idc_old > idc && idc != 0) ? CellArr[idc].neig[0] : CellArr[idc].neig[1];
							
				}

					else
		
					{		
						idc = (idc_old > idc && idc != 0) ? CellArr[idc].neig[0] : CellArr[idc].neig[1];
						rx0 = CellArr[idc].r*ni;
						ry0 = CellArr[idc].r*nj;
						rz0 = CellArr[idc].r*nk;
						r0 = sqrt( SQR(rx0) + SQR(ry0)+ SQR(rz0) );
					}

					px0 = rx0/r0;
					py0 = ry0/r0;
					pz0 = rz0/r0;

					vbulk0_x = vbulk_x;
					vbulk0_y = vbulk_y;
					vbulk0_z = vbulk_z;

					vbulk_x = px0 * CellArr[idc].vbulk[0];
					vbulk_y = py0 * CellArr[idc].vbulk[0];
					vbulk_z = pz0 * CellArr[idc].vbulk[0];
		
					xp += (ni * (vbulk0_x-vbulk_x) + nj*(vbulk0_y-vbulk_y) + nk*(vbulk0_z-vbulk_z) )/vth;

					flag_zero++;
					
				}

				H_x = voigt(HList, xp);
//				H_x = 0.;     // H_x disables H scattering. 
	
#ifdef TCONST 
				
				s_nH = TPar * H_x*CellArr[idc].nH;

#else	
				s_nH = sx_const * pow(CellArr[idc].T,-0.5)*H_x*CellArr[idc].nH;
#endif

#ifdef TAUGUIDERDONI
				s_nd = Ext_Ratio * pow(CellArr[idc].z/zstar,spar) * CellArr[idc].nH/NHConst;
#else
				s_nd =  ext_zstar * CellArr[idc].z * CellArr[idc].nH;
#endif

				s_sum = s_nH + s_nd;
				s = t_0 / s_sum;

                rxf = rx0 + s*ni;
 		        ryf = ry0 + s*nj;
       		    rzf = rz0 + s*nk;	

				radius = sqrt( SQR(rxf) + SQR(ryf)+ SQR(rzf) );

				switch(gtype)
				{
					case 1:
						
						r0 = sqrt( SQR(rx0) + SQR(ry0)+ SQR(rz0) );

			
						if (idc > 0)
						{

							px0 = rx0/r0;
							py0 = ry0/r0;
							pz0 = rz0/r0;

							garg = -px0*ni - py0*nj - pz0*nk;
							gcond = fabs(garg);

							g_   = (garg >= 1.0) ? 0. : ((garg <= -1.0) ? Pi : acos(garg));
							
							aarg = CellArr[idc].r/r0;
							acrit = (aarg >= 1.0) ? Pi/2. : ( (aarg <= -1.00) ? -Pi/2. : asin(aarg) );
							rE = (g_ >= acrit ) ? CellArr[idc].r + dr : CellArr[idc].r;


						}
						else
						{
		
							px0 = (r0 == 0 ) ? 0 : rx0/r0;
							py0 = (r0 == 0 ) ? 0 : ry0/r0;
							pz0 = (r0 == 0 ) ? 0 : rz0/r0;

							garg = -px0*ni - py0*nj - pz0*nk;
								
							g_   = (garg >= 1.0) ? 0. : ((garg <=-1.0) ? Pi : acos(garg));

							acrit	=-20;
							rE = CellArr[idc+1].r;

						}

						if  ((fabs(garg) < 1.001 && fabs(garg) > 0.999) || (idc == 0 && garg == 0))
						{	
//							printf("condition reached \n");

							a_ = 0.;
							b_ = 0.;
							s_ = (idc == 0) ? CellArr[1].r : dr;
						} 
						else 
						{

							aarg = r0 * sin(g_)/rE;
							aarg = (aarg >= 1) ? Pi/2. : ((aarg <=-1.0) ? -Pi/2. : asin(aarg));
							b_ = (g_ >= acrit ) ? aarg : Pi - aarg;
							a_ = Pi - (b_ + g_);
							a_ = (a_ < 0 ) ? 0. : a_;	
							
							s_ = r0 * sin(a_)/sin(b_);

						}

// strong debug check

						if ( isnan(xp) || isnan(t_0) || isnan(s_))	
						{
							printf("+++++++++++++++++++++++++++++\n");	
							printf("Something could be wrong here -- at the end of the loop\n");
							dprintd(xp);	
							dprintd(t_0);
							dprintd(radius);
							dprintd(r0);
							dprintd(rE);
							dprintd(px0);
							dprintd(py0);
							dprintd(pz0);
							dprintd(ni);
							dprintd(nj);
							dprintd(nk);
							dprintd(g_);
							dprintd(b_);
							dprintd(s_);
							dprintd(a_);
							dprintd(gcond);
							dprintd(garg);
							dprintd(acos(garg));
							dprintl(idc);
							dprintl(flag_zero);
							dprintl(nscat);
							dprintd(acrit);
							exit(0);
						}


						if (fabs(s_) > s)
						{		
							t_0 = -1;
							 break;
						}

						t_used = s_sum * s_;
						t_0 = t_0 - t_used;

						if (t_0 <= 0)
						{
							t_0 = -1;
							break;
						}

						rx0 += s_*ni;
						ry0 += s_*nj;
						rz0 += s_*nk;

						idc_old = idc;
						idc = (g_ >= acrit ) ? CellArr[idc].neig[1] : CellArr[idc].neig[0];
						r0 = sqrt( SQR(rx0) + SQR(ry0)+ SQR(rz0) );

						if (idc == -1)
						{
						//	fflush(stdout); exit(0);
							t_0 = -1;
							EscCond = 1;
							break;
						}

						px0 = rx0/r0;
						py0 = ry0/r0;
						pz0 = rz0/r0;

						vbulk0_x = vbulk_x;
						vbulk0_y = vbulk_y;
						vbulk0_z = vbulk_z;
						
						vbulk_x = px0 * CellArr[idc].vbulk[0];
						vbulk_y = py0 * CellArr[idc].vbulk[0];
						vbulk_z = pz0 * CellArr[idc].vbulk[0];
		
						xp += (ni * (vbulk0_x-vbulk_x) + nj*(vbulk0_y-vbulk_y) + nk*(vbulk0_z-vbulk_z) )/vth;

						if (strcmp(GeomName,"Biconical_Wind")==0 && idc != 0)
						{		
							th_c = acos(pz0);							
							EscCond = (pz0 > 0) ? ( (th_c > app_angle) ? 1 : 0) : \
									  ((th_c < Pi - app_angle) ? 1 : 0);

							if (EscCond == 1)
							{

/*								ObscuredCond = (pz0>0) ? ( (th0 > Pi/2.) ? 1 : 0) :\
												((th0 < Pi/2.) ? 1 : 0);
					
								if (ObscuredCond == 1)
								{
//									printf("Photon %d absorbed through perpendicular plane coming from %s side of cone:\ntheta=%f phi=%f r/RSphere=%f app=%f, nscat=%ld\n\n",\
									ip,(pz0 > 0) ? "+" : "-",th0*180./Pi,ph0*180./Pi,r0/RSphere,app_angle*180./Pi,nscat);
									idc = -1;
									goto absorbed;
								}

*/
								t_0 = -1;
//								printf("Photon %d escaped through %s side of cone:\ntheta=%f phi=%f r/RSphere=%f app=%f, nscat=%ld\n\n",\
								ip,(pz0 > 0) ? "+" : "-",th_c*180./Pi,atan2(py0,px0)*180./Pi,r0/RSphere,app_angle*180./Pi,nscat);
								idc = -1;
/*								find_edgecone(rx0,ry0,rz0,r0,th0,ph0,pos3d);
								rxf = pos3d[0];
								ryf = pos3d[1];
		 						rzf = pos3d[2];
								radius = sqrt(SQR(rxf) + SQR(ryf) + SQR(rzf));
//								printf("corrected position: %f,%f,%f\n",rx0,ry0,rz0); 
*/
								rxf = rx0;
								ryf = ry0;
								rzf = rz0;

								break;
							}

						}
/*						if (r0 < CellArr[idc].r || r0 > CellArr[idc].r+dr)
						{
							printf("ERROR: Something wrong happened, r0 is outside the limits of cell %d\n",idc);
							dprintd(CellArr[idc].r);		
							dprintd(CellArr[idc].r+dr);
*/
//							dprintd(r0);
/*							r0 = (r0 < CellArr[idc].r) ? CellArr[idc].r : CellArr[idc].r+dr;	
							rx0 = px0 * r0;
							ry0 = py0 * r0;
							rz0 = pz0 * r0;
//						}
*/
						break;
						
					case 0:	

						t_0 = -1;
						vbulk_x = 0;
						vbulk_y = 0;
						vbulk_z = 0;

						if (fabs(rzf) >= zSize/2.)
						{
							EscCond =  1;	
							idc		= -1;	
							break;
						}
						
						break;

					case 2:
						t_0 = -1;
						vbulk_x = 0;
						vbulk_y = 0;
						vbulk_z = 0;

						if (radius >= RSphere)
						{
							EscCond 	=  1;	
							idc		= -1;	
							break;
						}
						break;
				}
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
	
					x = xp + (ni*vbulk_x + nj*vbulk_y + nk*vbulk_z)/vth;

//				x = c*(g*nup - nu0)/(vth*nu0) - g*(nup/nu0)*(ni*vbulk_x + nj*vbulk_y + nk*vbulk_z)/vth;

//				dprintd(x);
//				printf("Photon has escaped in %ld scatterings\n",nscat);
// ................................

					inter = 4;
					(void) time(&t2);
					op_time = (float) (t2-t1)/60.;
//				printf("Total calculation took %f minutes.\n",op_time);
					if (strcmp(OutMode,"Long")==0)
					{
						printf("\n");
						record_data_long(nscat,ip,x,rxf,ryf,rzf,radius,upar,uper1,uper2,H_x,xp,inter);
					}
					else
					{
//						record_data_short(nscat,ip,x,rxf,ryf,rzf,radius,H_x,inter,op_time,flag_zero);
						XArr[ip] = x;
						X0Arr[ip]=xp0;
						InterArr[ip] = inter;
						NscatArr[ip] = nscat;
						
#ifdef GETPOSDIR
						PosArr[3*ip] 	= rxf;
						PosArr[3*ip+1]	= ryf;
						PosArr[3*ip+2]	= rzf;
						AngleArr[2*ip]	= th0;
						AngleArr[2*ip+1]= ph0;	
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

				rx0 = rxf;
				ry0 = ryf;
				rz0 = rzf;

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

					xi3 = gsl_rng_uniform (r);
					xi3 = (xi3 == 0.0) ? gsl_rng_uniform (r) : xi3;
					xi4 = gsl_rng_uniform (r);
					xi5 = gsl_rng_uniform (r);
					xi6 = gsl_rng_uniform (r);
					xi7 = gsl_rng_uniform (r);


#ifdef TEST_RII
	xp = x_test;
	xabs = fabs(xp);	
#endif

#ifdef USEREJECTION
					upar = -999.0;
					i = 0;
					do
					{
						xi0 = gsl_rng_uniform (r);
						xi1 = gsl_rng_uniform (r);
						xi2 = gsl_rng_uniform (r);
						upar = vp_rejection(xp,a_par,xi0,xi1,xi2);
						i++;
					} while(upar == -999.0);
//					if (i > 100000)
//						printf("VP_REJECTION TOOK %d tries to get upar = %f, x = %f,  nscatter = %d\n" \
						 ,i,upar,xp,nscat);
#endif

					cosph0 = cos(ph0);
					sinph0 = sin(ph0);
					costh0 = cos(th0);
					sinth0 = sin(th0);

					Arg1_uper = sqrt(SQR(xcrit) - log(xi3));	
					Arg2_uper = 2*Pi*xi4;
					
//					uper1 = sqrt(SQR(xcrit) -log(xi3))*cos(2*Pi*xi4);	
//					uper2 = sqrt(SQR(xcrit) -log(xi3))*sin(2*Pi*xi4);	

					uper1 = Arg1_uper * cos(Arg2_uper);
					uper2 = Arg1_uper * sin(Arg2_uper);

					n1_i = sinph0;
					n1_j = -cosph0;
					n1_k = 0;
		
					n2_i = cosph0 * costh0 ;
					n2_j = sinph0 * costh0; 
					n2_k = -sinth0;	

//					v_i = upar*sin(th0)*cos(ph0) + uper1 * n1_i + uper2 *n2_i ;
//					v_j = upar*sin(th0)*sin(ph0) + uper1 * n1_j + uper2 *n2_j ;
//					v_k = upar*cos(th0) + uper1 * n1_k + uper2 *n2_k;

					v_i = upar*sinth0*cosph0 + uper1 * n1_i + uper2 *n2_i ;
					v_j = upar*sinth0*sinph0 + uper1 * n1_j + uper2 *n2_j ;
					v_k = upar*costh0 + uper1 * n1_k + uper2 *n2_k;

					utot = sqrt(SQR(v_i) + SQR(v_j) + SQR(v_k));
					if (utot == 0.)
					{
						dprintd(utot);
						exit(0);
					}

					alpha = acos(upar/utot);

					cosalpha 		= upar/utot;
					sinalpha 		= sqrt(1 - SQR(cosalpha));
					Inv_sinalpha 	= 1./sinalpha; 

					v_i = v_i/utot;
					v_j = v_j/utot;
					v_k = v_k/utot;

//					v_pi = (1./sin(alpha))* (sin(th0)*cos(ph0) - cos(alpha)*v_i);					
//					v_pj = (1./sin(alpha))* (sin(th0)*sin(ph0) - cos(alpha)*v_j);					
//					v_pk = (1./sin(alpha))* (cos(th0) - cos(alpha)*v_k);					

					v_pi = Inv_sinalpha* (sinth0*cosph0 - cosalpha*v_i);					
					v_pj = Inv_sinalpha* (sinth0*sinph0 - cosalpha*v_j);					
					v_pk = Inv_sinalpha* (costh0 - cosalpha*v_k);					

					vpn = sqrt(	SQR(v_pi) + SQR(v_pj) + SQR(v_pk));
					v_pi = v_pi/vpn;
					v_pj = v_pj/vpn;
					v_pk = v_pk/vpn;

					g = 1./(sqrt(1 - SQR(utot*vth/c)));

					k_p_par = (c*cosalpha - utot*vth)/(1 - (upar*vth*Inv_c));
					k_p_per = (c*sinalpha)/((1 - (upar*vth*Inv_c)));


					i = 0;
					while (xi5 >= DipList[i] && i < ndip)
					{
						xilow = DipList[i];
						mulow = DipList[i + NTAB1];
						i++;
					}
					xihigh = DipList[i];
					muhigh = DipList[i+NTAB1];
					cosmu = (muhigh-mulow)/(xihigh-xilow) * (xi5 - xihigh) + muhigh;

					if (i == ndip)
						cosmu = mulow;

/*
					if (nscat <100)
					{
						dprintd(xi5);
						dprintd(cosmu);
					}
					else
						exit(0);	
*/
//					cosmu = 2*xi5 - 1;

					mu =  acos(cosmu);
					sinmu = sin(mu);

					mu2 = 2*Pi*xi6;
					cosmu2 = cos(mu2);
					sinmu2 = sin(mu2);

					ko1 = Inv_c*(cosmu * k_p_par + sinmu*(k_p_per)*cosmu2);
					ko2 = Inv_c*(cosmu * k_p_per - sinmu*(k_p_par)*cosmu2);
					ko3 = sinmu * sinmu2;

					nko0 = sqrt( SQR(ko1) + SQR(ko2) + SQR(ko3));

					ko1 = ko1/nko0;
					ko2 = ko2/nko0;
					ko3 = ko3/nko0;

					th_ = acos(ko1);
								
					costh_ = cos(th_);
					sinth_ = sin(th_);
	
					k_po_par = ko1;
					k_po_per = sqrt(SQR(ko2) + SQR(ko3));

					k_o_par = (c*costh_ + utot*vth)/(1 + (utot*costh_*vth)*Inv_c);
					k_o_per = (c*sinth_)/((1 + (utot*costh_*vth)*Inv_c));

					Inv_kpoper = 1./k_po_per;

					v_po_i = Inv_kpoper * (ko2*v_pi + ko3*(v_pk*v_j - v_k*v_pj));
					v_po_j = Inv_kpoper * (ko2*v_pj - ko3*(v_pk*v_i - v_k*v_pi));
					v_po_k = Inv_kpoper * (ko2*v_pk + ko3*(v_pj*v_i - v_j*v_pi));

					ko_x = Inv_c * (k_o_par * v_i + k_o_per *v_po_i);	
					ko_y = Inv_c * (k_o_par * v_j + k_o_per *v_po_j);	
					ko_z = Inv_c * (k_o_par * v_k + k_o_per *v_po_k);	
				
					n_ko = sqrt( SQR(ko_x) + SQR(ko_y) + SQR(ko_z));
	
					ko_x = ko_x/n_ko;
					ko_y = ko_y/n_ko;
					ko_z = ko_z/n_ko;

					v_ko = utot * (v_i*ko_x + v_j*ko_y + v_k*ko_z);

					th0 = acos(ko_z);
					ph0 = atan2(ko_y,ko_x);
					ph0 = (ph0 < 0) ? (ph0 + 2*Pi) : ph0;

						
					nup = xp * dnd + nu0;
					xp += (v_ko - upar); 

					xpabs = fabs(xp);
					xcrit = (xpabs < xcritval) ? xcritval : 0.;


/*	
					if (xpabs > 5000 || isnan(x) || n_ko > 1.01 || n_ko < 0.99)
					{
						printf("********************************************\n");
						printf("ERROR: Frequency out of range, NAN reached or norm of ko != 1.0\n");
						printf("Number of scatterings so far: %ld\n",nscat);
						printf("Photon's location %f %f %f\n",(rxf),(ryf),(rzf));
						exit(0);
						idc = -1;
					}
*/

#ifdef TEST_RII
	record_data_long(nscat,ip,xp,rxf,ryf,rzf,radius,upar,uper1,uper2,H_x,\
	xp,inter);
	if (nscat == nscatmax)
	{
		dprintl(nscat);
		dprintl(nscatmax);
		dprints(OutMode);
		printf("Test of redistribution function for x = %f is done.\n",x_test);
		exit(0);
	}

#endif

				break;
				
				case 2:
// Interaction with dust		

					xi1 = gsl_rng_uniform (r);
					xi2 = gsl_rng_uniform (r);
					xi3 = gsl_rng_uniform (r);

					interd = (xi1 <= Albedo) ? 1 : 0;
					if (interd == 1)	//scattering
					{

						i = 0;
						while (xi2 >= HGList[i] && i < nHG)
						{
							xilow = HGList[i];
							hglow = HGList[i + NTAB2];
							i++;
						}

						xihigh = HGList[i];
						hghigh = HGList[i+NTAB2];

						th0 = (hghigh-hglow)/(xihigh-xilow) * (xi2- xihigh) + hghigh;
						if (i == nHG)
							th0 = hglow;

        				ph0 = 2*Pi*xi3;
					}
					else 
					{
						absorbed:
	
						x = xp + (ni*vbulk_x + nj*vbulk_y + nk*vbulk_z)/vth;
//						printf("Photon has been absorbed by dust after %ld scatters.\n",nscat);
						inter = 3;
						(void) time(&t2);
						op_time = (float) (t2-t1)/60.;
//						printf("Total calculation took %f minutes.\n",op_time);

						if (strcmp(OutMode,"Long")==0)
						{
							record_data_long(nscat,ip,x,rxf,ryf,rzf,radius,upar,uper1,uper2,H_x,\
							xp,inter);
							dprintl(nscat);
							printf("\n");
						}
						else
//							record_data_short(nscat,ip,x,rxf,ryf,rzf,radius,H_x,inter,flag_zero,op_time);
//							record_data_short(nscat,ip,x,ph0,th0,radius,inter,op_time,flag_zero);

						XArr[ip] = x;
						X0Arr[ip]=xp0;
						InterArr[ip] = inter;
						NscatArr[ip] = nscat;
#ifdef GETPOSDIR
						PosArr[3*ip] 	= rxf;
						PosArr[3*ip+1]	= ryf;
						PosArr[3*ip+2]	= rzf;
						AngleArr[2*ip]	= th0;
						AngleArr[2*ip+1]= ph0;	

//						record_data_pos(nscat,ip,x,ph0,th0,rxf,ryf,rzf,inter,op_time,flag_zero);
//#else
//						record_data_short(nscat,ip,x,ph0,th0,radius,inter,op_time,flag_zero);
#endif
/*
#ifdef GETPOSDIR
						record_data_pos(nscat,ip,x,ph0,th0,rxf,ryf,rzf,inter,op_time,flag_zero);
#else
						record_data_short(nscat,ip,x,ph0,th0,radius,inter,op_time,flag_zero);
#endif
*/	
							idc = -1;
						goto end;
					
					}
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
		if (ip < 100)
		{
			x = xp + (ni*vbulk_x + nj*vbulk_y + nk*vbulk_z)/vth;
			record_data_long(nscat,ip,x,rxf,ryf,rzf,radius,upar,uper1,uper2,H_x,\
			xp,inter);
		}
#endif
		if (gtype ==1)
		{
			if (radius <= R_inner)
			{
				printf("scatter %d of photon %d occurs within empty zone. Something is not ok\n",nscat,ip);
				printf("rx/r_inner %f ry/r_inner %f rz/r_inner %f radius/r_inner %f\nidc %d\n",rxf/R_inner,\
				ryf/R_inner,rzf/R_inner,radius/R_inner,idc);
			}
		}
		nscat++;
		}				


	end:
//	printf("done\n");


	

	if (strcmp(Set_Tolerance,"yes") == 0)
	{
//		if ((ip >= np_min && nout >= nout_max) || \
//			(ip >= np_max))


		if ((ip >= np_min && nout >= nout_max) || \
			(ip >= np_max) )
//			(ip >= np_max && nout >= (int) nout_max/5) || \

		{
			printf("Max number of absorbed photons or limit on the number of photons reached,  exiting...\n");

#ifdef GETPOSDIR
	printf("Writing data\n");
	record_data_pos(ip);
	printf("File %s written\n",OutShort);
	free(PosArr);
	free(AngleArr);
#else
	record_data_short();
#endif

			free(HList);
			free(DipList);
			free(HGList);
			free(CellArr);
			exit(0);
		}
	}		

#ifdef TIMELIMIT
	(void) time(&t4);
	t4 = (float) (t4-t0)/60.;
	if (t4 > MAXTIME)
	{
		printf("MAXTIME reached, forcing output and exit\n");
#ifdef GETPOSDIR
	record_data_pos(ip);
	printf("File %s written\n",OutShort);
	free(PosArr);
	free(AngleArr);
#else
	record_data_short();
#endif
		exit(0);
	}
#endif

	}

	printf("Writing data \n");

#ifdef GETPOSDIR
	record_data_pos(ip);
	printf("File %s written\n",OutShort);
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
