#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

void scattering_hydrogen(photon *P, int ip)
{

	const gsl_rng_type * T;
	double k_p_par,k_p_per,mu,mu2,ko1,ko2,ko3,k_po_par,k_po_per,k_o_par,k_o_per,b_ang;
	double n1_i,n1_j,n1_k,n2_i,n2_j,n2_k,v_i,v_j,v_k,v_pi,v_pj,v_pk,v_po_i,v_po_j,v_po_k;
	double xilow,xihigh,mulow,muhigh,hglow,hghigh;
	double ko_x,ko_y,ko_z;
	double x0,y0,z0,th_;
	double xcw,th_aux,th_range[2],u0,u_par,g_u,accept,u_p1,u_p2;
	double a1,a2,a3,vpar_i,vpar_j,vpar_k,vper_i,vper_j,vper_k,v_ko,uper,vper;
	float g;
	long i;

	gsl_rng *r;
	gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
	
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
//	if (i > 100000)
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

	P[ip].th = th0;
	P[ip].ph = ph0;
	P[ip].xp = xp;

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
	record_data_long(nscat,ip,xp,rxf,ryf,rzf,r0,upar,uper1,uper2,H_x,\
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
}
