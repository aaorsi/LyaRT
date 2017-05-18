#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include "allvars.h"
#include "proto.h"

void dust_interaction(photon *P, int ip)
{
	double xilow,xihigh;
	long i;
    double cosph0 , sinph0 , costh0 , sinth0 ; 
    double Arg1_uper , Arg2_uper , uper1 , uper2 ;
    double n1_i , n1_j , n1_k , n2_i , n2_j , n2_k , v_i , v_j , v_k , utot ; 
    double alpha , cosalpha , sinalpha , Inv_sinalpha ;
    double v_pi , v_pj , v_pk , vpn ; 
    float g;
    double k_p_par , k_p_per;
    double mulow , muhigh , cosmu , sinmu , mu , cosmu2 , sinmu2 , mu2;
    double ko1 , ko2 , ko3 , nko0 ;
    double th_ , costh_ , sinth_ ; 
    double k_po_par , k_po_per , k_o_par , k_o_per , Inv_kpoper , v_po_i , v_po_j , v_po_k , ko_x , ko_y , ko_z , n_ko , v_ko; 


	interd = (xi1 <= Albedo) ? 1 : 0;
	if (interd == 1)	//scattering
	{

        // creamos estas viariables para no repetir el mismo calculo muchas veces
        cosph0 = cos(ph0);
        sinph0 = sin(ph0);
        costh0 = cos(th0);
        sinth0 = sin(th0);

        // Estos son los angulos de la velocidad del grano de polvo ( eq. 5.37 y eq. 5.38 ) 
        Arg1_uper = sqrt(SQR(xcrit) - log(xi3));       
        Arg2_uper = 2*Pi*xi4;

        uper1 = Arg1_uper * cos(Arg2_uper); // esta es la velocidad perpendicular 1 (eq. 5.37)
        uper2 = Arg1_uper * sin(Arg2_uper); // esta es la velocidad perpendicular 2 (eq. 5.38)

        // Estos son unos vectores muy convenientes para algo.
        // parte paralela ??
        n1_i = sinph0;
        n1_j = -cosph0;
        n1_k = 0;
    
        // parte transversal ??
        n2_i = cosph0 * costh0 ;
        n2_j = sinph0 * costh0; 
        n2_k = -sinth0; 

        // creo que esto es la velocidad del photon en el sistema de referencia del laboratio
        v_i = upar*sinth0*cosph0 + uper1 * n1_i + uper2 *n2_i ;
        v_j = upar*sinth0*sinph0 + uper1 * n1_j + uper2 *n2_j ;
        v_k = upar*costh0 + uper1 * n1_k + uper2 *n2_k;

        // Tal vez esto se deberia llamar vtot, pero que le vamos a hacer. Es el modulo de v.
        utot = sqrt(SQR(v_i) + SQR(v_j) + SQR(v_k));
        if (utot == 0.)
        {
            dprintd(utot);
            exit(0);
        }

        // Este alpha es el angulo entre la componente paralela del movimiento del grano de polvo la vel del photon. (eq. 5.54)
        alpha           = acos(upar/utot);
        cosalpha        = upar/utot;
        sinalpha        = sqrt(1 - SQR(cosalpha));
        Inv_sinalpha    = 1./sinalpha;

        // Aqui se normaliza v ... se podria hacer antes ...
        v_i = v_i/utot;
        v_j = v_j/utot;
        v_k = v_k/utot;

        // Sera esto algo asi como la velocidad proyectada???
        v_pi = Inv_sinalpha* (sinth0*cosph0 - cosalpha*v_i);
        v_pj = Inv_sinalpha* (sinth0*sinph0 - cosalpha*v_j);
        v_pk = Inv_sinalpha* (costh0 - cosalpha*v_k);

        // Este es el modulo de la velocidad esa, que vete tu a saber que es, y se normaliza
        vpn = sqrt( SQR(v_pi) + SQR(v_pj) + SQR(v_pk));
        v_pi = v_pi/vpn;
        v_pj = v_pj/vpn;
        v_pk = v_pk/vpn;
        
        //Esto si que no tengo ni idea de que es, pero no le hace daño a nadie... o eso espero :(
        g = 1./(sqrt(1 - SQR(utot*vth/c)));

        //Estas cosas son importantes, todas las lineas de antes es para poder calcularlo.
        k_p_par = (c*cosalpha - utot*vth)/(1 - (upar*vth*Inv_c));
        k_p_per = (c*sinalpha)/((1 - (upar*vth*Inv_c)));

        // Ahora vamos a obtener el coseno del angulo en el que sale disparado el photon en el frame del grano de polvo
		i = 0;
		while (xi2 >= HGList[i] && i < nHG)
		{
			xilow = HGList[i];
			mulow = HGList[i + NTAB2];
			i++;
		}
		xihigh = HGList[i];
		muhigh = HGList[i+NTAB2];

		cosmu = (muhigh-mulow)/(xihigh-xilow) * (xi2- xihigh) + muhigh;

		if (i == nHG)
            cosmu = mulow;

        // Este es el angulo polar en el sistema de referencia del atomo de hidrogeno
        mu =  acos(cosmu); 
        sinmu = sin(mu);

        // Este es el angulo azimutal en el sistema de referencia del atomo de hidrogeno
        mu2 = 2*Pi*xi6; 
        cosmu2 = cos(mu2);
        sinmu2 = sin(mu2);

        // Esto me supera, la verdad. Es posible que sea la direccion del photon en el sistema del grano de polvo.
        ko1 = Inv_c*(cosmu * k_p_par + sinmu*(k_p_per)*cosmu2);
        ko2 = Inv_c*(cosmu * k_p_per - sinmu*(k_p_par)*cosmu2);
        ko3 = sinmu * sinmu2;

        // Se calcula el modulo y se normaliza
        nko0 = sqrt( SQR(ko1) + SQR(ko2) + SQR(ko3));
        ko1 = ko1/nko0;
        ko2 = ko2/nko0;
        ko3 = ko3/nko0;

        // Este angulo hace llorar al ninio jesus...
        th_ = acos(ko1);
        costh_ = cos(th_);
        sinth_ = sin(th_);

        // Estoy empezando a pensar que tal vez la 'p' en el nombre se refiere a un sistema de referencia
        k_po_par = ko1;
        k_po_per = sqrt(SQR(ko2) + SQR(ko3));

        // Aqui hace magia...
        k_o_par = (c*costh_ + utot*vth)/(1 + (utot*costh_*vth)*Inv_c);
        k_o_per = (c*sinth_)/((1 + (utot*costh_*vth)*Inv_c));

        // Alvaro es un mago
        Inv_kpoper = 1./k_po_per;

        v_po_i = Inv_kpoper * (ko2*v_pi + ko3*(v_pk*v_j - v_k*v_pj));
        v_po_j = Inv_kpoper * (ko2*v_pj - ko3*(v_pk*v_i - v_k*v_pi));
        v_po_k = Inv_kpoper * (ko2*v_pk + ko3*(v_pj*v_i - v_j*v_pi));

        // Estas son las coordenadas finales de la V del photon en el Lab Sys
        ko_x = Inv_c * (k_o_par * v_i + k_o_per *v_po_i); 
        ko_y = Inv_c * (k_o_par * v_j + k_o_per *v_po_j); 
        ko_z = Inv_c * (k_o_par * v_k + k_o_per *v_po_k);

        // Calcula el modulo y normaliza
        n_ko = sqrt( SQR(ko_x) + SQR(ko_y) + SQR(ko_z));
        ko_x = ko_x/n_ko;
        ko_y = ko_y/n_ko;
        ko_z = ko_z/n_ko;

        // Pruedo prometer y prometo que no se que es esto, pero parace ser para el cambio en frecuencia, que aquí no se aplica.
        v_ko = utot * (v_i*ko_x + v_j*ko_y + v_k*ko_z);

        // Finalmente obtenemos los angulos en el sistema de referencia del laboratorio.
        th0 = acos(ko_z);
        ph0 = atan2(ko_y,ko_x);
        ph0 = (ph0 < 0) ? (ph0 + 2*Pi) : ph0; // si ph0 es menor que 0 le suma 2Pi???

        // Le asignamos los angulos que hemos sacado del cajon de los angulos buenos al photon.
		P[ip].th = th0;
		P[ip].ph = ph0;
		end_syg = 0;

        // Un aplauso para Alvaro por favor.

	}
	else 
	{
		absorbed:
		P[ip].xp += (ni*vbulk_x + nj*vbulk_y + nk*vbulk_z)/vth;
//			printf("Photon has been absorbed by dust after %ld scatters.\n",nscat);
		inter = 3;
//		(void) time(&t2);
//		op_time = (float) (t2-t1)/60.;
//				printf("Total calculation took %f minutes.\n",op_time);
		if (strcmp(OutMode,"Long")==0)
		{
			record_data_long(nscat,ip,x,rxf,ryf,rzf,radius,upar,uper1,uper2,H_x,\
			xp,inter);
			dprintl(nscat);
			printf("\n");
		}
								
//							record_data_short(nscat,ip,x,rxf,ryf,rzf,radius,H_x,inter,flag_zero,op_time);
//							record_data_short(nscat,ip,x,ph0,th0,radius,inter,op_time,flag_zero);

		XArr[ip] = P[ip].xp;
		X0Arr[ip]=xp0;
		InterArr[ip] = inter;
		NscatArr[ip] = nscat;
#ifdef GETPOSDIR
		PosArr[3*ip] 	= P[ip].x;//rxf;
		PosArr[3*ip+1]	= P[ip].y;//ryf;
		PosArr[3*ip+2]	= P[ip].z;//rzf;
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
		P[ip].idc = idc;
		end_syg = 1;
	}
}						
