/* Here different geometries are defined.
The grid of cells is filled according to given parameters
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

void define_geometry(char *GeomName, cell *CellArr)
{

	long i,j,k,id;
	FILE *fout;
	int n,GType;
	double r;
	double v_r;	
	double A_par;
	double x_nHII = 0;
	double x_nHII_shield,x_nHII_thin;	

	double Gamma_z_thin,beta_T,y,Gamma_z_shield,Gamma_z_from_spec,Gamma_z;
	double tau_0_nH,sig0,nh_r; 

#ifdef WRITEPROFILE
		
		printf("Writing density and velocity profiles \n");
		char nHProfile[NCHARMAX];
//		strcpy(nHProfile,"../data/test_nh/test_UV_z3.0");
		sprintf(nHProfile,"../data/test_nh/test_UV%d_z%g",IncShield,redshift);
		printf("%s\n",nHProfile);
		fout = fopen(nHProfile,"w");
		fprintf(fout,"# r	v[km s-1]	nH_0	x_nHII(shield) x_nHII(thin) \n");
#endif

	spec *minispec;
	float *uv_z, *uv_logenergy;

	if (IncUV == 1)
	{
		printf("Including ionization from the UV background\n");


		Gamma_z_thin = Gamma_HI(redshift);	// photoionization rate without self-shielding

		printf("Allocating memory for UV data\n");
		minispec 		= (spec*) malloc(NEN * sizeof(spec));
		uv_z	 		= (float*) malloc(NZZ * sizeof(float));
		uv_logenergy	= (float*) malloc(NEN * sizeof(float));
	 
		set_spec(redshift, uv_z, uv_logenergy, minispec);	
	}
		

	if (strcmp(GeomName,"HomSlab") ==0)			// Homogeneous, static slab
		GType = 0;
	if (strcmp(GeomName,"HomSphere") ==0)		// Homogeneous, linearly expanding sphere
		GType = 1;
	if (strcmp(GeomName,"Exp_Const_Sphere") ==0) // expanding sphere with constant velocity
		GType = 2;
	if (strcmp(GeomName,"Wind") ==0)	  // A wind expanding with constant velocity 
		GType = 3;
	if (strcmp(GeomName,"ThinShell") ==0)		  // A thin shell expanding with constant velocity
		GType = 4;
	if (strcmp(GeomName,"Wind2") ==0)			// A wind plus a linearly expanding sphere
		GType = 5;
	if (strcmp(GeomName,"Wind3") ==0)			// A wind plus a static sphere 
		GType = 6;
	if (strcmp(GeomName,"Wind4") ==0 )			// A wind plus a static sphere with a gap in between
		GType = 7;
	if (strcmp(GeomName,"ThinShell2") ==0)		// A shell plus a static sphere
		GType = 8;
	if (strcmp(GeomName,"Wind_DijkstraProfile") ==0)
		GType = 9;
	if (strcmp(GeomName,"Wind_VelProfile") ==0)
	  GType = 10;
	if (strcmp(GeomName,"Biconical_Wind") == 0)   // Similar to a wind, but with a biconical geometry.
	  GType = 11;

	id=0;


	switch(GType)
	{
		case 0:

			for(i=0;i<NCellX;i++)
			{
				for(j=0;j<NCellY;j++)
				{
					for(k=0;k<NCellZ;k++)
					{
						id = i + j*NCellX + k*NCellX*NCellY;

						CellArr[id].id 	 = i;
//						CellArr[id].neig = {id-1,id+1,id-NCellX,id+NCellX,id-NCellX*NCellY,id+NCellX*NCellY};
						CellArr[id].neig[0] = id-1;
						CellArr[id].neig[1] = id+1;
						CellArr[id].neig[2] = id-NCellX;
						CellArr[id].neig[3] = id+NCellX;
						CellArr[id].neig[4] = id- NCellX*NCellY;
						CellArr[id].neig[5] = id+ NCellX*NCellY;


//	Pointing to escape	
						if (i == 0) 
							CellArr[id].neig[0] = -1;
						if (i == NCellX-1)	
							CellArr[id].neig[1] = -1;
						if (j == 0) 
							CellArr[id].neig[2] = -1;
						if (j == NCellY-1)	
							CellArr[id].neig[3] = -1;
						if (k == 0) 
							CellArr[id].neig[4] = -1;
						if (k == NCellZ-1)	
							CellArr[id].neig[5] = -1;
					
						if (ZPeriodic)
						{
							if(k == 0)
							CellArr[id].neig[4] = id + NCellX*NCellY*(NCellZ-1);
							if(k == NCellZ-1)
								CellArr[id].neig[5] = id - NCellX*NCellY*(NCellZ-1);
						}
						if (YPeriodic)
						{

							if(j == 0)
							CellArr[id].neig[2] = id + NCellX*(NCellY-1);
							if(j == NCellY-1)
							CellArr[id].neig[3] = id - NCellX*(NCellY-1);
						}
						if (XPeriodic)
						{
						if(i == 0)
							CellArr[id].neig[0] = id + (NCellX-1);
						if(i == NCellX-1)
							CellArr[id].neig[1] = id - (NCellX-1);
						}

						CellArr[id].id 	     = i;
						CellArr[id].p_lya    = 1.;
						CellArr[id].T	     = pow(10,Temp);
						CellArr[id].nH	     = mean_nH;
						CellArr[id].z	     = mean_z;
						CellArr[id].vbulk[0] = 0.;
						CellArr[id].vbulk[1] = 0.;
						CellArr[id].vbulk[2] = 0.;

					}
				}
			}

		break;		

//	Grid with shells	
		case 1:
			
		H = vmax / RSphere;
		dr = (NCells > 1) ? RSphere / (NCells-1) : 0.0;

		for (i = 0; i<NCells; i++)
		{
			r = i * dr; 

			id = i;
			CellArr[id].id =  i;			
			CellArr[id].neig[0] = (id == 0 )? 0 : id-1;
			CellArr[id].neig[1] = id+1;
			if (i == NCells-1)
				CellArr[id].neig[1] = -1;
			
			CellArr[id].p_lya 	= 1.;
			CellArr[id].T		= pow(10,Temp);
			CellArr[id].nH		= mean_nH;
			CellArr[id].z		= mean_z;
			CellArr[id].vbulk[0] = H * r;
			CellArr[id].r		 = r;			
		}
		break;	

		case 2:
			
		dr = RSphere / NCells;
		H = vmax;

		for (i = 0; i<NCells; i++)
		{
			r = i * dr; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] = id-1;
			CellArr[id].neig[1] = id+1;
			if (i == NCells-1)
				CellArr[id].neig[1] = -1;
			
			CellArr[id].p_lya 	= 1.;
			CellArr[id].T		= pow(10,Temp);
			CellArr[id].nH		= mean_nH;
			CellArr[id].z		= mean_z;
			CellArr[id].vbulk[0] = (i == 0) ? 0 : vmax;
			CellArr[id].r		 = r;			
		}
		break;	
		
		
		
		case 3:

/* 	A shell with density decreasing as r^-2 and constant outflow velocity
	Here RSphere corresponds to the outer radius of the shell. The inner radius is given by R_inner
	The number density at a given radius r inside the shell is given by
		nH(r) = K/(v * r^2),
	where v[km s-1] and r [cm]. The constant K is an input of the code (for flexibility), and 
	in general it depends on the mass ejection rate Mdot and other numerical constants to make the units work.
	In summary, this configuration requires to define:
		- An inner radius R_inner, from which the number density is greater than zero.
		- An outer radius RSphere (which should be, in general, > 20*R_inner)
		- The velocity of the outflow vmax
		- The constant KnH, relating nH(r) to v and r

* REMEMBER: The units must be: Ndot [1e-40 s-1], r[cm] and v[kms-1]. If the units change
 then change the constant knH accordingly.

*/

		id = 0;
		CellArr[id].id = 0;			
		CellArr[id].neig[0]  = 0;
		CellArr[id].neig[1]  = id+1;
		CellArr[id].p_lya 	 = 1.;
		CellArr[id].T		 = pow(10,Temp);
		CellArr[id].z		 = 0;
		CellArr[id].r		 = 0;	
		CellArr[id].nH		 = 0;
		CellArr[id].vbulk[0] = 0;
	
		dr = (RSphere-R_inner) / (NCells-2);


//		H = 0.;
//		H = vmax / RSphere;
		for (i = 1; i<NCells; i++)
		{
			r = (i-1) * dr + R_inner; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] = id-1;
			CellArr[id].neig[1] = id+1;
			if (i == NCells-1)
			CellArr[id].neig[1] = -1;
			
			CellArr[id].p_lya 	= 1.;
			CellArr[id].T		= pow(10,Temp);
			CellArr[id].z		= mean_z;
			CellArr[id].r		= r;			
			CellArr[id].vbulk[0]= vmax;
			CellArr[id].nH		= KnH * Ndot / ( 4 * Pi * SQR(r) * vmax);
			
		}	

		break;	

//	Shell with constant velocity, constant density.	

		case 4:
		printf("Defining Configuration 4:%s\n",GeomName);
	
		id = 0;
		CellArr[id].id = 0;			
		CellArr[id].neig[0]  = 0;
		CellArr[id].neig[1]  = id+1;
		CellArr[id].p_lya 	 = 1.;
		CellArr[id].T		 = pow(10,Temp);
		CellArr[id].z		 = 0;
		CellArr[id].r		 = 0;			
		CellArr[id].nH		 = 0;
		CellArr[id].vbulk[0] = 0;
	
		dr = (RSphere-R_inner) / (NCells-2);
		for (i = 1; i<NCells; i++)
		{
			r = (i-1) * dr + R_inner; 
			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] 	= id-1;
			CellArr[id].neig[1] 	= id+1;
			if (i == NCells-1)
				CellArr[id].neig[1] = -1;
			CellArr[id].p_lya 		= 1.;
			CellArr[id].T			= pow(10,Temp);
			CellArr[id].nH			= mean_nH;
			CellArr[id].z			= mean_z;
			CellArr[id].vbulk[0] 	= vmax;
			CellArr[id].r		 	= r;			
		}
		break;	

		case 5:

//	A wind consisting of a homogeneous sphere expanding linearly outwards (like case 1) plus
//  a wind like the one described in case 3.
	

//		dr = (RSphere-R_inner) / (NCells-1);

		H 			= vmax / R_inner;
//		NCellSphere	= R_inner/dr + 1.0;

		for (i = 0; i<NCellSphere; i++)
		{
			r = i * dr; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0]  = (id == 0 )? 0 : id-1;
			CellArr[id].neig[1]  = id+1;
			CellArr[id].p_lya 	 = 1.;
			CellArr[id].T		 = pow(10,Temp);
			CellArr[id].nH		 = KnH * Ndot / (4 * Pi * SQR(R_inner) * vmax);
			CellArr[id].z		 = mean_z;
			CellArr[id].vbulk[0] = H * r;
			CellArr[id].r		 = r;			
		}

		for (i = NCellSphere ; i<NCells ; i++)
		{
			r = i * dr; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] = id-1;
			CellArr[id].neig[1] = id+1;
			if (i >= NCells-NCellSphere-1)
			CellArr[id].neig[1] = -1;
			
			CellArr[id].p_lya 	= 1.;
			CellArr[id].T		= pow(10,Temp);
			CellArr[id].z		= mean_z;
			CellArr[id].r		= r;			
			CellArr[id].vbulk[0]= vmax;
			CellArr[id].nH		= KnH * Ndot / ( 4 * Pi * SQR(r) * vmax);
			
		}	



//		NCells += NCellSphere;
		dprinti(NCells);
		dprinti(NCellSphere);
		dprintd(dr);
		dprintd(R_inner);
		break;	

	case 6:

//	A wind consisting of a static homogeneous sphere plus
//  a wind like the one described in case 3.
	

//		dr = (RSphere-R_inner) / (NCells-1);

		H 			= vmax / R_inner;
//		NCellSphere	= R_inner/dr + 1.0;

		for (i = 0; i<NCellSphere; i++)
		{
			r = i * dr; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0]  = (id == 0 )? 0 : id-1;
			CellArr[id].neig[1]  = id+1;
			CellArr[id].p_lya 	 = 1.;
			CellArr[id].T		 = pow(10,Temp);
			CellArr[id].nH		 = KnH * Ndot / (4 * Pi * SQR(R_inner) * vmax);
			CellArr[id].z		 = mean_z;
			CellArr[id].vbulk[0] = 0;
			CellArr[id].r		 = r;			
		}

		for (i = NCellSphere ; i<NCells ; i++)
		{
			r = i * dr; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] = id-1;
			CellArr[id].neig[1] = id+1;
			if (i >= NCells-NCellSphere-1)
			CellArr[id].neig[1] = -1;
			
			CellArr[id].p_lya 	= 1.;
			CellArr[id].T		= pow(10,Temp);
			CellArr[id].z		= mean_z;
			CellArr[id].r		= r;			
			CellArr[id].vbulk[0]= vmax;
			CellArr[id].nH		= KnH * Ndot / ( 4 * Pi * SQR(r) * vmax);
			
		}	



//		NCells += NCellSphere;
		dprinti(NCells);
		dprinti(NCellSphere);
		dprintd(dr);
		dprintd(R_inner);
		break;

	case 7:

//	A wind consisting of a static homogeneous sphere reaching a certain radius R_static < R_inner, and
//  a wind like the one described in case 3. This resembles (although is not equal to) a static slab with a
//  given column density surounded by a wind.
  
		
//		dr = (RSphere-R_inner) / (NCells-1);

// 		Inner sphere:

		for (i = 0; i<NCellSphere; i++)
		{
			r = i * dr;

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0]  = (id == 0 )? 0 : id-1;
			CellArr[id].neig[1]  = id+1;
			CellArr[id].p_lya 	 = 1.;
			CellArr[id].T		 = pow(10,Temp);
			CellArr[id].nH		 = mean_nH;
//			dprintd(CellArr[id].nH);
			CellArr[id].z		 = mean_z;
			CellArr[id].vbulk[0] = 0;
			CellArr[id].r		 = r;			
		}

//		Empty space in between

		i = NCellSphere;
		r = i * dr;

		id = i;
		CellArr[id].id = i;			
		CellArr[id].neig[0]  = (id == 0 )? 0 : id-1;
		CellArr[id].neig[1]  = id+1;
		CellArr[id].p_lya 	 = 1.;
		CellArr[id].T		 = pow(10,Temp);
		CellArr[id].nH		 = 1e-18;
		CellArr[id].z		 = 0;
		CellArr[id].vbulk[0] = 0;
		CellArr[id].r		 = r;			

//		Wind
		for (i = NCellSphere+1 ; i < NCells ; i++)
		{
			r = (i-(NCellSphere+1)) * dr + R_inner; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] = id-1;
			CellArr[id].neig[1] = id+1;
			if (i >= NCells)
			CellArr[id].neig[1] = -1;
			
			CellArr[id].p_lya 	= 1.;
			CellArr[id].T		= pow(10,Temp);
			CellArr[id].z		= mean_z;
			CellArr[id].r		= r;
			CellArr[id].vbulk[0]= vmax;
			CellArr[id].nH		= KnH * Ndot / ( 4 * Pi * SQR(r) * vmax);
			
		}	

//		NCells += NCellSphere;
		dprinti(NCells);
		dprinti(NCellSphere);
		dprintd(dr);
		dprintd(R_inner);
		break;	

	case 8:

//	An outflow consisting of a homogeneous sphere expanding linearly outwards (like case 1) plus
//  a shell like the one described in case 3.
	

//		dr = (RSphere-R_inner) / (NCells-1);

		H 			= vmax / R_inner;
//		NCellSphere	= R_inner/dr + 1.0;

		for (i = 0; i<NCellSphere; i++)
		{
			r = i * dr1; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0]  = (id == 0 )? 0 : id-1;
			CellArr[id].neig[1]  = id+1;
			CellArr[id].p_lya 	 = 1.;
			CellArr[id].T		 = pow(10,Temp);
			CellArr[id].nH		 = mean_nH_static;
			CellArr[id].z		 = mean_z;
			CellArr[id].vbulk[0] = 0.;
			CellArr[id].r		 = r;			
		}


//		dr = (RSphere-R_inner) / (NCells-1);
		for (i = NCellSphere; i<NCells-NCellSphere; i++)
		{
			r = (i-1) * dr2 + R_inner; 
			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] 	= id-1;
			CellArr[id].neig[1] 	= id+1;
			if (i == NCells-1)
				CellArr[id].neig[1] = -1;
			CellArr[id].p_lya 		= 1.;
			CellArr[id].T			= pow(10,Temp);
			CellArr[id].nH			= mean_nH;
			CellArr[id].z			= mean_z;
			CellArr[id].vbulk[0] 	= vmax;
			CellArr[id].r		 	= r;			
		}
		break;	

	case 9:
	
		id = 0;
		CellArr[id].id = 0;			
		CellArr[id].neig[0]  = 0;
		CellArr[id].neig[1]  = id+1;
		CellArr[id].p_lya 	 = 1.;
		CellArr[id].T		 = pow(10,Temp);
		CellArr[id].z		 = 0;
		CellArr[id].r		 = 0;	
		CellArr[id].nH		 = 0;
		CellArr[id].vbulk[0] = 0;
	
		dr = (RSphere-R_inner) / (NCells-2);

		for (i = 1;i<NCells-1; i++)
		{
			CellArr[i].nH = 0;
			CellArr[i].vbulk[0] = 0;	
		}


//		H = 0.;
//		H = vmax / RSphere;
		A_par = pow(2.5*3.08e21,alpha_vprof-1);
//		A_par = pow(2.5*R_inner,alpha_vprof-1);
		dprintd(A_par);
		dprintd(R_inner);
		dprintd(sigma_vprof);
		dprintd(alpha_vprof);
		for (i = 1; i<NCells; i++)
		{
			r = (i-1) * dr + R_inner; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] = id-1;
			CellArr[id].neig[1] = id+1;
			if (i == NCells-1)
			CellArr[id].neig[1] = -1;
			CellArr[id].p_lya 	= 1.;
			CellArr[id].T		= pow(10,Temp);
			CellArr[id].z		= mean_z;
			CellArr[id].r		= r;			
			v_r = 2*sigma_vprof*sqrt(log(R_inner/r)+ (A_par/(alpha_vprof-1)) * \
				  (-pow(r,1-alpha_vprof) + pow(R_inner,1-alpha_vprof)));
			if (isnan(v_r)) 
			{
				v_r = 0;
				printf("Maximum radius %f\n",r);
				break;
			}
			CellArr[id].vbulk[0]= v_r;
			CellArr[id].nH		= (v_r == 0) ? 0: KnH * Ndot / ( 4 * Pi * SQR(r) * v_r);	
			if (i < 10)
			{
				dprinti(i);  
			    dprintd(v_r);
/*				dprintd(log(R_inner/r));
				dprintd(pow(r,1-alpha_vprof));
				dprintd(pow(R_inner,1-alpha_vprof));
				dprintd(log(R_inner/r)+ (A_par/(alpha_vprof-1)) * \
                (-pow(r,1-alpha_vprof) + pow(R_inner,1-alpha_vprof)));
*/
				dprintd(r);
				dprintd(CellArr[id].nH);
			}
		}			
	
		break;	
	

  //wind with velocity profile
	case 10:
	  
	  id = 0;
	  CellArr[id].id = 0;
	  CellArr[id].neig[0]  = 0;
	  CellArr[id].neig[1]  = id+1;
	  CellArr[id].p_lya        = 1.;
	  CellArr[id].T            = pow(10,Temp);
	  CellArr[id].z            = 0;
	  CellArr[id].r            = 0;
	  CellArr[id].nH           = 0;
	  CellArr[id].vbulk[0] = 0;
	  //CellArr[NCells].idc=-1;
	  dr = (RSphere-R_inner) / (NCells-2);


	  //              H = 0.;                                                           
	  //              H = vmax / RSphere;                                               
	  
	  dprintd(R_inner);
	  dprintd(R_Peak);
	  

	  for (i = 1; i<NCells; i++)
	    {
	      r = (i-1) * dr + R_inner;
	      
	      id = i;
	      CellArr[id].id = i;
	      CellArr[id].neig[0] = id-1;
	      CellArr[id].neig[1] = id+1;
	      if (i == NCells-1)
		CellArr[id].neig[1] = -1;
	      CellArr[id].p_lya       = 1.;
	      CellArr[id].T           = pow(10,Temp);
	      CellArr[id].z           = mean_z;
	      CellArr[id].r           = r;

	      
	      if( r < R_Peak)
		{
		  v_r = vmax*(r)/(R_Peak)+1.;
		  
		}
	      else
		{
		  v_r=vmax*(RSphere - r)/(RSphere - R_Peak);
		  
		}
		printf("i: %d r = %f v_r = %f\n",i,r,v_r);
	      
	      if (isnan(v_r))
		{
		  v_r = 0;
		  printf("Maximum radius %f\n",r);
		  //	  break;
		}
	      CellArr[id].vbulk[0]= v_r;
	      CellArr[id].nH= (v_r == 0) ? 0: KnH * Ndot / ( 4* Pi * SQR(r) * v_r);
	      if (i < 3)
		{
		  dprintd(v_r);
		  dprintd(r);
		  dprintd(CellArr[id].nH);
		}
	    }
	  

	  break;

	case 11:

		app_angle = Pi/4.;  // By default.

		id = 0;
		CellArr[id].id = 0;			
		CellArr[id].neig[0]  = 0;
		CellArr[id].neig[1]  = id+1;
		CellArr[id].p_lya 	 = 1.;
		CellArr[id].T		 = pow(10,Temp);
		CellArr[id].z		 = 0;
		CellArr[id].r		 = 0;	
		CellArr[id].nH		 = 0;
		CellArr[id].vbulk[0] = 0;
	
		dr = (RSphere-R_inner) / (NCells-2);
		printf("dr for Biconical wind:%e\n",dr);
//		H = 0.;
//		H = vmax / RSphere;

		for (i = 1; i<NCells; i++)
		{
			r = (i-1) * dr + R_inner; 

			id = i;
			CellArr[id].id = i;			
			CellArr[id].neig[0] = id-1;
			CellArr[id].neig[1] = id+1;
			if (i == NCells-1)
			CellArr[id].neig[1] = -1;
			CellArr[id].p_lya 	= 1.;
			CellArr[id].T		= pow(10,Temp);
			CellArr[id].z		= mean_z;
			CellArr[id].r		= r;			
			CellArr[id].vbulk[0]= vmax;
			CellArr[id].nH		= KnH * Ndot / ( -4 * Pi * (cos(app_angle)-1) * SQR(r) * vmax);	
		}

		break;

	}

	if (IncUV == 1)
/* 		    Compute x_nHII from the UV background */
	{
		tau_0_nH = 0;
		sig0 = HI_photo_cs_analytic(1.00);
		for (i = NCells-1;i>=0;i--)			
		{
			id = i;
			beta_T = beta_HI(CellArr[id].T * 1e4);
			nh_r = 0;

			for (j = NCells-1; j>=i;j--)
			{ 
				nh_r += CellArr[j].nH * dr/2.;
			}


			tau_0_nH = (IncShield == 1) ? sig0 * nh_r : 0.;	
			Gamma_z_shield = gamma_shield(minispec, tau_0_nH);
			Gamma_z = (IncShield == 1) ? Gamma_z_shield : Gamma_z_thin;

			x_nHII_shield = get_xnHII(Gamma_z_shield,beta_T, id, CellArr);
			x_nHII_thin   = get_xnHII(Gamma_z_thin,beta_T, id, CellArr );

//			printf("%d\n",i);
#ifdef WRITEPROFILE
			fprintf(fout," %f\t %g\t %g\t %g\t %g\t\n",CellArr[id].r,CellArr[id].vbulk[0],CellArr[id].nH,x_nHII_shield, x_nHII_thin);
#endif

			x_nHII = (IncShield ==1 )? x_nHII_shield : x_nHII_thin;

			if (isnan(x_nHII) == 1)
			{
				printf("ERROR: NaN found after applying UV background ionization\n");
				dprinti(id);
				dprintd(x_nHII);
				dprinti(IncShield);
				dprintd(x_nHII_shield);
				dprintd(x_nHII_thin);
				dprintd(CellArr[id].nH);
				printf("Exiting...\n");
				exit(0);
			}

			CellArr[id].nH *= (1 - x_nHII);

	
//		printf("Gamma_shield %g Gamma_thin %g thin/shield %g\n",Gamma_z_shield,Gamma_z_thin,Gamma_z_thin/Gamma_z_shield);
//		exit(0);
						
		}
	}
#ifdef WRITEPROFILE

	fclose(fout);
	exit(0);

#endif

/*		for (i = 0;i<200;i++)
		{
			dprinti(i);
			dprintd(CellArr[i].nH);
		}
*/

}
