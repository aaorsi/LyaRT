
// Debugging
#define dprinti(expr) printf(#expr " = %d\n", expr)
#define dprintl(expr) printf(#expr " = %ld\n", expr)
#define dprintf(expr) printf(#expr " = %f\n", expr)
#define dprintd(expr) printf(#expr " = %g\n", expr)
#define dprints(expr) printf(#expr " = %s\n", expr)


// Functions
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MIN3(a,b,c) ((MIN(a,b))<(c)?(MIN(a,b)):(c))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define SQR(a) ((a)*(a))

// Numerical Constants
#define Pi 3.141592654
#define sqrtPi 1.7724539041519165
#define NCELLMAX 1000000		
#define NCHARMAX 500			
#define NTABARR  5000						
#define NTAB1	 5000
#define NTAB2	 50000
#define KnH		 1e35
#define MAXX	 20			// This is the maximum value of x in te lookup table to get H(x)
#define NLINES	 10		
#define MAXTIME  1440		// Minutes

// Physical constants
#define zstar 0.02				// Solar Metallicity
#define ext_zstar 8.87e-20		// Some ratio related to zstar
#define Albedo 0.39				// Albedo of dust grains at Ly-alpha centre
#define sx_const 5.8689e-14		// cross section at Ly-alpha freq.
#define alpha_A 4.18e-13		// case A recombination coefficient, T=10000K
#define Planck 6.6262e-27		// [erg s]
#define KB 1.3806503e-16
#define MP 1.67262158e-24
#define L_REST 1.215668e-5



#ifdef TAUGUIDERDONI

#define spar		1.35
#define Ext_Ratio	3.43
#define NHConst		2.1e21

#endif

// Variables
#ifdef HAPPROX

#define hc1  0.855
#define hc2  3.42
#define hc3  5.674
#define hc4  9.207
#define hc5  4.421
#define hc6  0.1117

#endif

typedef struct {
    long id;
	float p_lya;
	double T;
	double nH;
	float z;
	double vbulk[3];
	double r;
// neighbour: Each of the following neighbouring cells
// [-x,+x,-y,+y,-z,+z]
	long neig[6];
} cell;



#define NZZ	50
#define NEN 500

#define minlogryd 0
#define maxlogryd 100


typedef struct {
	double flux;
	double logryd;	//log energy [Ryd]
	double sigmaHI;
	double sigma_ratio;
} spec;

	long UV_Nz, UV_Ne,Nspec;
	int IncShield,IncUV;
	char GHIFile[NCHARMAX],UVzFile[NCHARMAX],UVeFile[NCHARMAX], \
		 UVspecFile[NCHARMAX];

	char OutDir[NCHARMAX],tabledir[NCHARMAX],sxfile[NCHARMAX],HFile[NCHARMAX];
        char IncDust[NCHARMAX];
        
	char GeomName[NCHARMAX],OutMode[NCHARMAX],OutLong[NCHARMAX],OutShort[NCHARMAX],Set_Tolerance[NCHARMAX];
    long NPhotons;
long NCells,NCellX, NCellY,NCellZ,NCells,nscatmax,nscat,NCellSphere,NCellsIn,NCellsOut;
	int XPeriodic,YPeriodic,ZPeriodic,Static;
double xSize,ySize,zSize,RSphere,H,vmax,R_Peak,dr,dr1,dr2;
double Temp,mean_nH,mean_z,xcritval,mean_nH_static,ColDens;
	double ni,nj,nk,vbulk_x,vbulk_y,vbulk_z,th0,ph0,radius,vbulk0_x,vbulk0_y,vbulk0_z;
	double c,x_test;
	long nHList, nxlist, nvplist, ndip, nHG,nout_max,np_max,np_min;		
	double R_inner,Ndot,b,R_Static; 
        float redshift, a_par;
        double dX;
        int StartAtCenter,GetVpRej,DefinePosAndDir;
        int VarXCrit;
	double Ix0,Iy0,Iz0,Ith0,Iph0;
	double Vrot, x_mean;
        float *XArr,*X0Arr,*PosArr,*AngleArr,*PosArrLong;
	int *InterArr, *NscatArr;
	double alpha_vprof, sigma_vprof,app_angle;
