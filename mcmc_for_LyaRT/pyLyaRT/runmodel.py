import os 
import sys



#import matplotlib
# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')

#import matplotlib.pyplot as P
from subprocess import call
from struct import unpack
import numpy as np
from numpy.random import randn as randn





#### Define Geometry ####

#Geom='ThinShell'
#Geom='Wind'
Geom='Biconical_Wind'

#### Define Geometry ####


nspec=2

#### Deninig input and output paths, MUST be redifined by each user ####

source='harold'                             # currently choose between  'harold and 'mike_LAE' 
in_spectrum_dir="./useful_spectra/"         #input spectra are are inside this folder and their names have the structure: source + '_'  + nspec + '.data'
user_path="/Users/jemejia/"             #user main path
LyaRT_path=user_path+"LyaRT/src/"  #LyaRT source code

#### Deninig input and output paths, MUST be redifined by each user ####



####### pick which parameters will be minimized by mcmc, be aware of  uncommenting  one of the posibilities  below #######

parameters=['NH','Vmax','Z','theta']
parameters=['NH','Vmax','theta']
#parameters=['NH','Vmax']               # By default
#parameters=['NH','Vmax','Z']

####### pick which parameters will be minimized by mcmc, be aware of  uncommenting  one of the posibilities  above #######


#### Define ndim accordingly ####

if parameters==['NH','Vmax']:
	ndim=2
elif parameters==['NH','Vmax','Z']:
	ndim=3
elif parameters==['NH','Vmax','Z','theta']:
	ndim=4
	Geom='Biconical_Wind'
elif parameters==['NH','Vmax','theta']:
	ndim=3
	Geom='Biconical_Wind'
else:
	ndim=2

#### Define ndim accordingly ####




### Choose number of walkers and iterations to be performed by emceee ####

nwalkers,niter = 6, 200

### Choose number of walkers and iterations to be performed by emceee ####



#### CRUCIAL, DEFINE SEEDS FOR INPUT PARAMETERS!!!!!!! #####

logNH=18.5 #log(cm-2)
DlogNH=0.5 #Delta ini, for seeds

Vmax=600.0 #km/s
DVmax=300.0 #Delta ini, for seeds

logZ=-2.0 # log(Z/Zsun)
DlogZ=0.3

theta=0
Dtheta=np.pi/2
#np.rand

gauss_width=90.0 #km/s

scaling=np.arange(0.7,1.31,0.05) # after re-calibrating to line peak flux. 
shift=np.arange(-7.0,7.0,0.2) # \AA



NH=logNH+DlogNH*randn(nwalkers)
NH=np.abs(NH)
V=Vmax+DVmax*randn(nwalkers)
V=np.abs(V)
Z=logZ+DlogZ*randn(nwalkers)
TH=theta+Dtheta*np.random.random(nwalkers)

#### CRUCIAL, DEFINE SEEDS FOR INPUT PARAMETERS!!!!!!! #####

#GW=gauss_width+20.0*randn(nwalkers)
#GW=np.abs(GW)
#F=scaling + 0.15*randn(nwalkers)
#F=np.abs(F)
#S=shift + 1.0*randn(nwalkers)
#S=np.abs(S)





def gaussian(x,lambda_0=1215.57,sigma_v=90.0,lambda_up=1220.0, lambda_low=1210.0):
    

    """ definition of the guassian function"""
    minbins=20
    sigma_l= ( 1.0*sigma_v/(1.0*3e5) ) * lambda_0
    y = 1 / ( np.sqrt(2.0 * np.pi) * sigma_l )  * np.exp( -((x - lambda_0 + 0.0) / (2.0*sigma_l) )**2 )
    return y

 
def read_wave(filename):
    """
    Reads the output file of LyaRT and transforms it from bynary to arrays.
    This needs major update because the output has changed since 2015.
    
    PARAMS:
    filename: Name of the file to be read
    RETURNS:
    wl, wl0,interact,nscat, x,y,z,theta,phi
    """
    
    C=299792458
    KB=1.3806503e-23
    MP=1.67262158e-27
    T=10000.0
    Vth=np.sqrt(2*KB*T/MP);
    Kth=Vth/C;
    lambda_0=1215.668    
    d_lambda=0.136

    f=open(filename, "rb")
    f.read(8) #int Nphotons
    Np=unpack('l', f.read(8))[0]
    # print Np
    interact=np.empty(Np)
    for i in range(Np):  interact[i]= unpack('i', f.read(4))[0] 
    # reading interact array. (int array). 4 is the size of an int. and Np is 
    #If interact[i] is equal to: 
    #1 Last interaction with Hydrogen
    #2 Last interaction with Dust
    #3 Absorbed by dust at last interaction
    #4 No interaction 
    #f.read(4*Np) can be used in  the case that you dont need this array
    # and you want to save a bit of time when reading. It skips the next 
    #4*Np bytes which is the length of the interact array.  

    

    
    
    wl=np.empty(Np)
    for i in range(Np):  wl[i]= unpack('f', f.read(4))[0] 
    
    #f.read(4*Np) #nscat size (int array)
   
    nscat=np.empty(Np)
    for i in range(Np):  
        nscat[i]= unpack('f', f.read(4))[0] 
    


    #f.read(20*Np) #posarr (3Np) and angarr (2Np) size (int array )
    

    #pos array
    x=np.empty(Np)
    y=np.empty(Np)
    z=np.empty(Np)
    for i in range(Np):  
        x[i]= unpack('f', f.read(4))[0] 
        y[i]= unpack('f', f.read(4))[0] 
        z[i]= unpack('f', f.read(4))[0] 
        #print i
        #print Np
    
    #ang array
    phi=np.empty(Np)
    theta=np.empty(Np)
    for i in range(Np):  
        theta[i]= unpack('f', f.read(4))[0]
        phi[i]= unpack('f', f.read(4))[0] 




    
    wl0=np.empty(Np)
    for i in range(Np):  wl0[i]= unpack('f', f.read(4))[0] 
    f.close()    
    
    wl=lambda_0/(1+ wl*Kth)
    wl0=lambda_0/(1+ wl0*Kth)
    
    return wl, wl0,interact,nscat, x,y,z,theta,phi


def gen_par_file_LyaRT(logNH,Vmax,logZ,Geom,NPhotons=20000,xcrit=4.0,Temp = 0.0,dlambda=11.0,user="/Users/jemejia/",in_dir='LyaRT/data/Params/grid/',out_dir='mcmcrun/'):
    """
    Returns the parameter file to run LyaRT given the properties logNH, Vmax,logZ and Geom.
    
    PARAMS:
    logNH: log10 of the column density of Neutral hidrogen (expressed in units of cm**(-2))
    Vmax: Expansion velocity of the Thin Shell and maximal expansion velocity for the Wind. (expressed in km/s) 
    logZ: Metallicity in log10 scale expressed in Solar units: log10(Z/Zsun)
    Geom: Geometry to use: 'THIN' for the Thin Shell, 'WIND' for the simetric WIND and 'BICONE' of the Biconical Wind.
    NPhotons=20000
    xcrit=4.0
    Temp = 0.0
    user="/home/Users/jemejia/" --> Is the user path and the root for the directories below
    in_dir="LyaRT/data/Params/grid/" --> Is path were the input parameter files will be saved
    out_dir="mcmcrun/" --> In th path were the model files will be eventually saved after running the LyaRT code. 
    """


    
    #----------------Physical constants--------------
    Zsun=0.02
    C= 299792458
    KB= 1.3806503e-23 #Boltzman constant in MKS
    MP= 1.67262158e-27 #Hidrogen mass in Kg
    mH= 1.67262158e-24 #Hidrogen mass in g
    MSun=1.9891e33 # Sun mass in g
    yr_in_s=60*60*24*365.25 #year in seconds
    kpc=3.08e21
    km_s_in_cm_s=1e5
    T=10**(4+Temp) #Temp is in log units log10( T/(10**4K) )
    lambda_0=1215.67e-10
    dlambda=dlambda*1e-10 #in cm (1cm =1e-8A)
    Vth=np.sqrt(2*KB*T/MP)
    Kth=Vth/C
    #---------------------------------------------

    #----------Transforming into Lyart units-----------
    X_min = -0.5*dlambda/((lambda_0+0.5*dlambda)*Kth)
    X_max=0.5*dlambda/((lambda_0 - 0.5*dlambda)*Kth)
    dX=X_max-X_min
    #------------------------------------------------

    param_dir=user+in_dir
    out_dir=user+out_dir
    
    if not os.path.exists(param_dir):
        os.mkdir(param_dir)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
            
    #----Opening list file----------------
    file=open(param_dir+"file_list","w")
            
    #-----Converting Units------------------
    vmax_1=Vmax*km_s_in_cm_s #vmax in cm/s
    NH=np.power(10,logNH) #NH from log to normal
    Z=Zsun*np.power(10,1.0*logZ) #fractional abundace of metals
    """
    if logNH < 15:
	    xcrit = 2.0
    if logNH < 17 and logNH >= 15:
	    xcrit = 3.0
    if logNH < 19 and logNH >= 17:
	    xcrit = 5.0
    if logNH > 19:
	    xcrit = 6.0 + (logNH-19.0)
    """
		

    # ---------------------------Defining Radius and Ndot----------------------------------------
    if (Geom=="ThinShell"):
        RSphere=kpc
        R_inner=0.9*RSphere
        Ndot=8*3.141592*NH*RSphere*vmax_1*1e-40 # (particles/sec in units of 1e-40s)
        mean_nH=10.0*NH/(RSphere)
    if (Geom=="Wind"):
        R_inner=0.1*kpc
        RSphere=500.0*R_inner
        Ndot=4*3.141592*NH*R_inner*vmax_1*1e-40 # (particles/sec in units of 1e-40s)
        mean_nH=NH/R_inner
    if (Geom=="Biconical_Wind"):
        R_inner=0.1*kpc
        RSphere=500.0*R_inner
        Ndot=2.0*3.141592*NH*R_inner*vmax_1*1e-40 # (particles/sec in units of 1e-40s)
        mean_nH=NH/R_inner
                            
    model_name=''
    #---------------Naming and opening the parameter files---------------------------------------
 
    model_name="%s_log_nH%.4f_vmax%.2f_log_Z%.2f"  % ( Geom ,  logNH , Vmax ,  logZ )
 
    filename=param_dir+model_name+".data"
                        
    outfile=model_name+".data"
    f=open(filename,"w")
    #-------------------------------------------------------------------------------------    
    file.write(outfile+"\n") # Writing into a list the names of the parameter files
    #--------------------Very Important for running the models---------------------------#
    
    #------------------Writing the parameter file-----------------
    f.write("GeomName = "+Geom+"\n")
    f.write("Set_Tolerance = no\n")
    f.write("nout_max        ="+ str(int(NPhotons))+"  \n")
    f.write("np_max  =        "+ str(NPhotons) + "\n")
    f.write("np_min  =         " +str(int(1.0*NPhotons/100.0)) +" \n")
    f.write("NPhotons = "+str(int(NPhotons))+"\n")
    f.write("NCellX = 1000\n")
    f.write("NCellY = 1000\n")
    f.write("NCellZ = 1000\n")
    f.write("NCells =  1000\n")
    f.write("XPeriodic = 1\n")
    f.write("YPeriodic = 1\n")
    f.write("ZPeriodic = 1\n")
    f.write("IncUV = 0\n")
    f.write("xSize          =         1077639037059072.\n")
    f.write("ySize          =         1077639037059072.\n")
    f.write("zSize          =         1077639037059072.\n")
    f.write("Z = " + str(Z) + "\n")
    f.write("R_inner = "+str(R_inner)+"\n")
    f.write("RSphere = "+str(RSphere)+"\n")
    f.write("mean_nH = "+str(mean_nH)+"\n")
    f.write("ColDens = "+str(NH)+"\n")
    f.write("VarXCrit = 0\n")
    f.write("Temp = "+str(Temp)+"\n")
    f.write("vmax = "+str(Vmax)+"\n")
    f.write("Ndot = "+str(Ndot)+"\n")
    f.write("xcritval  = "+ str(xcrit)+ "\n")
    f.write("OutShort = "+out_dir+outfile+"\n")
    f.write("IncDust = Yes \n")
    f.write("OutMode = Short \n")
    f.write("dX = "+str(dX)+"\n")
    f.close()
#    print 'saved parfile   ',filename
    #----------------------------------------------------------
    return 0
    

def runmodel(logNH,Vmax,logZ,Geom,NPhotons=20000,xcrit=3.0,Temp = 0.0,user=user_path,in_dir='LyaRT/data/Params/grid/',out_dir='mcmcrun/',LyaRT_path="LyaRT/src/"):
    """
    Run an LyaRT model to be saved in "in out_dir". LyaRT.out code should be placed inside "LyaRT_path". "in_dir" is the directory where the input parameter
    files are saved. 
    
    PARAMS:

    logNH: log10 of the column density of Neutral hidrogen (expressed in units of cm**(-2))
    Vmax: Expansion velocity of the Thin Shell and maximal expansion velocity for the Wind. (expressed in km/s) 
    logZ: Metallicity in log10 scale expressed in Solar units: log10(Z/Zsun)
    Geom: Geometry to use: 'THIN' for the Thin Shell, 'WIND' for the simetric WIND and 'BICONE' of the Biconical Wind.
    NPhotons=20000 by default. Photons to be used by the LyaRT code for the MCMC calculations.
    xcrit=3.0 by default. Critical frequency to speed up the code. This moves resonant photons to the wings. See appendix A from Orsi et al 2012
    Temp = 0.0 by default. In log10 units where 10000K is the reference temparature.
    user="/home/Users/jemejia/" --> Is the user path and the root for the directories below
    in_dir="LyaRT/data/Params/grid/" --> Is path were the input parameter files will be saved
    out_dir="mcmcrun/" --> In th path were the model files will be eventually saved after running the LyaRT code. 
    LyaRT_path="LyaRT/src/"

    RETURN:

    Output file in bynary with this structure:
    
    "%s_log_nH%.4f_vmax%.2f_log_Z%.2f.data"  % ( Geom ,  logNH , Vmax ,  logZ )

    """
    LyaRT_path=user_path+LyaRT_path
    
    model_name="%s_log_nH%.4f_vmax%.2f_log_Z%.2f"  % ( Geom ,  logNH , Vmax ,  logZ )

    
    outname=user+out_dir+model_name+".data"
    filename=user+in_dir+model_name+".data"
    
    
    gen_par_file_LyaRT(logNH,Vmax,logZ,Geom,model_name=model_name)
    os.chdir(LyaRT_path)
    print(filename)
    call(['./LyaRT',filename,str(np.random.randint(100000))])
    

	    

        


