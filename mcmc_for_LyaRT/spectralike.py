import os 
import sys
from scipy.stats import chisquare
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integrate
from scipy.interpolate import UnivariateSpline

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import matplotlib.pyplot as P
from subprocess import call
from struct import unpack
import numpy as np
import emcee
from numpy.random import randn as randn
import triangle
from emcee.utils import MPIPool


use_ptsampler = False

if use_ptsampler:
	from emcee import PTSampler
	print 'running emcee PTSampler'


#### Define Geometry ####

#Geom='ThinShell'
#Geom='Wind'
Geom='Biconical_Wind'

#### Define Geometry ####


nspec=2

#### Deninig input and output paths, MUST be redifined by each user ####

source='harold'                             # currently choose between  'harold and 'mike_LAE' 
in_spectrum_dir="./useful_spectra/"         #input spectra are are inside this folder and their names have the structure: source + '_'  + nspec + '.data'
user_path="/home/CEFCA/aaorsi/"             #user main path
LyaRT_path="/home/CEFCA/aaorsi/LyaRT/src/"  #LyaRT source code

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
        minbins=20

        """ definition of the guassian function"""
	sigma_l= ( 1.0*sigma_v/(1.0*3e5) ) * lambda_0
	y = 1 / ( np.sqrt(2.0 * np.pi) * sigma_l )  * np.exp( -((x - lambda_0 + 0.0) / (2.0*sigma_l) )**2 )
	return y

 
def read_wave(filename):
    
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


def gen_par_file_LyaRT(logNH,Vmax,gauss_width,logZ,Geom,dlambda=11.0,user=user_path,in_dir='LyaRT/data/Params/grid/',out_dir='mcmcrun/',NPhotons=20000,xcrit=4.0,Temp = 0.0):

    
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
    Z=Zsun*np.power(10,1.0*logZ) #abundace of metals (%)
    
    if logNH < 15:
			xcrit = 2.0
    if logNH < 17 and logNH >= 15:
			xcrit = 3.0
    if logNH < 19 and logNH >= 17:
			xcrit = 5.0
    if logNH > 19:
			xcrit = 6.0 + (logNH-19.0)

		

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
                            
                            
    #---------------Naming and opening the parameter files---------------------------------------
    
    model_name="%s_log_nH%.4f_vmax%.2f_log_Z%.2f"  % ( Geom ,  logNH , Vmax ,  logZ )
    if (Geom=="Biconical_Wind"):
        model_name="%s_log_nH%.4f_vmax%.2f_log_Z%.2f"  % ( 'BWind' ,  logNH , Vmax ,  logZ )
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
    

def spectralike(logNH,Vmax,gauss_width,logZ,scaling,shift,theta,Geom,wl_obs,f_obs,f_err,specnumber,wl_min=1195.0,wl_max=1237.0,wl_0=1215.668,theta_bins=5,user="/home/CEFCA/aaorsi/",in_dir="LyaRT/data/Params/grid/",out_dir="mcmcrun/"):
    if Vmax>1500.0 or logNH>21.5 or Vmax<10.0 or logNH < 14 or logZ<-5.0 if theta<0 if theta>np.pi/2:
        chi2=np.inf
    else:
        out_dir=user+out_dir
        in_dir=user+in_dir
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
    
        d_wl=np.mean(wl_obs[1:]-wl_obs[:-1])
    
 
        arg_fobs_max=np.argmax(f_obs)
        fmax_obs=f_obs[arg_fobs_max]
        centroid_obs=np.sum( np.multiply(f_obs,wl_obs) )/np.sum(f_obs)
        
        model_name="%s_log_nH%.4f_vmax%.2f_log_Z%.2f"  % ( Geom ,  logNH , Vmax ,  logZ )
        if (Geom=="Biconical_Wind"):
            model_name="%s_log_nH%.4f_vmax%.2f_log_Z%.2f"  % ( 'BWind' ,  logNH , Vmax ,  logZ )
    
        outname=out_dir+model_name+".data"
        filename=in_dir+model_name+".data"
        no=1
        while no:
            if os.path.exists(outname) and os.path.isfile(outname):
                no=0
#                print 'saved codefile   ',outname
            else:
                gen_par_file_LyaRT(logNH,Vmax,gauss_width,logZ,Geom)
#                call(['cd '+LyaRT_path])
                os.chdir(LyaRT_path)
                call([LyaRT_path+'LyaRT',filename,str(np.random.randint(100000))])
#                print 'abrio'        
#       sys.exit(0)
                                    
        wavelt,wavel0t,interact,nscat, x,y,z,polar,azim=read_wave(outname)
	if Geom=='Biconical_Wind':
		delta_th=np.pi/(2.0*theta_bins)
		theta_values=np.arange(0,np.pi/2,delta_th)
		argth=np.argmin(np.abs(theta-theta_values))
		theta_out=theta_values[argth]
		if argth< len(theta_values)-1:
			thl=theta_values[argth]
			thm=theta_values[argth+1]
			
		else:
			thl=theta_values[argth-1]
			thm=theta_values[argth]
		wl=polar>thl
		wm=polar[wl]<thm
		wavelt1=wavelt[wl][wm]
		wavel0t1=wavel0t[wl][wm]
		interact1=interact[wl][wm]
		nscat1=nscat[wl][wm]
		x1=x[wl][wm]
		y1=y[wl][wm]
		z1=z[wl][wm]
		polar1=polar[wl][wm]
		azim1=azim[wl][wm]


		theta_values=np.arange(np.pi/2,np.pi,delta_th)
		argth=np.argmin(np.abs(np.pi-theta-theta_values))
		if argth> 0:
			thl=theta_values[argth-1]
			thm=theta_values[argth]
			
		else:
			thl=theta_values[argth]
			thm=theta_values[argth+1]
		
		thl=theta_values[argth-1]
		thm=theta_values[argth]
		wl=polar>thl
		wm=polar[wl]<thm
		wavelt2=wavelt[wl][wm]
		wavel0t2=wavel0t[wl][wm]
		interact2=interact[wl][wm]
		nscat2=nscat[wl][wm]
		x2=x[wl][wm]
		y2=y[wl][wm]
		z2=z[wl][wm]
		polar2=polar[wl][wm]
		azim2=azim[wl][wm]
	
	
		wavelt=np.append(wavelt1,wavelt2)
		wavel0t=np.append(wavel0t1,wavel0t2)
		interact=np.append(interact1,interact2)
		nscat=np.append(nscat1,nscat2)
		x=np.append(x1,x2)
		y=np.append(y1,y2)
		z=np.append(z1,z2)
		polar=np.append(polar1,polar2)
		azim=np.append(azim1,azim2)
	
		
		
	
	

        dust_abs=interact==3
        wavel=wavelt[~dust_abs]
        wavel0=wavel0t[~dust_abs]
        if len(wavel) == 0:
		print 'zero-sized output array found in '+outname
		sys.stdout.flush()
		chi2=1000000000
		return -1*chi2
	
        wl_max=np.amax(wavel)
        wl_min=np.amin(wavel)
        nbins=np.int( (wl_max - wl_min) / d_wl )
        wl_teo=np.array([wl_min + d_wl*x for x in range(nbins+1)])
        wl_hist= wl_teo  - 1.0*d_wl/2.0
        wl_hist=np.append(wl_hist, wl_hist[-1] + d_wl )
        out_weights=gaussian(wavel0,lambda_0=1215.668,sigma_v=gauss_width,lambda_up=np.amax(wavel0), lambda_low=np.amin(wavel0))
	out_weights_all=gaussian(wavelt,lambda_0=1215.668,sigma_v=gauss_width,lambda_up=np.amax(wavelt), lambda_low=np.amin(wavelt))
        flux_teo, wl=np.histogram(wavel,bins=wl_hist,normed=False,weights=out_weights)
        
        centroid_teo=np.sum(np.multiply(flux_teo,wl_teo))/np.sum(flux_teo)
        fmax_teo=np.amax(flux_teo)
        f_teo_err = np.sqrt( flux_teo)*fmax_obs/fmax_teo              
        f_teo=flux_teo*fmax_obs/fmax_teo              
    
        arg_max_teo=np.argmax(flux_teo)
        arg_cent_teo=np.argmin( np.abs(wl_teo-centroid_teo) )
        arg_cent_obs=np.argmin(np.abs(wl_obs-centroid_obs))
    
        #delta_low=np.amin((arg_cent_teo,arg_cent_obs))
        delta_low=np.amin((arg_max_teo,arg_fobs_max))
#	f_teo_func = interp1d(wl_teo, f_teo, kind='cubic',bounds_error=False,fill_value=0)                         
	if len(wl_teo) < 2:
		f_teo_func = interp1d([wl_min,wl_max],[0.0,0.0],bounds_error=False,fill_value=0)
	else:	
		f_teo_func = interp1d(wl_teo, f_teo, kind='linear',bounds_error=False,fill_value=0)                         

	chi2=1000000000
	for scale in scaling:
		for dx in shift:
			wl_obs_model=wl_obs + wl_teo[arg_max_teo]- wl_obs[arg_fobs_max] + dx  
			len_teo=len(wl_teo);len_obs=len(wl_obs_model)
			f_teo=f_teo_func(wl_obs_model)
			f_teo=scale*f_teo
			num_chi=f_teo - f_obs            
			den_chi=np.power(f_err,2)
			deg_free=np.float(len(f_obs))
			chi2a=1.0*np.sum(np.divide(np.power(num_chi,2),den_chi))/(deg_free-3)
			chi2=np.min([chi2,chi2a])
			if np.isnan(chi2): chi2=np.inf
	fn = "spec_output_"+str(specnumber)+".out"
	f = open(fn, "a")
	escape_fraction=len(wavel)/(1.0*len(wavelt))
	if Geom=='Biconical_Wind':
		f.write(" ".join([str(logNH),"\t"+str(Vmax),"\t"+str(logZ),"\t"+str(theta_out),"\t"+str(escape_fraction),"\t"+str(chi2)]))
	else:
		f.write(" ".join([str(logNH),"\t"+str(Vmax),"\t"+str(logZ),"\t"+str(escape_fraction),"\t"+str(chi2)]))
	f.write("\n")
	f.close()

    print 'chi2 = ',chi2
    


    sys.stdout.flush()
    return -1.0*chi2
                                    

def lnlike(params,Geom,wl_obs,f_obs,f_err,gauss_width,logZ,theta,scaling,shift,specnumber):
	logNH,Vmax=params
	return spectralike(logNH,Vmax,gauss_width,logZ,scaling,shift,theta,Geom,wl_obs,f_obs,f_err,specnumber)

def lnlikeZ(params,Geom,wl_obs,f_obs,f_err,gauss_width,theta,scaling,shift,specnumber):
	logNH,Vmax,logZ=params
	return spectralike(logNH,Vmax,gauss_width,logZ,scaling,shift,theta,Geom,wl_obs,f_obs,f_err,specnumber)



###### NEED TO BE REDEFINED TO ACCOUNT FOR THE ANGLE ######
def lnlikeZTH(params,Geom,wl_obs,f_obs,f_err,gauss_width,logZ,theta,scaling,shift,specnumber):
	logNH,Vmax,logZ,theta=params
	return spectralike(logNH,Vmax,gauss_width,logZ,scaling,shift,theta,Geom,wl_obs,f_obs,f_err,specnumber)

def lnlikeTH(params,Geom,wl_obs,f_obs,f_err,gauss_width,theta,logZ,scaling,shift,specnumber):
	logNH,Vmax,theta=params
	return spectralike(logNH,Vmax,gauss_width,logZ,scaling,shift,theta,Geom,wl_obs,f_obs,f_err,specnumber)
###### NEED TO BE REDEFINED TO ACCOUNT FOR THE ANGLE ######




wl_obs,f_obs,f_err=np.genfromtxt(in_spectrum_dir+ source + "_" + str(nspec) + ".data",unpack=1)

#----rescaling flux and error. Keep this order-----#
f_err=f_err/np.sum(f_obs)
f_obs=f_obs/np.sum(f_obs)
#----rescaling flux and error. Keep this order-----#
            
#chi2=spectralike(logNH,Vmax,gauss_width,logZ,scaling,shift,Geom,wl_obs,f_obs,f_err)
#print chi2






pool = MPIPool(loadbalance = True)

if not pool.is_master():
	pool.wait()
	sys.exit(0)



#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike, args=(Geom,wl_obs,f_obs,f_err,gauss_width,logZ,scaling,shift,nspec),threads=nwalkers)
if use_ptsampler:
	# use a flat prior
	def logp(x):
		return 0.0

	ntemps = 20
	
	#SAMPLER AND SEEDS ARE DEFINED ACCORDING TO THE MCMC PARAMETERS CHOSEN
	if parameters==['NH','Vmax']:
		pos = [[NH[i],V[i]] for i in range(nwalkers)]
		sampler = PTSampler(ntemps,nwalkers,ndim,lnlike,logp,loglargs=(Geom,wl_obs,f_obs,f_err,gauss_width,theta,logZ,scaling,shift,nspec),pool=pool)
	elif parameters==['NH','Vmax','Z']:
		pos = [[NH[i],V[i],Z[i]] for i in range(nwalkers)]
		sampler = PTSampler(ntemps,nwalkers,ndim,lnlikeZ,logp,loglargs=(Geom,wl_obs,f_obs,f_err,gauss_width,theta,scaling,shift,nspec),pool=pool)
	elif parameters==['NH','Vmax','Z','theta']:
		pos = [[NH[i],V[i],Z[i],TH[i]] for i in range(nwalkers)]
		sampler = PTSampler(ntemps,nwalkers,ndim,lnlikeZTH,logp,loglargs=(Geom,wl_obs,f_obs,f_err,gauss_width,scaling,shift,nspec),pool=pool)
	elif parameters==['NH','Vmax','theta']:
		pos = [[NH[i],V[i],TH[i]] for i in range(nwalkers)]
		sampler = PTSampler(ntemps,nwalkers,ndim,lnlikeTH,logp,loglargs=(Geom,wl_obs,f_obs,f_err,gauss_width,logZ,scaling,shift,nspec),pool=pool)
	else:
		pos = [[NH[i],V[i],Z[i]] for i in range(nwalkers)]
		sampler = PTSampler(ntemps,nwalkers,ndim,lnlikeZ,logp,loglargs=(Geom,wl_obs,f_obs,f_err,gauss_width,theta,logZ,scaling,shift,nspec),pool=pool)

	
else:
	#SAMPLER AND SEEDS ARE DEFINED ACCORDING TO THE MCMC PARAMETERS CHOSEN
	if parameters==['NH','Vmax']:
		pos = [[NH[i],V[i]] for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike, args=(Geom,wl_obs,f_obs,f_err,gauss_width,logZ,scaling,shift,nspec),pool=pool)
	elif parameters==['NH','Vmax','Z']:
		pos = [[NH[i],V[i],Z[i]] for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlikeZ, args=(Geom,wl_obs,f_obs,f_err,gauss_width,scaling,shift,nspec),pool=pool)
	elif parameters==['NH','Vmax','Z','theta']:
		pos = [[NH[i],V[i],Z[i],TH[i]] for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlikeZTH, args=(Geom,wl_obs,f_obs,f_err,gauss_width,scaling,shift,nspec),pool=pool)
	elif parameters==['NH','Vmax','theta']:
		pos = [[NH[i],V[i],TH[i]] for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlikeTH, args=(Geom,wl_obs,f_obs,f_err,gauss_width,scaling,shift,nspec),pool=pool)
	else:
		pos = [[NH[i],V[i]] for i in range(nwalkers)]
		sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike, args=(Geom,wl_obs,f_obs,f_err,gauss_width,scaling,shift,nspec),pool=pool)



print 'starting burn-in iterations'
sys.stdout.flush()
pos0, prob, state = sampler.run_mcmc(pos, 5)
print 'done'
sys.stdout.flush()

sampler.reset()


#f = open("chain.dat", "w")
#f.close()


#for result in sampler.sample(pos_end, iterations=10, storechain=False):
#    position = result[0]
#    f = open("chain.dat", "a")
#    for k in range(position.shape[0]):
#        f.write("{0:4d} {1:s}\n".format(k, " ".join(position[k])))
#    fc.lose()




fn = "spec_mcmc_new_"+str(nspec)+".out"
f = open(fn, "w")
f.close()
print 'now running sampler with '+np.str(niter * nwalkers)+' total steps'
sys.stdout.flush()
for pos, prob, rstate in sampler.sample(pos0, prob, state, iterations=niter):
    # Write the current position to a file, one line per walker                                                                                                                                                                                                                
    f = open(fn, "a")
    f.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
    #f.write("\n".join(["\t".join( pos.tolist()[i] +[prob.tolist()[i]] )  for i in range( len(prob) ) ] ) )
    f.write("\n")
    f.close()

pool.close()

print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))


#Plotting triangle plot that shows the model density  over each posible tuples of parameter projections
samples = sampler.chain.reshape((-1, ndim))
if parameters==['NH','Vmax']:
	fig = triangle.corner(samples, labels=[r"$\log(N_H)$", r"$V_{\rm exp}$"])
elif parameters==['NH','Vmax','Z']:
	fig = triangle.corner(samples, labels=[r"$\log(N_H)$", r"$V_{\rm exp}$",r"$\log(Z/Z_{\odot})$"])
elif parameters==['NH','Vmax','Z','theta']:
	fig = triangle.corner(samples, labels=[r"$\log(N_H)$", r"$V_{\rm exp}$",r"$\log(Z/Z_{\odot})$",r"$\theta$"])
elif parameters==['NH','Vmax','theta']:
	fig = triangle.corner(samples, labels=[r"$\log(N_H)$", r"$V_{\rm exp}$",r"$\theta$"])
else:
	fig = triangle.corner(samples, labels=[r"$\log(N_H)$", r"$V_{\rm exp}$"])
	
fig.savefig("triangle.png")





"""
try:
    import matplotlib.pyplot as pl
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    print("Try installing matplotlib to generate some sweet plots...")
else:
    pl.figure()
    for k in range(nwalkers):
        pl.plot(sampler.chain[k, :, 0])
    pl.xlabel("time")
    pl.savefig("eggbox_time.png")

    pl.figure(figsize=(8,8))
    x, y = sampler.flatchain[:,0], sampler.flatchain[:,1]
    pl.plot(x, y, "ok", ms=1, alpha=0.1)
    pl.savefig("eggbox_2d.png")

    fig = pl.figure()
    ax = fig.add_subplot(111, projection="3d")

    for k in range(nwalkers):
        x, y = sampler.chain[k,:,0], sampler.chain[k,:,1]
        z = sampler.lnprobability[k,:]
        ax.scatter(x, y, z, marker="o", c="k", alpha=0.5, s=10)
    pl.savefig("eggbox_3d.png")
"""



