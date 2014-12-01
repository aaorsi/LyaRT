from struct import unpack
import numpy as np



def read_wave(filename):
    C=299792458
    KB=1.3806503e-23
    MP=1.67262158e-27
    T=10000.0
    Vth=np.sqrt(2*KB*T/MP);
    Kth=Vth/C;
    lambda_0=1215.668    
    
    
    
    f=open(filename, "rb")
    Nout=unpack('l', f.read(8) )[0]
   #int  N photons that scaped the galaxy. 8 bites is the size of a longint
    Np=unpack('l', f.read(8))[0] # Number of input photons
        

    
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

    
    freq=np.empty(Np)
    for i in range(Np):  freq[i]= unpack('f', f.read(4))[0] 
    #Output frequency of the photon (in x units). 4 is the lenght of a float 


    
    nscat=np.empty(Np)
    for i in range(Np):  nscat[i]= unpack('i', f.read(4))[0] 
    # nscat[i] is number of scatterings before last interaction    
    
    
    
    x=np.empty(Np)
    y=np.empty(Np)
    z=np.empty(Np)
    for i in range(Np):  
        x[i]= unpack('f', f.read(4))[0] 
        y[i]= unpack('f', f.read(4))[0] 
        z[i]= unpack('f', f.read(4))[0] 
    
    phi=np.empty(Np)
    theta=np.empty(Np)
    for i in range(Np):  
        theta[i]= unpack('f', f.read(4))[0]
        phi[i]= unpack('f', f.read(4))[0] 
    
    #Reading x,y,z, theta and phi of the photons at last scattering    
    #f.read(20*Np) can be used  if you dont need the position (3Np*4) and angles at last scattering (2Np*4) size (int array )
    freq0=np.empty(Np)
    for i in range(Np):  freq0[i]= unpack('f', f.read(4))[0] 
    f.close()    
    #initial frequency of the photons in x units. They are supposed to be generated at the center of the galaxy
    
    wl=lambda_0/(1+ freq*Kth)
    wl0=lambda_0/(1+ freq0*Kth) 
    #converting frequency in x units to lambda in angstrons
    #if you  need the frequency in the x units just return freq and freq0.
    
    return wl, wl0,x,y,z,theta,phi,Np,Nout,nscat, interact


#filename="/Users/jemejia/new_grid/Wind_log_nH22.0_vmax10.0_SFR10.0_log_Z-4.0.data"
#filename='/home/aaorsi/work/LyaRT/Output/Bicone_test/test_bicone.data'
#wl, wl0,x,y,z,theta,phi,Np,Nout,nscat, interact=read_wave(filename)
#print wl, wl0,x,y,z,np.max(theta),np.max(phi),Np,Nout,nscat, interact
