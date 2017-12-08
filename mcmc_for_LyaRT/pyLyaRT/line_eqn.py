import sys
import os
from pylab import *
import scipy 
from scipy.stats import *
import scipy.stats
import numpy as np
#==========================================================================================================#
#                                            FUNCIONES UTILES                                              #
#==========================================================================================================#
def func_laplace( x , e , w, A ):
    return scipy.stats.laplace.pdf( x , e , w ) * A

def skew_func( x , e , w , a ):
    return skew( x , e , w , a)

def powerlaw( x , c , k ):
    return c * (x**k)

def lineal( x , c , k ):
    return c + x * k 

def func_fermi( x , a , b , c ):
    return( 1.5 / (np.e ** ( (a - x * 1. / c ) * 1. / b )   + 1) )

def pdf(x):                           
    return 1/np.sqrt(2*np.pi) * np.exp(-x**2/2) 
                                       
def cdf(x):                           
    return (1 + scipy.special.erf(x/np.sqrt(2))) / 2    
                                       
def skew(x,e=0,w=1,a=0):               
    # e = location                     
    # w = scale                       
    # a = shape                        
    t = (x-e) / w                      
    return 2. / w * pdf(t) * cdf(a*t)
#==========================================================================================================#
#                                FUNCIONES QUE DEVUELVEN LA FORMA DE LAS CURVAS                            #
#==========================================================================================================#
def line_MODELO_THIN_SHELL( x , V , logNH  ):  

        cAe_1 =  -3757.61118257
        cAe_2 =  222.31619595
        cAw_1 =  3352.46655007
        cAw_2 =  -147.488962509
        cAA_1 =  16228.7258176
        cAA_2 =  -730.922923683

        cBe_1  =  10.7499981248
        cBe_2  =  3.72021278078
        cBw_1  =  21.3684490763
        cBw_2  =  2.34548010086
        cBA_1  =  15.6408108851
        cBA_2  =  3.0346431816

        cCe_1  =  1.91744811355
        cCe_2  =  4.1947423119
        cCw_1  =  12.5179329655
        cCw_2  =  2.55604844932
        cCA_1  =  25.8674951863
        cCA_2  =  2.62103597907

        auxA1 = lineal( logNH , cAe_1 , cAe_2) 
        auxA2 = lineal( logNH , cAw_1 , cAw_2) 
        auxA3 = lineal( logNH , cAA_1 , cAA_2) 

        A = func_laplace( V , auxA1 , auxA2 , auxA3  )    

        auxB1 = powerlaw( logNH - 18 , cBe_1 , cBe_2) 
        auxB2 = powerlaw( logNH - 18 , cBw_1 , cBw_2) 
        auxB3 = powerlaw( logNH - 18 , cBA_1 , cBA_2)
 
        B = func_laplace( V , auxB1 , auxB2 , auxB3  )    

        auxC1 = powerlaw( logNH - 18 , cCe_1 , cCe_2) 
        auxC2 = powerlaw( logNH - 18 , cCw_1 , cCw_2) 
        auxC3 = powerlaw( logNH - 18 , cCA_1 , cCA_2)
 
        C = func_laplace( V , auxC1 , auxC2 , auxC3  )    


        e1 = -0.15       * V 
        e2 = -0.26521279 * V 
    
        w0 = (10**-29.71395395) * logNH ** 23.57131848  
        w1 =  1.3273309158902669e-20 * logNH ** 15.72863486  
        w2 = (10**-26.6994039) * logNH ** 21.20551954
    
        dist =  scipy.stats.cauchy( 0. , w0 )
    
        line = A * dist.pdf( x ) + B * scipy.stats.norm.pdf(x, e1, w1) + C * scipy.stats.norm.pdf(x, e2, w2) 
    
        line[ x>0 ] = 0 

        return line



def line_MODELO_BICONE( x , V , logNH  ):  

    e1 = -0.15234952462716664 * V 

    w0 = 2.16209054425e-32 * logNH  ** 25.5391295025

    w1_m =  0.00295447 * logNH **2 + (-0.12665552) * logNH + 1.38446258
    w1   = w1_m * V 

    AuxA1 =  10**( -37.6186895847 * np.log10( logNH ) + np.log10( 1.02587491746e+47 )  )
    AuxA2 = -2.17971807e-01 * logNH ** 3 +  1.31380212e+01 * logNH ** 2 +  (-2.63538905e+02) * logNH + 1.75998809e+03
    A = powerlaw( V , AuxA1 , AuxA2 )

    AuxB1 =  0.40084952 * logNH ** 2 + (-16.64101343) * logNH + 171.83117539
    AuxB2 = -1.95586868  * logNH ** 2 +   81.59370138  * logNH + (-846.67235739)
    AuxB3 =  2.29002188   * logNH ** 2 + (-95.99263695 ) * logNH +  1001.49132297
    B  = AuxB1 * np.log10(V) ** 2 + AuxB2 * np.log10(V) + AuxB3

    ca_s1 =   0.14260217 * logNH ** 2 + (-6.01342823) * logNH +   63.48248556
    ca_s2 =  -0.06754816 * logNH ** 2 +   2.79164084  * logNH + (-28.27873284)
    a = powerlaw( V , ca_s1 , ca_s2)

    line = A * scipy.stats.laplace.pdf( x , 0 , w0 ) + B * skew_func( x , e1 , w1 , a ) 

    line[ x>0 ] = 0 

    return line

def line_MODELO_WIND( x , V , logNH  ):

    e1 = -0.15234952462716664 * V
    e2 = -0.06302857 * V
    e3 = -0.26521279 * V

    w2 = logNH ** 24.0531610395 * 3.71483556417e-31

    w1 = logNH ** 20.19 * 4.1e-26

    w0 = logNH ** 24.9951464064 * 7.42889429055e-32

    w3 = logNH ** 24.0531610395 * 3.71483556417e-31

    cDe_1 =  3.81535580405e-28
    cDe_2 =  22.3870363827
    cDw_1 =  1.68910345201e-26
    cDw_2 =  21.1688985271
    cDa_1 =  1.7784258
    cDA_1 =  7.71573781884e-26
    cDA_2 =  20.8655698539

    CDe = powerlaw( logNH , cDe_1 , cDe_2 )
    CDw = powerlaw( logNH , cDw_1 , cDw_2 )
    CDa = cDa_1
    CDA = powerlaw( logNH , cDA_1 , cDA_2 )

    D = skew( V , CDe , CDw , CDa   ) * CDA

    cCe_1 =  2.66559427835e-29
    cCe_2 =  23.6807531568
    cCw_1 =  7.0084912688e-22
    cCw_2 =  17.9108165532
    cCa_1 =  6.0
    cCA_1 =  4.11593936549e-26
    cCA_2 =  20.8422256546

    CCe = powerlaw( logNH , cCe_1 , cCe_2 )
    CCw = powerlaw( logNH , cCw_1 , cCw_2 )
    CCa = cCa_1
    CCA = powerlaw( logNH , cCA_1 , cCA_2 )

    C = skew( V , CCe , CCw , CCa   ) * CCA

    cBe_1 =  1.8822597655e-30
    cBe_2 =  24.3724521539
    cBw_1 =  2.0892604197e-22
    cBw_2 =  18.4016637192
    cBa_1 =  8.85294425002
    cBA_1 =  1.29151062692e-23
    cBA_2 =  19.3676286499

    CBe = powerlaw( logNH , cBe_1 , cBe_2 )
    CBw = powerlaw( logNH , cBw_1 , cBw_2 )
    CBa = cBa_1
    CBA = powerlaw( logNH , cBA_1 , cBA_2 )

    B = skew( V , CBe , CBw , CBa   ) * CBA

    cAm_1 =  2.51987345608e-36
    cAm_2 =  28.5338257734

    CAm = powerlaw( logNH , cAm_1 , cAm_2 )

    A = func_fermi( V , 8 , 1 , CAm  )

    line = A * scipy.stats.laplace.pdf( x , 0 , w0 ) + B * scipy.stats.cauchy.pdf( x , e1 , w1 ) + C * scipy.stats.norm.pdf( x , e2 , w2 ) + D * scipy.stats.norm.pdf( x , e3 , w3 )

    line[ x>0 ] = 0


    return line

#==========================================================================================================#
#                             Inicializa algunos parametros importantes                                    #
#==========================================================================================================#
def RT_parameters():
    T4 = 1. # = 10000. / 1e4
    nu0 = 2.46777 * 1.e15 #3. * 10.**8 / (1215.67 * (10**(-10)))
    Vth = 12.85 * np.sqrt(T4) 
    Dv = Vth * nu0 *1. / ( 3 * (10**5))

    return T4 , nu0 , Vth , Dv 
#==========================================================================================================#
#                                                LINE PROFILES                                             #
#==========================================================================================================#
def line_profile( X , V_exp , logNH , GEONAME ):
    '''
    Returns an analitycal prediction of the Ly-alpha emission profile. This is just an approximation, for
    accuracy you need to run the LyaRT code. See Orsi et al 2012 and Gurung et al. 2018 (In prep).

    PARAMS:
    X: Shift in frequency around nu_Lya in units of the termal width: x (nu - nu_Lya)/dnu_th. Where dnu_th=(v_th/c)*nu_Lya and v_th=(2*k_B*T/mp)***(1/2)
    V_exp: Expansion velocity of the Thin Shell and maximal expansion velocity for the Wind. (expressed in km/s) 
    LogNH: log10 of the column density of Neutral hidrogen (expressed in units of cm**(-2))
    GEOMNAME: Geometry to be used,  'THIN' for the Thin Shell, 'WIND' for the simetric outflowing Wind and 'BICONE' of the Biconical outflowing Wind.
    '''

    T4 , nu0 , Vth , Dv = RT_parameters()

    #X = np.linspace(-500,50,Nx) 


    if GEONAME == 'THIN'   :
        distri = line_MODELO_THIN_SHELL( X , V_exp , logNH )

    if GEONAME == 'BICONE'   :
        distri = line_MODELO_BICONE( X , V_exp , logNH )

    if GEONAME == 'WIND'   :
        distri = line_MODELO_WIND( X , V_exp , logNH )

    Area = np.trapz( np.absolute(distri) , X )
    distri = np.absolute( distri ) *1. / Area

    return distri 

