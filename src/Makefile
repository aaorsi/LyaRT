
# --Choose your system

#SYSTYPE = "Quintor"
#SYSTYPE = "Titania"
#SYSTYPE = "Quintor"
#SYSTYPE = "Miranda"

#SYSTYPE	= 'PUC-geryon'
SYSTYPE	= 'local'

# The following are specific for the code:
DEFFLAGS	= -DUSEREJECTION  -DTCONST -DHAPPROX -DTIMELIMIT -DGETPOSDIR  -DWRITEALL # -DXCONTINUUM # -DXINDIST  # -DTEST_RII  # -DWRITEPROFILE # -DTAUGUIDERDONI	# -DUSEREJECTION # -DWRITEPROFILE # -DWRITEALL # -DTEST_RII  

# USERREJECTION: Use rejection method to compute vp
# TCONST	   : Assume constant temperature --it simplifies the calculations
# HAPPROX	   : Use approximations for the function H(x,a) when neccessary
# TIMELIMIT	   : Defines a Time Limit at which the program stops 
# XINDIST	   : Initial x is drawn from a given distribution
# XCONTINUUM       :Initial flat distrubtion of photons in the range  [x-5,x+5]
# GETPOSDIR	   : In the output file, get additionally the position and direction of the escaping photon

ifeq ($(SYSTYPE),'PUC-geryon')
	CC        = gcc
#	OPTIMIZE = -O3 -m64 
	OPTIMIZE = -g
	CFLAGS  = $(OPTIMIZE)  $(GSL64INCL) $(GSL64LIB) $(DEFFLAGS) # -I/usr/local/include
	LIBS = -lm -lgsl -lgslcblas
endif

ifeq ($(SYSTYPE),'local')
	CC        = cc
#	OPTIMIZE = -O3 -m64 
	OPTIMIZE = -g -m64 
	CFLAGS  = $(OPTIMIZE)  $(GSL64INCL) $(GSL64LIB) $(DEFFLAGS) # -I/usr/local/include
	LIBS = -lm -lgsl -lgslcblas
endif


ifeq ($(SYSTYPE), "Quintor_DEBUG")
	CC        = cc -g -C
	OPTIMIZE = \
		-m64 -xarch=sparcvis2 
	CFLAGS  = $(OPTIMIZE) $(GSL64INCL) $(GSL64LIB) $(DEFFLAGS)
	LIBS = -lm -lgsl -lgslcblas
endif

ifeq ($(SYSTYPE), "Quintor")
	CC        = cc -C
	OPTIMIZE = -fast -xO3 \
		-dalign \
		-xbuiltin=%all \
		-m64 -xarch=sparcvis2 \
		-xtarget=ultra3i -xcache=64/32/4:1024/64/4 
	CFLAGS  = $(OPTIMIZE) $(GSL64INCL) $(GSL64LIB) $(DEFFLAGS)
	LIBS = -lm -lgsl -lgslcblas
endif

ifeq ($(SYSTYPE), "Miranda")
	CC        = cc -C
	OPTIMIZE = -fast -xO3 \
		-dalign \
		-xbuiltin=%all \
		-m64 -xarch=sparcvis2 \
		-xtarget=ultra3i -xcache=64/32/4:1024/64/4 
	CFLAGS  = $(OPTIMIZE) $(GSL64INCL) $(GSL64LIB) $(DEFFLAGS)
	LIBS = -lm -lgsl -lgslcblas
endif

ifeq ($(SYSTYPE), "Titania")
	CC        = cc -g -C
	OPTIMIZE = 	-fast \
		-xO5 \
		-xlibmil \
		-fsimple=2 \
		-dalign \
		-fns \
		-fsingle \
		-ftrap=%none \
		-xbuiltin=%all \
		-xtarget=ultra3 -xcache=64/32/4:8192/512/1 \
		-xopenmp=parallel
	CFLAGS  = $(OPTIMIZE)  
	LIBS    = -lm 
endif

ifeq ($(SYSTYPE), "CDM")
	CC        = gcc -g 
	OPTIMIZE = -03 -xarch=native
	CFLAGS   = -I/opt/local/gsl/32/current/include -L/opt/local/gsl/32/current/lib \
	-R/opt/local/gsl/32/current/lib \
	-DUSEREJECTION -DUVBACKGROUND # -DWRITEPROFILE # -DWRITEALL # -DTEST_RII  
	LIBS = -lm -lgsl -lgslcblas
endif


#---------------------------------------------------------

EXEC   = LyaRT

OBJS   = main.o  \
         	default_parameters.o read_parameters.o define_geometry.o \
			record_data_long.o record_data_short.o record_data_pos.o \
			get_H.o get_HG.o get_vp.o get_dipolar.o vp_rejection.o voigt.o lorentz.o \
			ionization_eqns.o numerical.o find_edgecone.o
			

INCL   =  proto.h allvars.h Makefile

#---------------------------------------------------------

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o  $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC) core

#---------------------------------------------------------


