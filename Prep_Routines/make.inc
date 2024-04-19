
ifeq ($(MPI),YES)
	CF90 = ifort
	PARO = -Mprof=dwarf -Minform=inform
	CPPFLAGS = -DQMPI 
else	
	CF90 = ifort
	PAR0 = 
endif
CF77 = $(CF90)
CC = icc
LD =  $(CF90)

# CPP, used internally by compilers - note that some compilers
# need special options for this, consult compiler manual if in trouble
CPP= /usr/bin/cpp -traditional-cpp -P 
CPPARCH = -DIA32 -DREA8 -DSERIA_IO -DTIMER -DIARGC -DFFTW -DIFPORT 
CPPFLAGS += $(CPPARCH) 
#CPPFLAGS += $(CPPARCH) -DF90_NOFLUSH 
#CPPFLAGS = -P $(CPPARCH) -std=gnu9x


# Flags for fixed/free format
#F90FLG= -free
#F77FLG = -nofree


FFLAGS = -g -convert big_endian -assume byterecl -cm -vec_reprot0 -w  -O2 -r8 -g -traceback 

CFLAGS= -O -shared-intel -mcmodel=medium
#CFLAGS= -xCORE-AVX2 -O2 -no-prec-div -qopt-prefetch -auto-p32 -no-ansi-alias -qopt-mem-layout-trans=2 


# Linker flags

LINKFLAGS= $(FFLAGS) 


#Libraries, Netcdf, FFTW, LAPACK

#LIBS = -lfftw3 -lgsl -lnetcdff -lnetcdf
LIBS = -lfftw3 -lnetcdff -lnetcdf 

