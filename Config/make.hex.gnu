INC_NETCDF = -I/home/nersc/pavelsa/local/include
LIB_NETCDF = /home/nersc/pavelsa/local/lib/libnetcdf.a
LIB_LAPACK = /home/nersc/pavelsa/local/lib/lapack.a /home/nersc/pavelsa/local/lib/tmglib.a /home/nersc/pavelsa/local/lib/blas.a 
LIB_FFT = /opt/fftw/3.2.1/lib/libfftw3.a

INCS = $(INC_NETCDF) -I/opt/fftw/3.2.2/include
LIBS = $(LIB_LAPACK) $(LIB_NETCDF) $(LIB_FFT)

ifeq ($(MPI),YES)
	CF90 = /home/nersc/pavelsa/local/bin/mpif90
	PARO =
	CPPFLAGS = -DQMPI
else
	CF90 = /home/nersc/pavelsa/local/bin/g95
	PAR0 =
endif
CF77 = $(CF90)
LD = $(CF90)
CPP = /usr/bin/cpp -traditional-cpp
CC = gcc

CPPARCH = -DIA32 -DFFTW -DNOMPIR8
CPPFLAGS += -P $(CPPARCH) -DF90_NOFLUSH -D_G95_

SIZEO = -r8
#OPTO = -O2 -Wall
OPTO = -Wall
#ARCHO = -fno-second-underscore
INLO =
DIVO =
DEBUG_FLAGS = -g

FFLAGS = $(SIZEO) $(OPTO) $(ARCHO) $(PARO) $(INLO) $(DIVO) $(DEBUG_FLAGS) $(INCS)
CFLAGS = $(FFLAGS) -Df2cFortran
LINKFLAGS = $(SIZEO) $(OPTO) $(ARCHO) $(PARO) $(INLO) $(DIVO) $(DEBUG_FLAGS)

# language-specific flags
#
F77FLG =
F90FLG =
