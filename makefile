MPI = YES

include ./make.inc

SHELL = /bin/bash
VPATH = .:TMP

PROGS = EnKF

all: $(PROGS)

ENKF_SRC_F90 = \
qmpi.F90\
m_parameters.F90\
m_Generate_element_Si.F90\
mod_analysisfields.F90\
m_confmap.F90\
mod_measurement.F90\
m_oldtonew.F90\
m_random.F90\
m_spherdist.F90\
distribute.F90\
m_bilincoeff.F90\
m_get_mod_fld.F90\
m_get_mod_grid.F90\
m_get_mod_nrens.F90\
m_insitu.F90\
m_local_analysis.F90\
m_obs.F90\
m_parse_blkdat.F90\
m_pivotp.F90\
m_point2nc.F90\
m_prep_4_EnKF.F90\
m_put_mod_fld.F90\
m_set_random_seed2.F90\
m_uobs.F90\
nfw.F90\
EnKF.F90

ENKF_SRC_F77 = mod_raw_io.F

ENKF_SRC_C = order.c

ENKF_OBJ = $(ENKF_SRC_C:.c=.o) $(ENKF_SRC_F77:.F=.o) $(ENKF_SRC_F90:.F90=.o)

# some fine tuning; add more dependancies when/if required 
#
m_obs.o: m_uobs.o
m_Generate_element_Si.o: m_parse_blkdat.o mod_measurement.o m_get_mod_fld.o m_insitu.o m_obs.o
m_insitu.o: nfw.o mod_measurement.o
m_local_analysis.o: mod_measurement.o m_point2nc.o m_parameters.o

EnKF: $(ENKF_OBJ)
	@echo "->EnKF"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../EnKF $(ENKF_OBJ) $(LIBS) 

$(ENKF_OBJ): makefile make.inc MODEL.CPP

clean:
	@rm -f TMP/*.*  $(PROGS)

%.o: %.F90
	@echo "  $*".F90
	@rm -f ./TMP/$*.f90
	@cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS) > ./TMP/$*.f90
	@cd ./TMP ; $(CF90) -c $(FFLAGS) $(F90FLG) -o $*.o $*.f90

%.o: %.F
	@echo "  $*".F
	@rm -f ./TMP/$*.f
	@cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS)  > ./TMP/$*.f
	@cd ./TMP ; $(CF77) -c $(FFLAGS) $(F77FLG) -o $*.o $*.f  

%.o: %.c
	@echo "  $*".c
	@cd ./TMP ; $(CC) -c $(CFLAGS) -o $*.o -I.. ../$*.c
