MPI = NO
include ../make.inc

VPATH = ../:./:TMP

PROGS = EnKF_assemble\
checkice\
checkice_en\
consistency\
fixhycom\
fixhycom_eco\
obsstats\
oldtonewobs

PROGS0 = consistency

all: $(PROGS)

ENKF_ASSEMBLE_OBJ = qmpi.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_put_mod_fld.o\
p_EnKF_assemble.o

CONSISTENCY_OBJ = qmpi.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_put_mod_fld.o\
m_get_mod_grid.o\
m_get_mod_fld.o\
nfw.o\
mod_testinfo.o\
p_consistency.o

CHECK_ICE_OBJ = qmpi.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_get_mod_grid.o\
p_check_ice.o

CHECK_ICE_OBJ_EN = qmpi.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_get_mod_grid.o\
m_get_mod_fld.o\
p_check_ice_en.o


FIXHYCOM_OBJ = qmpi.o\
m_random.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_put_mod_fld.o\
m_get_mod_grid.o\
m_get_mod_fld.o\
p_fixhycom.o

FIXHYCOM_ECO_OBJ = qmpi.o\
m_fixhycom_eco_metno.o\
mod_raw_io.o\
m_parse_blkdat.o\
m_put_mod_fld.o\
m_get_mod_grid.o\
m_get_mod_fld.o\
distribute.o\
p_fixhycom_eco.o

OBSSTATS_SRC_C = order.c

OBSSTATS_SRC_F90 = qmpi.F90\
m_parameters.F90\
mod_measurement.F90\
mod_raw_io.F90\
m_parse_blkdat.F90\
m_put_mod_fld.F90\
m_get_mod_grid.F90\
m_spherdist.F90\
m_Generate_element_Si.F90\
m_get_mod_nrens.F90\
m_insitu.F90\
m_obs.F90\
m_set_random_seed2.F90\
mod_sphere_tools.F90\
m_get_mod_fld.F90\
m_random.F90\
m_prep_4_EnKF.F90\
m_uobs.F90\
nfw.F90\
p_obsstats.F90

TESTX4_SRC_F90 = qmpi.F90\
m_parse_blkdat.F90\
m_get_mod_nrens.F90\
testX4.F90

OBSSTATS_OBJ = $(OBSSTATS_SRC_C:.c=.o) $(OBSSTATS_SRC_F77:.F=.o) $(OBSSTATS_SRC_F90:.F90=.o)

OLDTONEWOBS_OBJ = mod_measurement_oldnew.o\
p_oldtonewobs.o

TESTX4_OBJ = $(TESTX4_SRC_F90:.F90=.o)

m_Generate_element_Si.o: m_parse_blkdat.o mod_measurement.o m_get_mod_fld.o m_insitu.o m_obs.o
m_insitu.o: nfw.o
m_get_mod_nrens.o: qmpi.o
m_obs.o : m_parameters.o m_uobs.o

$(ENKF_ASSEMBLE_OBJ) $(CONSISTENCY_OBJ) $(FIXHYCOM_ECO_OBJ) $(FIXHYCOM_OBJ) $(OBSSTATS_OBJ) $(OLDTONEWOBS_OBJ): MODEL.CPP makefile

EnKF_assemble: $(ENKF_ASSEMBLE_OBJ)
	@echo "->EnKF_assemble"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../EnKF_assemble $(ENKF_ASSEMBLE_OBJ) 

consistency: $(CONSISTENCY_OBJ)
	@echo "->consistency"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../consistency $(CONSISTENCY_OBJ) $(LIBS)

checkice: $(CHECK_ICE_OBJ)
	   @echo "->checkice"
	   @cd ./TMP ; $(LD) $(LINKFLAGS) -o ../checkice $(CHECK_ICE_OBJ) $(LIBS)

checkice_en: $(CHECK_ICE_OBJ_EN)
	   @echo "->checkice_en"
	   @cd ./TMP ; $(LD) $(LINKFLAGS) -o ../checkice_en $(CHECK_ICE_OBJ_EN) $(LIBS)

fixhycom: $(FIXHYCOM_OBJ)
	@echo "->fixhycom"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../fixhycom $(FIXHYCOM_OBJ)

fixhycom_eco: $(FIXHYCOM_ECO_OBJ)
	@echo "->fixhycom_eco"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../fixhycom_eco $(FIXHYCOM_ECO_OBJ)

obsstats: $(OBSSTATS_OBJ)
	@echo "->obsstats"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../obsstats $(OBSSTATS_OBJ) $(LIBS)  $(MATLIBS)

oldtonewobs: $(OLDTONEWOBS_OBJ)  
	@echo "->oldtonewobs"
	@cd ./TMP ; $(LD) $(LINKFLAGS) -o ../oldtonewobs $(OLDTONEWOBS_OBJ) 

clean:
	@rm -f *.o *.mod TMP/*.* $(PROGS)

%.o: %.F90
	@echo "  $*".F90
	@rm -f ./TMP/$*.f90
	@cat MODEL.CPP $< | $(CPP) $(CPPFLAGS) > ./TMP/$*.f90
	@cd ./TMP; $(CF90) -c $(FFLAGS) $(F90FLG) -o $*.o $*.f90

%.o: %.F
	@echo "  $*".F
	@rm -f ./TMP/$*.f
	@cat MODEL.CPP $< | $(CPP) $(CPPFLAGS) > ./TMP/$*.f
	@cd ./TMP; $(CF77) -c $(FFLAGS) $(F77FLG) -o $*.o $*.f 2> /dev/null

%.o: %.c
	@echo "  $*".c
	@cd ./TMP ; $(CC) -c $(CFLAGS) -o $*.o -I../.. ../../$*.c
