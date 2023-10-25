MPI = NO
include ./make.inc

SHELL = /bin/bash
VPATH = ..:.:TMP

PROGS=prep_obs_hice
CHK_OBS=chk_obs

all: $(PROGS) $(CHK_OBS)

PREP_SRC_F90 = mod_angles.F90\
byteswapper.F90\
m_confmap.F90\
m_bilincoeff.F90\
m_get_def_wet_point.F90\
m_nf90_err.F90\
m_oldtonew.F90\
m_pivotp.F90\
m_read_CERSAT_data.F90\
m_read_OSISAF_data.F90\
m_read_CLS_SLA.F90\
m_read_CLS_TSLA.F90\
m_read_MET_SST.F90\
m_read_CLS_SSH.F90\
m_read_CLS_SST.F90\
m_read_CLS_SST_grid.F90\
m_read_MET_SST_grid.F90\
m_read_CLS_TSLA_grid.F90\
m_read_CLS_data.F90\
m_read_CLS_header.F90\
m_read_FFI_glider.F90\
m_read_ifremer_argo.F90\
m_read_ifremer_skim.F90\
m_read_amsr_norsex.F90\
m_read_metno_icec.F90\
m_read_SMOS_HICE.F90\
m_read_jpl_hice.F90\
m_superobs.F90\
m_uobs.F90\
m_write_wet_file.F90\
mod_grid.F90\
nfw.F90\
qmpi.F90

PREP_OBS_SRC_F90 = $(PREP_SRC_F90)\
mod_measurement.F90\
m_spherdist.F90\
m_parse_blkdat.F90\
m_get_mod_grid.F90\
p_prep_obs.F90

CHK_OBS_SRC_F90 = $(CHK_SRC_F90)\
mod_measurement.F90\
p_chk_obs.F90

PREP_OBS_SRC_F = mod_raw_io.F

PREP_OBS_SRC_C = superobs.c\
superobs3d.c

PREP_OBS_OBJ = $(PREP_OBS_SRC_C:.c=.o) $(PREP_OBS_SRC_F90:.F90=.o) $(PREP_OBS_SRC_F:.F=.o)
CHK_OBS_OBJ = $(CHK_OBS_SRC_F90:.F90=.o) $(PREP_OBS_SRC_C:.c=.o) $(PREP_OBS_SRC_F:.F=.o) 

m_bilincoeff.o: m_oldtonew.o
m_get_def_wet_point.o: m_pivotp.o m_confmap.o m_spherdist.o mod_measurement.o mod_grid.o
m_get_mod_grid.o: mod_raw_io.o
m_read_CERSAT_data.o: nfw.o
m_read_CLS_TSLA.o: nfw.o
m_read_MYO_TSLA.o: nfw.o
p_prep_obs.o: m_uobs.o
m_uobs.o: qmpi.o

chk_obs: $(CHK_OBS_OBJ)
	@echo "->chk_obs"
	@cd ./TMP; $(LD) $(LINKFLAGS) -o ../chk_obs $(CHK_OBS_OBJ) $(LIBS)

prep_obs_hice: $(PREP_OBS_OBJ)
	@echo "->prep_obs_hice"
	@cd ./TMP; $(LD) $(LINKFLAGS) -o ../prep_obs_hice $(PREP_OBS_OBJ) $(LIBS)

$(PREP_OBS_OBJ): makefile make.inc MODEL.CPP

clean:
	@rm -f TMP/*.f  TMP/*.o TMP/*.f90 TMP/*.h TMP/*.mod $(PROGS)

%.o: %.F90
	@echo "  $*".F90
	@rm -f ./TMP/$*.f90
	@cat MODEL.CPP $< | $(CPP) $(CPPFLAGS)  > ./TMP/$*.f90
	@cd ./TMP; $(CF90) -c $(FFLAGS) $(F90FLG) -o $*.o $*.f90

%.o: %.F
	@echo "  $*".F
	@rm -f ./TMP/$*.f
	@cat MODEL.CPP $< | $(CPP) $(CPPFLAGS) > ./TMP/$*.f
	@cd ./TMP; $(CF77) -c $(FFLAGS) $(F77FLG) -o $*.o $*.f 2> /dev/null

%.o: %.c
	@echo "  $*".c
	@cd ./TMP ; $(CC) -c $(CFLAGS) -o $*.o -I.. -I../.. ../$*.c
