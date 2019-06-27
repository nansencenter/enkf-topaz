! File:          EnKF.F90
!
! Created:       ???
!
! Last modified: 20/04/2010
!
! Purpose:       Main program for EnKF analysis
!
! Description:   The workflow is as follows:
!                -- read model parameters
!                -- read obs
!                -- conduct necessary pre-processing of obs (superobing)
!                -- calculate ensemble observations
!                -- calculate X5
!                -- update the ensemble
!
! Modifications:
!                20/9/2011 PS:
!                  Modified code to allow individual inflations for each of
!                  `NFIELD' fields updated in a batch - thanks to Ehouarn Simon
!                  for spotting this inconsistency
!                6/8/2010 PS:
!                  Small changes in calls to calc_X5() and update_fields() to
!                  reflect changes in interfaces.
!                6/7/2010 PS:
!                  Moved point output to a separate module m_point2nc.F90
!                25/5/2010 PS:
!                  Added inflation as a 4th command line argument
!                20/5/2010 PS:
!                  Set NFIELD = 4. This requires 4 GB per node in TOPAZ and
!                  "medium" memory model on Hexagon (a single allocation for a
!                   variable over 2GB)
!                20/4/2010 PS:
!                  Set NFIELD = 4. This will require 2 GB per node in TOPAZ.
!                  Thanks to Alok Gupta for hinting this possibility.
!                10/4/2010 PS:
!                  Moved variable `field' from real(8) to real(4);
!                  set NFIELD = 2.
!                Prior history:
!                  Not documented.

program EnKF
#if defined(QMPI)
  use qmpi
#else
  use qmpi_fake
#endif
  use m_parameters
  use distribute
  use mod_measurement
  use m_get_mod_grid
  use m_get_mod_nrens
  use m_obs
  use m_local_analysis
  use m_prep_4_EnKF
  use m_set_random_seed2
  use m_get_mod_fld
  use m_put_mod_fld
  use mod_analysisfields
  use m_parse_blkdat
  use m_random
  use m_point2nc
  implicit none

  character(*), parameter :: ENKF_VERSION = "2.11"

  integer, external :: iargc

  ! NFIELD is the number of fields (x N) passed for the update during a call to
  ! update_fields(). In TOPAZ4 NFIELD = 2 if there is 1 GB of RAM per node, and
  ! NFIELD = 4 if there are 2 GB of RAM. Higher value of NFIELD reduces the
  ! number of times X5tmp.uf is read from disk, which is the main bottleneck
  ! for the analysis time right now.
  !
  integer, parameter :: NFIELD = 8

  character(512) :: options

  integer :: nrens
  real, allocatable, dimension(:,:) :: modlon, modlat, depths, readfld
  real, allocatable, dimension(:,:) :: S ! ensemble observations HE
  real, allocatable, dimension(:)   :: d ! d - Hx

  integer k, m

  ! "New" variables used in the parallelization 
  integer, dimension(:,:), allocatable :: nlobs_array
  real(4), allocatable :: fld(:,:)
  real(8) rtc, time0, time1, time2

  ! Additional fields
  character(len=3) :: cmem
  character(len=80) :: memfile
  integer :: fieldcounter

  character(100) :: text_string

  real :: rdummy
  integer :: idm, jdm, kdm

  real :: mindx
  real :: meandx
  integer :: m1, m2, nfields
  real :: infls(NFIELD)

#if defined(QMPI)
  call start_mpi()
#endif

  ! Read the characteristics of the assimilation to be carried out.

  if (iargc() /= 1) then
     print *, 'Usage: EnKF <parameter file>'
     print *, '       EnKF -h'
     print *, 'Options:'
     print *, '  -h -- describe parameter fie format'
     call stop_mpi()
  else
    call getarg(1, options)
    if (trim(options) == "-h") then
       call prm_describe()
       call stop_mpi()
    end if
  end if

  if (master) then
     print *
     print '(a, a)', ' EnKF version ', ENKF_VERSION
     print *
  end if

  call prm_read()
  call prm_print()

  ! get model dimensions
  !
   call parse_blkdat('idm   ', 'integer', rdummy, idm)
   call parse_blkdat('jdm   ', 'integer', rdummy, jdm)
   call parse_blkdat('kdm   ', 'integer', rdummy, kdm)

   allocate(modlon(idm, jdm))
   allocate(readfld(idm, jdm))
   allocate(modlat(idm, jdm))
   allocate(depths(idm, jdm))
   allocate(nlobs_array(idm, jdm))

   ! get model grid
   !
   call get_mod_grid(modlon, modlat, depths, mindx, meandx, idm, jdm)

   ! set a variable random seed
   !
!   call set_random_seed2
!   print *, 'Test Grid read'

   ! initialise point output
   !
   call p2nc_init

   time0 = rtc()

   ! read measurements
   !
   if (master) then
      print *, 'EnKF: reading observations'
   end if
   call obs_readobs
   if (master) then
      print '(a, i6)', '   # of obs = ', nobs
      print '(a, a, a, e10.3, a, e10.3)', '   first obs = "', trim(obs(1) % id),&
           '", v = ', obs(1) % d, ', var = ', obs(1) % var
      print '(a, a, a, e10.3, a, e10.3)', '   last obs = "', trim(obs(nobs) % id),&
           '", v = ', obs(nobs) % d, ', var = ', obs(nobs) % var
   end if
   if (master) then
      print *
   end if

   ! read ensemble size and store in A
   !
   nrens = get_mod_nrens(idm, jdm)
   if (master) then
      print '(a, i4, a)', ' EnKF: ', nrens, ' ensemble members found'
   end if
   if (ENSSIZE > 0) then
      ENSSIZE = min(nrens, ENSSIZE)
   else
      ENSSIZE = nrens
   end if
   if (master) then
      print '(a, i4, a)', ' EnKF: ', ENSSIZE, ' ensemble members used'
   end if
   if (master) then
      print *
   end if

   ! PS - preprocess the obs using the information about the ensemble fields
   ! here (if necessary), before running prep_4_EnKF(). This is necessary e.g.
   ! for assimilating in-situ data because of the dynamic vertical geometry in
   ! HYCOM
   !
   call obs_prepareobs

   allocate(S(nobs, ENSSIZE), d(nobs))
   call prep_4_EnKF(ENSSIZE, d, S, depths, meandx / 1000.0, idm, jdm, kdm)
   if (master) then
      print *, 'EnKF: finished initialisation, time = ',  rtc() - time0
   end if

   ! (no parallelization was required before this point)

   time1 = rtc()

   allocate(X5(ENSSIZE, ENSSIZE, idm))
   allocate(X5check(ENSSIZE, ENSSIZE, idm))
   call calc_X5(ENSSIZE, modlon, modlat, depths, mindx, meandx, d, S,&
        LOCRAD, RFACTOR2, nlobs_array, idm, jdm)
   deallocate(d, S, X5check)
   if (master) then
      print *, 'EnKF: finished calculation of X5, time = ', rtc() - time0
   end if

   allocate(fld(idm * jdm, ENSSIZE * NFIELD))

#if defined(QMPI)
   call barrier()
#endif

   ! get fieldnames and fieldlevels
   !
   call get_analysisfields()

   call distribute_iterations(numfields)
#if defined(QMPI)
   call barrier() !KAL - just for "niceness" of output
#endif
   time2 = rtc()
   do m1 = my_first_iteration, my_last_iteration, NFIELD
      m2 = min(my_last_iteration, m1 + NFIELD - 1)
      nfields = m2 - m1 + 1

      do m = m1, m2
         print '(a, i2, a, i3, a, a6, a, i3, a, f11.0)',&
              "I am ", qmpi_proc_num, ', m = ', m, ", field = ",&
              fieldnames(m), ", k = ", fieldlevel(m), ", time = ",&
              rtc() - time2
         do k = 1, ENSSIZE
            write(cmem, '(i3.3)') k
            memfile = 'forecast' // cmem
            call get_mod_fld_new(trim(memfile), readfld, k, fieldnames(m),&
                 fieldlevel(m), 1, idm, jdm)
            ! reshaping and conversion to real(4)
            fld(:, ENSSIZE * (m - m1) + k) = reshape(readfld, (/idm * jdm/))
         end do
         call p2nc_storeforecast(idm, jdm, ENSSIZE, numfields, m, fld(:, ENSSIZE * (m - m1) + 1 : ENSSIZE * (m + 1 - m1)))
         infls(m - m1 + 1) = prm_getinfl(trim(fieldnames(m)));
      end do

      call update_fields(idm, jdm, ENSSIZE, nfields, nlobs_array, depths,&
              fld(1,1), infls)

      do m = m1, m2
         fieldcounter = (m - my_first_iteration) + 1
         do k = 1, ENSSIZE
            write(cmem,'(i3.3)') k
            memfile = 'analysis' // cmem
            ! reshaping and conversion to real(8)
            readfld = reshape(fld(:, ENSSIZE * (m - m1) + k), (/idm, jdm/))
            write(text_string, '(a, i3.3)') '_proc', qmpi_proc_num
            call put_mod_fld(trim(memfile) // trim(text_string), readfld, k,&
                 fieldnames(m), fieldlevel(m), 1, fieldcounter, idm, jdm)
         end do
      end do
   end do
   deallocate(X5)
   deallocate(fld)

   call p2nc_writeforecast

   ! Barrier only necessary for timings
#if defined(QMPI)
   call barrier()
#endif
   if (master) then
      print *, 'EnKF: time for initialization = ', time1 - time0
      print *, 'EnKF: time for X5 calculation = ', time2 - time1
      print *, 'EnKF: time for ensemble update = ', rtc() - time2
      print *, 'EnKF: total time = ', rtc() - time0
   end if
   print *, 'EnKF: Finished'
   call stop_mpi()
 end program EnKF

#if defined(_G95_)
 ! not tested! - PS
 !
 real function rtc()
   integer :: c

   call system_clock(count=c)
   rtc = dfloat(c)
 end function rtc
#endif
