! File:          m_obs.F90
!
! Created:       6 Feb 2008
!
! Last modified: 21 Feb 2008
!
! Author:        Pavel Sakov*
!                NERSC
!
! Purpose:       Generic code to deal with observations.
!
! Description:   This module contains the following functions and subroutines:
!                  - obs_setobs
!                      reads the observations into allocatable array obs(nobs)
!                      of type(measurement)
!                  - obs_prepareobs
!                      conducts state-dependent pre-processing of observations
!                  - obs_prepareuobs
!                      conducts state-dependent pre-processing of observations
!                      of a given type
!                It also contains the following data:
!                  - obs
!                      allocatable array of type(measurement)
!                  - nobs
!                      number of observations (may differ from the size of the
!                      array)
!
!                * This file contains some modified code of unknown origin
!                  from EnKF package. In particular, the code here supersedes
!                  the code from:
!                    m_get_nrobs_d.F90
!                    m_get_obs_d.F90
!
! Modifications:
!                09/11/2012 Geir Arne Waagbo:
!                -- Added support for OSISAF ice drift obs
!                29/07/2010 PS:
!                -- modified obs_QC(). The maximal increment now does not go to
!                   0 as the innovation increases, but rather is limited by 
!                   KMAX * sigma_ens
!                29/06/2010 PS:
!                 -- added obs_QC()
!                26/02/2008 PS: 
!                 -- put "obs" and "nobs" as public data in this module

module m_obs
#if defined (QMPI)
    use qmpi
#else
    use qmpi_fake
#endif
  use mod_measurement
  use m_uobs
  use m_insitu
  implicit none

  !
  ! public stuff
  !

  integer, public :: nobs = -1
  type(measurement), allocatable, dimension(:), public :: obs

  ! index of observation types used by the radius tuning of localization
  real,           allocatable, dimension(:), public :: Typobs

  public obs_readobs
  public obs_prepareobs
  public obs_QC

  !
  ! private stuff
  !

  private obs_testrange

  integer, parameter, private :: STRLEN = 512

  real, parameter, private :: TEM_MIN = -2.0d0
  real, parameter, private :: TEM_MAX = 50.0d0
  real, parameter, private :: SAL_MIN = 2.0d0
  real, parameter, private :: SAL_MAX = 40.0d0
  real, parameter, private :: SSH_MIN = -3.0d0
  real, parameter, private :: SSH_MAX = 3.0d0
  real, parameter, private :: ICEC_MIN = 0.0d0
  real, parameter, private :: ICEC_MAX = 0.996d0
  real, parameter, private :: UVICE_MIN = -100.0d0
  real, parameter, private :: UVICE_MAX = 100.0d0
  real, parameter, private :: HICE_MIN = 0.0d0
  real, parameter, private :: HICE_MAX = 5.99d0
  real, parameter, private :: SKIM_MIN = -2.0d5
  real, parameter, private :: SKIM_MAX = 2.0d5

  private obs_prepareuobs, obs_realloc

contains

  ! Obtain observations to be used for assimilation from the file
  ! "observation.uf". Store the number of observations in "nobs" and the data
  ! in the array "obs".
  !
  subroutine obs_readobs
    use m_parameters

    logical :: exists = .false.
    type(measurement) :: record
    integer :: rsize
    integer :: ios
    integer :: o

    if (nobs >= 0) then
       return
    end if

    inquire(file = 'observations.uf', exist = exists)
    if (.not. exists) then
       if (master) then
          print *, 'ERROR: obs_getnobs(): file "observations.uf" does not exist'
       end if
       stop
    end if
    inquire(iolength = rsize) record
    open(10, file = 'observations.uf', form = 'unformatted',&
         access = 'direct', recl = rsize, status = 'old')

    ! I guess there is no other way to work out the length other than read the
    ! file in fortran - PS
    !
    o = 1
    do while (.true.)
       read(10, rec = o, iostat = ios) record
       if (ios /= 0) then
          nobs = o - 1
          exit
       end if
       o = o + 1
    enddo

    allocate(obs(nobs),Typobs(nobs))
    Typobs=1

    ! PS - there were problem with using rewind(): g95 reported:
    ! "Cannot REWIND a file opened for DIRECT access". Therefore reopen.
    !
    close(10)
    open(10, file = 'observations.uf', form = 'unformatted',&
         access = 'direct', recl = rsize, status = 'old')
    do o = 1, nobs
       read(10, rec = o) obs(o)
       call ucase(obs(o) % id)
    enddo
    close(10)

    if (RFACTOR1 /= 1.0d0) then
       do o = 1, nobs
          obs(o) % var = obs(o) % var * RFACTOR1
       end do
    end if

   call  uobs_get(obs % id, nobs, Typobs, master)
   call obs_testrange
  end subroutine obs_readobs


  subroutine obs_testrange
    integer :: o, uo, nbad
    real :: dmin, dmax
       
    if (master) then
       print '(a)', ' EnKF: testing range for each type of obs '
    end if
    do uo = 1, nuobs
       if (trim(unique_obs(uo)) == 'SST' .or. trim(unique_obs(uo)) == 'TEM'&
            .or. trim(unique_obs(uo)) == 'GTEM') then
          dmin = TEM_MIN
          dmax = TEM_MAX
       elseif (trim(unique_obs(uo)) == 'SAL'&
            .or. trim(unique_obs(uo)) == 'SSS'&
            .or. trim(unique_obs(uo)) == 'GSAL') then
          dmin = SAL_MIN
          dmax = SAL_MAX
       elseif (trim(unique_obs(uo)) == 'SLA'&
            .or. trim(unique_obs(uo)) == 'TSLA'&
            .or. trim(unique_obs(uo)) == 'SSH') then
          dmin = SSH_MIN
          dmax = SSH_MAX
       elseif (trim(unique_obs(uo)) == 'ICEC') then
          dmin = ICEC_MIN
          dmax = ICEC_MAX
       elseif (trim(unique_obs(uo)) == 'HICE') then
          dmin = HICE_MIN
          dmax = HICE_MAX
       elseif (trim(unique_obs(uo)) == 'VICE'&
            .or. trim(unique_obs(uo)) == 'UICE') then
          dmin = UVICE_MIN
          dmax = UVICE_MAX
       elseif (trim(unique_obs(uo)) == 'SKIM') then
          dmin = SKIM_MIN
          dmax = SKIM_MAX
       elseif ((index(trim(unique_obs(uo)),'DX') .gt. 0) &
            .or. (index(trim(unique_obs(uo)),'DY') .gt. 0)) then
          ! The type can be DX1,DX2,..,DX5,DY1,..DY5
          dmin = UVICE_MIN
          dmax = UVICE_MAX
       else
          dmin = -1.0d6
          dmax = 1.0d6
          print *, 'ERROR: obs_testrange(): "', trim(unique_obs(uo)), '": unknown type'
          stop
       end if
       
       nbad = 0
       do o = uobs_begin(uo), uobs_end(uo)
          if (obs(o) % status .and.&
               (obs(o) % d < dmin .or. obs(o) % d > dmax)) then
             obs(o) % status = .false.
             nbad = nbad + 1
          end if
       end do
       if (master) then
          print '(a, a, a, i6, a)', '   ', trim(unique_obs(uo)), ': ', nbad, ' outliers'
       end if
    end do

    if (master) then
       print *
    end if
  end subroutine obs_testrange


  ! Prepare observations before allocating matrices S, D, and A in EnKF().
  ! This invloves mainly thinning, superobing, or sorting.
  !
  ! Note that generically this processing can not be completely outsourced
  ! to the preprocessing stage, at least for in-situ data, because its thinning
  ! involves reading ensemble members for layer depth information.
  !
  subroutine obs_prepareobs()
    implicit none

    integer :: iuobs

    if (master) then
       print '(a)', ' EnKF: preparing observations'
    end if
    do iuobs = 1, nuobs
       call obs_prepareuobs(trim(unique_obs(iuobs)))
    end do

   ! calculate again the number of observation of each type (that could change
   ! in prepare_obs)
    call  uobs_get(obs % id, nobs, Typobs, master)
  end subroutine obs_prepareobs


  ! Prepare (thin, superob) observations of type "obstag".
  !
  subroutine obs_prepareuobs(obstag)
    character(*), intent(in) :: obstag

    character(STRLEN) :: fname

    if (trim(obstag) == 'SAL' .or. trim(obstag) == 'TEM' .or.&
         trim(obstag) == 'GSAL' .or. trim(obstag) == 'GTEM') then
       call insitu_prepareobs(trim(obstag), nobs, obs)
       if (master) then
          write(fname, '(a, ".nc")') trim(obstag)
          print *, 'Writing "', trim(obstag), '" obs to be assimilated to "',&
               trim(fname), '"'
          call insitu_writeprofiles(fname, trim(obstag), nobs, obs);
       end if
    else
       ! do nothing for obs of other types for now
    end if
    call obs_realloc
  end subroutine obs_prepareuobs

  
  subroutine obs_realloc()
    type(measurement), allocatable :: newobs(:)
    real,           allocatable :: newTypobs(:)
    
    if (nobs < 0 .or. nobs == size(obs)) then
       return
    end if

    allocate(newobs(nobs))
    allocate(newTypobs(nobs))

    newobs = obs(1 : nobs)
    newTypobs = Typobs(1:nobs)

    deallocate(obs,Typobs)
    allocate(obs(nobs),Typobs(nobs))

    obs = newobs
    Typobs=newTypobs

    deallocate(newobs,newTypobs)
  end subroutine obs_realloc


  subroutine obs_QC(m, S)
    use m_parameters
    implicit none

    integer :: m
    real :: S(nobs, m)

    integer :: nmodified(nuobs)
    real :: so(m), smean, svar
    integer :: o, uo
    real :: ovar, inn, newovar

    if (master) then
       print *, 'Starting generic observation QC'
    end if

    nmodified = 0

    do uo = 1, nuobs
       do o = uobs_begin(uo), uobs_end(uo)
          so = S(o, :);
          smean = sum(so) / m ! must be 0...
          svar = sum((so - smean) ** 2) / real(m - 1)
          ovar = obs(o) % var

          inn = obs(o) % d - smean
          obs(o) % var = sqrt((svar + ovar) ** 2 +&
               svar * (inn / KFACTOR) ** 2) - svar

          if (svar > 0 .and. obs(o) % var / ovar > 2.0d0) then
             nmodified(uo) = nmodified(uo) + 1
          end if
       end do
    end do

    if (master) then
       do uo = 1, nuobs
          print *, '  ', trim(unique_obs(uo)), ':'
          print *, '    # of observations:', uobs_end(uo) - uobs_begin(uo) + 1
          print *, '    (of them) substantially modified:', nmodified(uo)
       end do
    end if
  end subroutine obs_QC

end module m_obs
