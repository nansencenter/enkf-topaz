! File:          m_uobs.F90
!
! Created:       11 August 2010
!
! Last modified: 11.8.2010
!
! Author:        Pavel Sakov
!                NERSC
!
! Purpose:       Handle different observation types.
!
! Description:   This module is in charge of sorting of observations by types
!                and storing the results
!
! Modifications: None

module m_uobs
#if defined (QMPI)
    use qmpi
#else
    use qmpi_fake
#endif
  use mod_measurement
  implicit none

  public uobs_get
  
  integer, parameter, private :: MAXNUOBS = 19

  integer, public :: nuobs
  character(OBSTYPESTRLEN), public :: unique_obs(MAXNUOBS)
  integer, public :: nobseach(MAXNUOBS)
  integer :: uobs_begin(MAXNUOBS), uobs_end(MAXNUOBS)

contains

  subroutine uobs_get(tags, nrobs, Trobs,master)
    implicit none
    integer , intent(in) :: nrobs
    logical , intent(in) :: master
    character(OBSTYPESTRLEN), intent(in) :: tags(nrobs)
    real,intent(out) :: Trobs(nrobs)       ! used for Typobs by localization

    logical :: obsmatch
    integer :: o, uo

    nobseach = 0

    ! check for unique obs
    if (master) then
       print '(a)', ' EnKF: getting unique observations '
    end if
    nuobs = 0
    unique_obs = ''
    Trobs=1.0
    do o = 1, nrobs
       obsmatch = .false.
#if defined(PROF_RADIUS)
       if (trim(tags(o)) == 'SAL' .or. trim(tags(o)) == 'TEM' .or.&
          trim(tags(o)) == 'GSAL' .or. trim(tags(o)) == 'GTEM') then
          Trobs(o)=2.0
  !     else
  !        Trobs(o)=1
       endif 
#endif
       do uo = 1, nuobs
          if (trim(tags(o)) == trim(unique_obs(uo))) then
             obsmatch = .true.
             nobseach(uo) = nobseach(uo) + 1
             exit
          end if
       end do
       if (.not. obsmatch) then
          nuobs = nuobs + 1
          nobseach(nuobs) = 1
          if (nuobs > MAXNUOBS) then
             if (master) then
                print *, 'ERROR: uobs_get(): # of unique obs = ', nuobs,&
                     ' > MAXNUOBS = ', MAXNUOBS
                print *, '  obs # = ', o, ', tag = ', trim(tags(o))
             end if
             stop
          end if
          unique_obs(nuobs) = trim(tags(o))
       end if
    end do
    if (master) then
       do uo = 1, nuobs
          print '(a, i2, a, a, a, i7, a)', '   obs variable  ', uo, ' -- ',&
               trim(unique_obs(uo)), ',', nobseach(uo), ' observations'
       end do
    end if
    uobs_begin(1) = 1
    uobs_end(1) = nobseach(1)
    do uo = 2, nuobs
       uobs_begin(uo) = uobs_end(uo - 1) + 1
       uobs_end(uo) = uobs_begin(uo) + nobseach(uo) - 1
    end do
    if (master) then
       do uo = 1, nuobs
          do o = uobs_begin(uo), uobs_end(uo)
             if (trim(tags(o)) /= trim(unique_obs(uo))) then
                print *, 'ERROR: uobs_get(): uinique observations not ',&
                     'continuous in observation array'
                stop
             end if
          end do
       end do
    end if
    if (master) then
       print *
    end if
  end subroutine uobs_get

end module m_uobs
