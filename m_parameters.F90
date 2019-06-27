! File:          m_parameters.F90
!
! Created:       6 August 2010
!
! Last modified: 6/8/2010
!
! Author:        Pavel Sakov
!                NERSC
!
! Purpose:       Provide a simpl nml list-based parameter input into EnKF.
!
! Description:   Provides code for reading parameters from a specified 
!                parameter file.
!
! Modifications: none

module m_parameters
#if defined(QMPI)
  use qmpi
#else
  use qmpi_fake
#endif
  implicit none

  integer, parameter, private :: STRLEN = 512
  integer, parameter, private :: FID = 101

  character(STRLEN), public :: PRMFNAME = "NONE"

  integer, public :: ENSSIZE = 0
  namelist /ensemble/ ENSSIZE

  character(STRLEN), public :: METHODTAG = "NONE"
  namelist /method/ METHODTAG

  real, public :: LOCRAD = 0.0d0
  character(STRLEN), public :: LOCFUNTAG = "GASPARI-COHN"
  namelist /localisation/ LOCRAD, LOCFUNTAG

  real, public :: INFL = 1.0d0
  real, public :: RFACTOR1 = 1.0d0
  real, public :: RFACTOR2 = 1.0d0
  real, public :: KFACTOR = 2.0d0
  namelist /moderation/ INFL, RFACTOR1, RFACTOR2, KFACTOR

  character(STRLEN), public :: JMAPFNAME = "NONE"
  character(STRLEN), public :: POINTFNAME = "NONE"
  character(STRLEN), public :: MEANSSHFNAME = "NONE"
  namelist /files/ JMAPFNAME, POINTFNAME, MEANSSHFNAME

  integer, parameter, private :: NPRMESTMAX = 10
  integer :: nprmest = 0
  character(STRLEN), dimension(NPRMESTMAX), public :: PRMESTNAME
  real, dimension(NPRMESTMAX), public :: PRMINFL
  namelist /prmest/ PRMESTNAME, PRMINFL

  public prm_read, prm_describe, prm_print, prm_getinfl, prm_prmestexists, ucase

contains

  subroutine prm_read
    integer :: ios, i

    call getarg(1, PRMFNAME)

    if (master) then
       print *, 'EnKF: reading parameters from "', trim(PRMFNAME), '":'
    end if

    open(unit = FID, file = trim(PRMFNAME), form = "formatted",&
         status = "old", iostat = ios)
    if (ios /= 0) then
       if (master) then
          print *,  'ERROR: could not open "', trim(PRMFNAME), '", iostatus =', ios
          stop
       end if
    end if

    read(unit = FID, nml = method, iostat = ios)
    if (ios /= 0) then
       if (master) then
          print *, 'ERROR: "', trim(PRMFNAME), '": could not read namelist "method"'
       end if
       stop
    end if
    rewind(FID)

    read(unit = FID, nml = ensemble, iostat = ios)
    if (ios /= 0) then
       if (master) then
          print *, 'ERROR: "', trim(PRMFNAME), '": could not read namelist "ensemble"'
       end if
       stop
    end if
    rewind(FID)

    read(unit = FID, nml = localisation, iostat = ios)
    if (ios /= 0) then
       if (master) then
          print *, 'ERROR: "', trim(PRMFNAME), '": could not read namelist "localisation"'
       end if
       stop
    end if
    rewind(FID)

    read(unit = FID, nml = moderation, iostat = ios)
    if (ios /= 0) then
       if (master) then
          print *, 'ERROR: "', trim(PRMFNAME), '": could not read namelist "moderation"'
       end if
       stop
    end if
    rewind(FID)

    read(unit = FID, nml = files, iostat = ios)
    if (ios /= 0) then
       if (master) then
          print *, 'ERROR: "', trim(PRMFNAME), '": could not read namelist "files"'
       end if
       stop
    end if
    rewind(FID)

    do i = 1, NPRMESTMAX
       PRMESTNAME(i) =  ""
    end do
    read(unit = FID, nml = prmest, iostat = ios)
    if (ios /= 0) then
       if (master) then
          print *, 'ERROR: "', trim(PRMFNAME), '": could not read namelist "prmest"'
       end if
       stop
    end if
    do i = 1, NPRMESTMAX
       if (PRMESTNAME(i) ==  "") then
          nprmest = i - 1
          exit
       end if
    end do
    rewind(FID)

    close(FID)

    call ucase(METHODTAG)
    call ucase(LOCFUNTAG)
  end subroutine prm_read


  subroutine prm_describe
    if (.not. master) then
       return
    end if

    print '(a)', ' Example of EnKF parameter file:'
    print *
    print '(a)', '&method'
    print '(a)', '     methodtag    = "DEnKF"'
    print '(a)', '/'
    print '(a)', '&ensemble'
    print '(a)', '     enssize      = 0'
    print '(a)', '/'
    print '(a)', '&localisation'
    print '(a)', '     locfuntag    = "Gaspari-Cohn"'
    print '(a)', '     locrad       = 300.0'
    print '(a)', '/'
    print '(a)', '&moderation'
    print '(a)', '     infl         = 1.01 (<number>)'
    print '(a)', '     rfactor1     = 1.0 (<number>)'
    print '(a)', '     rfactor2     = 2.0 (<number>)'
    print '(a)', '     kfactor      = 2.0 (<number>)'
    print '(a)', '/'
    print '(a)', '&files'
    print '(a)', '     jmapfname    = "jmap.txt" (<file name>)'
    print '(a)', '     pointfname   = "point2nc.txt" (<file name>)'
    print '(a)', '     meansshfname = "meanssh.uf" (<file name>)'
    print *
    print '(a)', 'Parameter options:'
    print '(a)', '  method          = "EnKF" | "DEnKF"*'
    print '(a)', '  enssize         = <number> (0* to use all available states)'
    print '(a)', '  locfuntag       = "Gaspari-Cohn"* | "Step" | "None"'
    print '(a)', '  locrad          = <support radius in km>'
    print '(a)', '  infl            = <multiple, for ensemble anomalies> (* 1.0)'
    print '(a)', '  rfactor1        = <obs. error variance multiple> (* 1.0)'
    print '(a)', '  rfactor2        = <additional multiple for updating ens. anomalies> (* 1.0)'
    print '(a)', '  kfactor         = <max. allowed increment in terms of ensemble spread> (* 2.0)'
    print '(a)', '  jmapfname*      = <file with j remapping> (* none)'
    print '(a)', '  pointfname*     = <file with point coordinates> (* none)'
    print '(a)', '  meansshfname*   = <file with mean SSH> (* none)'
  end subroutine prm_describe


  subroutine prm_print
    integer :: i

    if (.not. master) then
       return
    end if

    print '(a)', ' EnKF parameters:'
    print '(a)', '   method:'
    print '(a, a, a)',  '     methodtag   = "', trim(METHODTAG), '"'
    print '(a)', '   ensemble:'
    print '(a, i0)',    '     enssize     = ', ENSSIZE
    print '(a)', '   localisation:'
    print '(a, f5.0)',  '     locrad      = ', LOCRAD
    print '(a, a, a)',  '     locfuntag   = "', trim(LOCFUNTAG), '"'
    print '(a)', '   moderation:'
    print '(a, f5.3)',  '     infl        = ', INFL
    print '(a, f3.1)',  '     rfactor1    = ', RFACTOR1
    print '(a, f3.1)',  '     rfactor2    = ', RFACTOR2
    print '(a, f3.1)',  '     kfactor     = ', KFACTOR
    print '(a)', '   files:'
    print '(a, a, a)', '     jmapfname    = "', trim(JMAPFNAME), '"'
    print '(a, a, a)', '     pointfname   = "', trim(POINTFNAME), '"'
    print '(a, a, a)', '     meansshfname = "', trim(MEANSSHFNAME), '"'
    print '(a, i0, a)', '   prmest: ', nprmest, ' fields'
    do i = 1, nprmest
       print '(a, a, a, f5.3)', '     prmestname = "', trim(PRMESTNAME(i)), '", infl = ', PRMINFL(i)
    end do
    print *
  end subroutine prm_print


  function prm_getinfl(fldname)
    real :: prm_getinfl
    character(*), intent(in) :: fldname
    integer :: i
    
    prm_getinfl = INFL
    do i = 1, nprmest
       if (trim(fldname) == PRMESTNAME(i)) then
          prm_getinfl = PRMINFL(i)
          print '(a, a, a, f5.3)', ' "', trim(fldname), '": using inflation = ', prm_getinfl
          return
       end if
    end do
  end function prm_getinfl


  function prm_prmestexists(varname)
    logical :: prm_prmestexists
    character(*), intent(in) :: varname
    integer :: i
    
    prm_prmestexists = .false.
    do i = 1, nprmest
       if (trim(varname) == PRMESTNAME(i)) then
          prm_prmestexists = .true.
          return
       end if
    end do
  end function prm_prmestexists


  ! Shift a character string to upper case.
  !
  subroutine ucase(string)
    character(*) :: string
    integer :: i

    do i = 1, len(string)
       if (string(i:i) >= 'a' .and. string(i:i) <= 'z') then
          string(i:i) = achar (iachar ( string(i:i) ) - 32)
       end if
    end do
  end subroutine ucase

end module m_parameters
