










!Fanf A short program to ensure that all member have the ice present.
!
! PS: if there are corrupted ic fields - then
!       (i) the member IDs of these fields will be wrtten to the file
!           "missing_icerecords.txt"
!       (ii) this also will be reported to stdout
!       (iii) and also will be reflected in icevolume.txt

program checkice_en
  use mod_raw_io
  use m_parse_blkdat
  use m_get_mod_grid
  use m_get_mod_fld
  implicit none
  integer*4, external :: iargc
  integer iens
  real, dimension(:,:), allocatable :: modlon,modlat,depths
  logical, allocatable, dimension(:, :) :: iswater
  integer :: idm,jdm,kdm
  integer :: ios
  integer ::  nens
  real, allocatable, dimension(:,:) :: ficem,hicem,hsnwm,ticem,tsrfm
  real, allocatable, dimension(:,:) :: fld
  integer :: reclICE
  real :: mindx,meandx,rdummy
  character(len=80) :: icerestart,filename
  character(len=3) :: ctmp
  integer :: nmissing
  real, allocatable, dimension(:) :: icevolume, icearea
  real :: meanicevolume, maxvalue_hicem, maxvalue_ticem, maxvalue_tsrfm
  logical :: ex

  if ( iargc()==2 ) then
     call getarg(1,icerestart)
     call getarg(2,ctmp)
     read(ctmp,*) nens
  else
     print *,'"check_ice_en" -- A routine to check that no ice records are missing'
     print *
     print *,'Usage: checkice_en <surname of ice_file> <ensemble_size>'
     call exit(1)
  endif

  !Get model dimensions
  call parse_blkdat('idm   ','integer',rdummy,idm)
  call parse_blkdat('jdm   ','integer',rdummy,jdm)

  open(11, file = 'icevolume.txt', status = 'replace')
  close(11)

  allocate(modlon (idm,jdm))
  allocate(modlat (idm,jdm))
  allocate(depths (idm,jdm))
  allocate(fld (idm,jdm))
  
  call get_mod_grid(modlon, modlat, depths, mindx, meandx, idm, jdm)

  allocate(iswater(idm, jdm))
  iswater = depths > 1.0d0 .and. depths < 1.0e25

  allocate(icevolume(nens))
  allocate(icearea(nens))
  icevolume = 0.0d0
  icearea = 0.0d0

  allocate(ficem(idm,jdm))
  allocate(hicem(idm,jdm))
  allocate(hsnwm(idm,jdm))
  allocate(ticem(idm,jdm))
  allocate(tsrfm(idm,jdm))
!  inquire(iolength = reclICE) ficem, hicem, hsnwm, ticem, tsrfm
   
!  open(20, file = trim(icerestart), form = 'unformatted', access = 'direct',&
!       recl = reclICE, status = 'old', iostat = ios)
!  if (ios /= 0) then
!     print *, 'ERROR: problem reading "', trim(icerestart), '"'
!     call exit(1)
!  end if
! to check the existing for the mem-file
  do iens=1,nens
    write(ctmp,'(i3.3)') iens  
    filename=trim(icerestart)//trim(ctmp)
    print *,iens
    inquire(exist=ex,file=trim(filename)//'.b')
    if (.not.ex) then
      write(*,*) 'Can not find '//trim(filename)//'.b'
      stop '(EnKF_postprocess)'
    end if
    call get_mod_fld_new(trim(filename),fld(:,:),1,'ficem   ',0,1,idm,jdm); ficem=fld;
    call get_mod_fld_new(trim(filename),fld(:,:),1,'hicem   ',0,1,idm,jdm); hicem=fld;
    icevolume(iens) = sum(ficem * hicem, mask = iswater)
    icearea(iens) = sum(ficem, mask = iswater)
  end do
  open(11, file = 'icevolume.txt', status = 'old', position = 'append')
    do iens=1,nens
      write(11, '(i4, f14.0, f14.0)') iens, icevolume(iens), icearea(iens)
    end do
  close(11)
  meanicevolume = sum(icevolume) / real(nens)

  nmissing = 0
  do iens=1,nens
    write(ctmp,'(i3.3)') iens  
    filename=trim(icerestart)//trim(ctmp)
    call get_mod_fld_new(trim(filename),fld(:,:),1,'hicem   ',0,1,idm,jdm); hicem=fld;
    call get_mod_fld_new(trim(filename),fld(:,:),1,'ticem   ',0,1,idm,jdm); ticem=fld;
    call get_mod_fld_new(trim(filename),fld(:,:),1,'tsrfm   ',0,1,idm,jdm); tsrfm=fld;
     maxvalue_hicem = maxval(hicem, mask = iswater) ! In meters
     maxvalue_ticem = maxval(ticem, mask = iswater) ! In Kelvin 
     maxvalue_tsrfm = maxval(tsrfm, mask = iswater) ! In Kelvin 
     if (maxvalue_hicem < 0.1  .or. maxvalue_hicem > 100.0 .or. &
         maxvalue_ticem < 10.0 .or. maxvalue_tsrfm < 10.0) then
        nmissing = nmissing + 1
        print '(A, $)', '-'
        open(10, file = 'missing_icerecords.txt', position = 'append')
        write(10, '(i4)') iens
        close(10)
     elseif (icevolume(iens) /= icevolume(iens) .or. (meanicevolume - icevolume(iens)) / meanicevolume > 0.35) then
        nmissing = nmissing + 1
        print '(A, $)', '*'
        print *, 'member ', iens, ': icevolume = ', icevolume(iens),&
             ', meanicevolume = ', meanicevolume
        open(10, file = 'missing_icerecords.txt', position = 'append')
        write(10, '(i4)') iens
        close(10)
     else
        print '(A, $)', '.'
     end if
  end do
! close(20)
  print *, ''
  if (nmissing > 0) then
     print *, 'ERROR: ice field is missing for', nmissing, ' member(s)',&
          ' check "missing_icerecords.txt" for member IDs'
  end if

end program checkice_en
