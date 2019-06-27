










! File:          m_read_FFI_glider.F90
!
! Created:       November 2009
!
! Author:        Pavel Sakov
!                NERSC
!
! Purpose:       Read glider data from text files by FFI into TOPAZ system
!
! Description:   Data file(s) are defined by the string in the 4th line of
!                "infile.data". It should have the following format:
!                <BEGIN>
!                FFI
!                GSAL | GTEM
!                <obs. error variance>
!                <File name>
!                <END>
!
!                This is a very beta code, just to make an initial assessment.
!
! Modifications: none

module m_read_FFI_glider
  implicit none

  integer, parameter, private :: STRLEN = 512

  public read_ffi_glider

  private grid_readxyz

contains

  subroutine read_ffi_glider(fname, obstype, variance, nx, ny, data)
    use mod_measurement
    use m_confmap
    use m_oldtonew
    use m_pivotp
    use m_bilincoeff

    character(*), intent(in) :: fname
    character(*), intent(in) :: obstype
    real, intent(in) :: variance
    integer, intent(in) :: nx, ny
    type(measurement), allocatable, intent(out) :: data(:)

    real*8, dimension(nx, ny) :: modlat, modlon, depths
    real :: latnew, lonnew

    character(STRLEN) :: record
    integer :: ios
    integer :: r, nr, o, nobs, oo

    real :: tmp
    type(measurement) :: obs
    type(measurement), allocatable :: tmpdata(:)

    ! count number of records
    !
    open(10, file = trim(fname), access = 'sequential', status = 'old', iostat = ios)
    if (ios /= 0) then
       print *, 'ERROR: read_FFI_glider(): could not open "', fname, '"'
    end if
    nr = 1
    do while(.true.)
       read(10, *, iostat = ios) record
       if (ios /= 0) then
          exit
       end if
       nr = nr + 1
    end do

    print *, trim(fname), ': ', nr, ' lines'
    if (nr == 0) then
       print *, 'ERROR: read_FFI_glider(): "', fname, '": empty file?'
       stop
    end if

    allocate(data(nr))

    close(10)
    open(10, file = trim(fname), access = 'sequential', status = 'old')
    nobs = 0
    do r = 1, nr
       if (trim(obstype) == 'GSAL' .or. trim(obstype) == 'SAL') then
          read(10, *, iostat = ios) obs % date, obs % lat, obs % lon, obs % depth, tmp, tmp, obs % d
       elseif (trim(obstype) == 'GTEM' .or. trim(obstype) == 'TEM') then
          read(10, *, iostat = ios) obs % date, obs % lat, obs % lon, obs % depth, tmp, obs % d, tmp
       else
          print *, trim(fname), ': unknown data type "', trim(obstype), '"'
          stop
       end if
       if (obs % date <= 0) then
          cycle
       end if
       nobs = nobs + 1
       data(nobs) = obs
    end do
    close(10)

    allocate(tmpdata(1 : nobs))
    tmpdata = data(1 : nobs)
    deallocate(data)
    allocate(data(nobs))
    data = tmpdata
    deallocate(tmpdata)

    if (nobs == 0) then
       print *, 'ERROR: read_FFI_glider(): "', trim(fname),&
            '": no meaningful data for ', trim(obstype), ' found'
       stop
    end if
    print *, trim(fname), ': ', nobs, ' records for ', trim(obstype)

    data % id = obstype
    data % var = variance
    data % status = .true.
    data % ns = 0
    data % i_orig_grid = 0
    ! convert seconds since 1/1/1970 to days since 1/1/1950
    !
    ! data(1 : nobs) % date = data(1 : nobs) % date / 86400 + 7305

    call confmap_init(nx, ny)
    call grid_readxyz(nx, ny, modlat, modlon, depths)
    do o = 1, nobs
       call oldtonew(data(o) % lat, data(o) % lon, latnew, lonnew)
       call pivotp(lonnew, latnew, data(o) % ipiv, data(o) % jpiv)
       if (data(o) % ipiv < 1 .or. data(o) % jpiv < 1&
            .or. data(o) % ipiv > nx - 1 .or. data(o) % jpiv > ny - 1) then
          data(o) % status = .false.
       else
          call bilincoeff(real(modlon), real(modlat), nx, ny, data(o) % lon, data(o) % lat,&
               data(o) % ipiv, data(o) % jpiv, data(o) % a1, data(o) % a2,&
               data(o) % a3, data(o) % a4)
       end if
    end do

    ! some basic QC
    where (data % depth < 0.0d0 .or. data % depth > 6000.0d0)
       data % status = .false.
    end where
    if (trim(obstype) == 'TEM') then
       where (data % d < -3.0d0 .or. data % d > 40.0d0)
          data % status = .false.
       end where
    elseif (trim(obstype) == 'SAL') then
       where (data % d < 30.0d0 .or. data % d > 40.0d0)
          data % status = .false.
       end where
    end if

    allocate(tmpdata(1 : count(data % status)))
    oo = 0
    do o = 1, nobs
       if (data(o) % status) then
          oo = oo + 1
          tmpdata(oo) = data(o)
       end if
    end do
    nobs = oo
    deallocate(data)
    allocate(data(nobs))
    data = tmpdata
    deallocate(tmpdata)

  end subroutine read_ffi_glider


  ! Copied from m_read_ifremer_argo.
  !
  subroutine grid_readxyz(nx, ny, lat, lon, depth)
    integer, intent(in) :: nx, ny
    real(8), dimension(nx, ny), intent(inout) :: lat, lon, depth

    logical :: exists
    character(len = 128) :: fname
    
    fname = 'newpos.uf'
    inquire(file = fname, exist = exists)
    if (.not. exists) then
       print *, 'grid_readxyz(): ERROR: "', trim(fname), '" does not exist'
       stop
    end if
    open(10, file = fname, form = 'unformatted', status = 'old')
    print *, '  grid_readxyz(): reading "', trim(fname), '"...'
    read(10) lat, lon
    close(10)

    write(fname, '(a, i3.3, a, i3.3, a)') 'depths', nx, 'x', ny, '.uf'
    inquire(file = fname, exist = exists)
    if (.not. exists) then
       print*, 'grid_readxyz(): ERROR: "', trim(fname), '" does not exist'
       stop
    end if
    open (unit = 10, file = fname, status = 'old', form = 'unformatted')
    print *, '  grid_readxyz(): reading "', trim(fname), '"...'
    read(10) depth
    close(10)
  end subroutine grid_readxyz

end module m_read_FFI_glider
