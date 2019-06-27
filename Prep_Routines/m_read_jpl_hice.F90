module m_read_jpl_hice

contains

  subroutine read_jpl_hice(fname, obstype, variance, nx, ny, data)
    use mod_measurement
    use m_oldtonew
    use m_confmap
    use m_bilincoeff
    use m_pivotp
    implicit none

    character(*), intent(in) :: fname
    character(*), intent(in) :: obstype
    real, intent(in) :: variance
    integer, intent(in) :: nx, ny
    type(measurement), allocatable, intent(out) :: data(:)

    type(measurement), allocatable :: tmpdata(:)
    real(8), dimension(nx, ny) :: modlat, modlon
    real(8), dimension(nx, ny) :: depths

    integer :: npoints
    integer :: i
    integer :: nobs
    real    :: lat, lon, tmp1, tmp2, h
    real    :: latnew, lonnew
    integer :: ipiv, jpiv

    open(101, file = trim(fname))
    read(101, *) npoints
    print *, '  ',trim(fname), ': ', npoints, ' data points'

    allocate(tmpdata(npoints))
    call confmap_init(nx, ny)
    call grid_readxyz(nx, ny, modlat, modlon, depths)

    nobs = 0
    do i = 1, npoints
       read(101, *) lat, lon, tmp1, tmp2, h
       if (h > 0.0d0 .and. h < 9990.0d0) then
          call oldtonew(lat, lon, latnew, lonnew)
          call pivotp(lonnew, latnew, ipiv, jpiv)

          if (ipiv < 1 .or. jpiv < 1 .or. ipiv > nx - 1 .or. jpiv > ny - 1) then
             cycle
          end if
          if (depths(ipiv, jpiv) < 10) then
             cycle
          end if

          nobs = nobs + 1
          tmpdata(nobs) % d = h / 100.0d0
          tmpdata(nobs) % id = obstype
          tmpdata(nobs) % var = variance
          tmpdata(nobs) % lon = lon
          tmpdata(nobs) % lat = lat
          tmpdata(nobs) % ipiv = ipiv
          tmpdata(nobs) % jpiv = jpiv
          tmpdata(nobs) % ns = 0 ! for a point (not gridded) measurement
          tmpdata(nobs) % date = 0 ! assimilate synchronously

          call bilincoeff(real(modlon), real(modlat), nx, ny, lon, lat, ipiv,&
               jpiv, tmpdata(nobs) % a1, tmpdata(nobs) % a2, tmpdata(nobs) % a3,&
               tmpdata(nobs) % a4)

          tmpdata(nobs) % status = .true. ! (active)
          tmpdata(nobs) % i_orig_grid = -1 ! not used
          tmpdata(nobs) % j_orig_grid = -1 ! not used
       end if
    end do
    close(101)

    print *, '  ', nobs, ' valid observations'

    allocate(data(nobs))
    do i = 1, nobs
       data(i) = tmpdata(i)
    end do
    deallocate(tmpdata)

  end subroutine read_jpl_hice


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

end module m_read_jpl_hice
