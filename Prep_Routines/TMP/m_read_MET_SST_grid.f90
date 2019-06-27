










module m_read_MET_SST_grid
  ! Reads the CLS SST NetCDF dimensions
contains
  subroutine read_MET_SST_grid(filename,gr)
    !use mod_dimensions
    use mod_grid
    use nfw_mod
    implicit none

    character(len=80), intent(in) :: filename
    type(grid),        intent(out) :: gr
    logical :: ex
    !dimension ids
    integer :: lon_ID,lat_ID

    ! Array dimensions
    integer :: nblon,nblat

    integer :: ncid
    real*8, allocatable :: lat(:), lon(:)

    print *, 'read_MET_SST_grid():'

    gr = default_grid
    inquire(file=trim(filename),exist=ex)
    if(ex) then
         call nfw_open(filename, nf_nowrite, ncid)
         print *, '  found "', trim(filename), '"...'
          ! Get dimension id in netcdf file ...
          call nfw_inq_dimid(filename, ncid, 'lon', lon_ID)
          call nfw_inq_dimid(filename, ncid, 'lat', lat_ID)
          ! Get dimension length from id
          call nfw_inq_dimlen(filename, ncid, lon_ID, nblon)
          call nfw_inq_dimlen(filename, ncid, lat_ID, nblat)
          print*, 'Dimensions lon,lat:', nblon, nblat
          allocate(lon(nblon), lat(nblat))
          call nfw_inq_varid(filename, ncid, 'lon', lon_ID)
          call nfw_inq_varid(filename, ncid, 'lat', lat_ID)
          call nfw_get_var_double(filename, ncid, lon_ID, lon)
          call nfw_get_var_double(filename, ncid, lat_ID, lat)
          call nfw_close(filename, ncid)


          gr%nx=nblon
          gr%ny=nblat
          gr%x0=lon(1)
          gr%y0=lat(1)
          gr%dx=lon(2)-lon(1)
          gr%dy=lat(2)-lat(1)

          gr%reg = .true.
          gr%order = 2
          gr%ux = 'm'
          gr%uy = 'm'
          gr%set = .true.
          deallocate(lon,lat)
    endif
  end subroutine read_MET_SST_grid
  
end module m_read_MET_SST_grid
