










module m_read_SMOS_HICE


contains


  subroutine read_smos_hice(fname,data, gr)
    use mod_measurement
    use mod_grid
    use netcdf
    use nfw_mod
    implicit none

    character(*), intent(in) :: fname
    type (measurement), allocatable, intent(out) :: data(:)
    type(grid), intent(out) :: gr
    real, parameter :: Tkmax=0.4
    logical :: ex
    integer :: ncid
    integer :: xc_id, yc_id
    integer :: nx, ny
    integer :: lon_id, lat_id, hice_id, hvar_id
    real*4, allocatable :: lon(:,:), lat(:,:), hice(:,:)
    real*4, allocatable :: hvar(:,:)
!    integer, allocatable :: flag(:,:),cflag(:,:)

    integer :: i, j, nobs

    print *, 'reading "', trim(fname), '"...'

    inquire(file = trim(fname), exist = ex)
    if (.not. ex) then
       print *, 'ERROR: file "', trim(fname), '" not found'
       stop
    end if


    call nfw_open(fname, nf_nowrite, ncid)
    call nfw_inq_dimid(fname, ncid, 'x', xc_id)
    call nfw_inq_dimid(fname, ncid, 'y', yc_id)
    call nfw_inq_dimlen(fname, ncid, xc_id, nx)
    call nfw_inq_dimlen(fname, ncid, yc_id, ny)
    print *, '  nx = ', nx
    print *, '  ny = ', ny
    allocate(lon(nx, ny))
    allocate(lat(nx, ny))
    allocate(hice(nx, ny))
    allocate(hvar(nx, ny))
    call nfw_inq_varid(fname, ncid, 'longitude', lon_id)
    call nfw_inq_varid(fname, ncid, 'latitude', lat_id)
    call nfw_inq_varid(fname, ncid, 'sea_ice_thickness', hice_id)
    call nfw_inq_varid(fname, ncid, 'ice_thickness_uncertainty', hvar_id)
    call nfw_get_var_real(fname, ncid, lon_id, lon)
    call nfw_get_var_real(fname, ncid, lat_id, lat)
    call nfw_get_var_real(fname, ncid, hice_id, hice)
    call nfw_get_var_real(fname, ncid, hvar_id, hvar)
    call nfw_close(fname, ncid)

    print *, 'filling the measurements array...'
    allocate(data(nx * ny))

    nobs=0
    do j=1, ny
      do i=1,nx
        nobs = nobs + 1
        if (hvar(i, j) <= 0 .or.hvar(i,j)>=5.or.hice(i,j)>Tkmax.or.hice(i,j)<= 0) then
          data(nobs) % status = .false.
          cycle
        end if
        data(nobs) % id = 'HICE'
        data(nobs) % d = hice(i, j) 
        data(nobs) % var =  (hvar(i,j) * 1.0) ** 2 ! Exaggerate, factor 2 
 
        data(nobs) % ipiv = i
        data(nobs) % jpiv = j
        data(nobs) % lon = lon(i, j)
        data(nobs) % lat = lat(i, j)
        data(nobs) % a1 = 1e10
        data(nobs) % a2 = 1e10
        data(nobs) % a3 = 1e10
        data(nobs) % a4 = 1e10
        data(nobs) % ns = 1
        data(nobs) % date = 0
        data(nobs) % depth = 0.0
        data(nobs) % status = .true.
 
      end do
    end do

    print *, '  ', nobs, 'primary SMOS HICE observations'
    print *, '  ', minval(data % d), ' <= hice <= ', maxval(data % d)

     gr = default_grid
    gr % nx = nx
    gr % ny = ny
    gr%reg = .true.
    gr % order = 2
    gr%ux = '10 km'
    gr%uy = '10 km'
    gr%set = .true.

    deallocate(lat, lon, hice, hvar)
 

  end subroutine read_smos_hice
 





  subroutine read_cysmos_hice(fname,data, gr)
    use mod_measurement
    use mod_grid
    use netcdf
    use nfw_mod
    implicit none

    character(*), intent(in) :: fname
    type (measurement), allocatable, intent(out) :: data(:)
    type(grid), intent(out) :: gr
    real, parameter :: Tkmax=10
    logical :: ex
    integer :: ncid
    integer :: xc_id, yc_id
    integer :: nx, ny
    integer :: lon_id, lat_id, fice_id,hice_id, hvar_id
    real*8, allocatable :: lon(:,:), lat(:,:)
    real*4, allocatable :: fice(:,:),hice(:,:),hvar(:,:)
    integer, allocatable :: mask(:,:)

    integer :: i, j, nobs
    integer :: i1,j1

    print *, 'reading "', trim(fname), '"...'

    inquire(file = trim(fname), exist = ex)
    if (.not. ex) then
       print *, 'ERROR: file "', trim(fname), '" not found'
       stop
    end if


    call nfw_open(fname, nf_nowrite, ncid)
    call nfw_inq_dimid(fname, ncid, 'xc', xc_id)
    call nfw_inq_dimid(fname, ncid, 'yc', yc_id)
    call nfw_inq_dimlen(fname, ncid, xc_id, nx)
    call nfw_inq_dimlen(fname, ncid, yc_id, ny)
    print *, '  nx = ', nx
    print *, '  ny = ', ny
    allocate(lon(nx, ny))
    allocate(lat(nx, ny))
    allocate(hice(nx, ny))
    allocate(hvar(nx, ny))
    allocate(fice(nx, ny),mask(nx,ny))
    call nfw_inq_varid(fname, ncid, 'longitude', lon_id)
    call nfw_inq_varid(fname, ncid, 'latitude', lat_id)
    call nfw_inq_varid(fname, ncid, 'analysis_thickness', hice_id)
    call nfw_inq_varid(fname, ncid, 'analysis_thickness_unc', hvar_id)
    call nfw_inq_varid(fname, ncid, 'ice_concentration', fice_id)
    print *, lon_id,lat_id,hice_id,hvar_id,fice_id

    call nfw_get_var_double(fname, ncid, lon_id, lon)
    call nfw_get_var_double(fname, ncid, lat_id, lat)
    call nfw_get_var_real(fname, ncid, hice_id, hice)
    call nfw_get_var_real(fname, ncid, hvar_id, hvar)
    call nfw_get_var_real(fname, ncid, fice_id, fice)
    call nfw_close(fname, ncid)

    print *, 'filling the measurements array...'
    allocate(data(nx * ny))
    mask=0;
    where(hvar>0) mask=1
    where(hvar<0) mask=1

    nobs=0
    do j=1, ny
      do i=1,nx
        nobs = nobs + 1
        i1=min(nx,i+1);
        j1=min(ny,j+1);
     !   if (hvar(i, j) <= 0 .or.hvar(i,j)>1e10.or.hice(i,j)>Tkmax.or.hice(i,j)<= 0) then
    !      data(nobs) % status = .false.
    !      cycle
    !    end if
        if (mask(i,j)==0) then
          data(nobs) % status = .false.
          cycle
        end if 

    !    print '(5f40.1)', hvar(i, j), hice(i,j),lon(i,j),lat(i,j),fice(i,j)

        data(nobs) % id = 'HICE'
        data(nobs) % d = real(hice(i, j)) 
        data(nobs) % var =  (hvar(i,j) * 1.0) ** 2 ! Exaggerate, factor 2 
        data(nobs) % ipiv = 0 
        data(nobs) % jpiv = 0 
        data(nobs) % lon = lon(i, j)
        data(nobs) % lat = lat(i, j)
        data(nobs) % a1 = 1e10
        data(nobs) % a2 = 1e10
        data(nobs) % a3 = 1e10
        data(nobs) % a4 = 1e10
        data(nobs) % ns = 1
        data(nobs) % date = 0
        data(nobs) % depth = 0.0
        data(nobs) % status = .true.
 
      end do
    end do

    print *, '  ', nobs, 'primary SMOS HICE observations'
    print *, '  ', minval(data % d), ' <= hice <= ', maxval(data % d)

     gr = default_grid
    gr % nx = nx
    gr % ny = ny
    gr%reg = .true.
    gr % order = 2
    gr%ux = '12 km'
    gr%uy = '12 km'
    gr%set = .true.

    deallocate(lat, lon, hice, hvar,fice)
 

  end subroutine read_cysmos_hice



end module m_read_SMOS_HICE
