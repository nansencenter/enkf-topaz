module m_read_metno_icec

contains

  subroutine read_metno_icec_repro(fname, data, gr)
    use nfw_mod
    use mod_measurement
    use mod_grid
    implicit none

    character(*), intent(in) :: fname
    type (measurement), allocatable, intent(out) :: data(:)
    type(grid), intent(out) :: gr

    logical :: ex
    integer :: ncid
    integer :: xc_id, yc_id
    integer :: nx, ny
    integer :: lon_id, lat_id, icec_id, flag_id,cflag_id
    real*4, allocatable :: lon(:,:), lat(:,:), icec(:,:)
    integer, allocatable :: flag(:,:),cflag(:,:)
    integer*1, allocatable :: tmp_flag(:,:),tmp_cflag(:,:)
    integer*2,allocatable :: tmp_icec(:,:)

    integer :: i, j, nobs

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
    allocate(icec(nx, ny))
    allocate(tmp_icec(nx, ny),tmp_flag(nx,ny),tmp_cflag(nx,ny))
    allocate(flag(nx, ny))
    allocate(cflag(nx, ny))
    call nfw_inq_varid(fname, ncid, 'lon', lon_id)
    call nfw_inq_varid(fname, ncid, 'lat', lat_id)
    call nfw_inq_varid(fname, ncid, 'ice_conc', icec_id)
!    call nfw_inq_varid(fname, ncid, 'confidence_level', cflag_id)
    call nfw_inq_varid(fname, ncid, 'status_flag', flag_id)
    call nfw_get_var_real(fname, ncid, lon_id, lon)
    call nfw_get_var_real(fname, ncid, lat_id, lat)
    call nfw_get_var_int2(fname, ncid, icec_id, tmp_icec)
!    call nfw_get_var_int(fname, ncid, cflag_id, cflag)
    call nfw_get_var_int1(fname, ncid, flag_id, tmp_flag)
    call nfw_close(fname, ncid)

    print *, 'filling the measurements array...'

    allocate(data(nx * ny))
    do j=1,ny
      do i=1,nx
        flag(i,j)=int(tmp_flag(i,j))
        icec(i,j)=real(tmp_icec(i,j))
      end do
   end do
    ! 0.995 is the max allowed by the model
    where (9950.0d0 <= icec .and. icec <= 10000.d0)
       icec = 9950.0d0
    end where   


    nobs = 0
    do j = 1, ny
       do i = 1, nx
          nobs = nobs + 1
!          if (flag(i, j) /= 0.or.cflag(i,j)<3) then
          if (flag(i, j) < 0.or.flag(i,j)>10) then
             data(nobs) % status = .false.
             cycle
          end if
          data(nobs) % id = 'ICEC'
          data(nobs) % d = icec(i, j) * 1d-4
!          data(nobs) % var = max(1d-8 * std(i, j) ** 2, 0.01d0 + (0.5d0 - abs(0.5d0 - data(nobs) % d)) ** 2)
          data(nobs) % var = 0.01+(0.5-abs(0.5-data(nobs)%d))**2 

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
    print *, '  ', nobs, 'primary ICEC observations'
    print *, '  ', minval(data % d), ' <= icec <= ', maxval(data % d)

    gr = default_grid
    gr % nx = nx
    gr % ny = ny
    gr%reg = .true.
    gr % order = 2
    gr%ux = '10 km'
    gr%uy = '10 km'
    gr%set = .true.

    deallocate(lat, lon, icec, flag, cflag)
    deallocate(tmp_icec,tmp_flag,tmp_cflag)
  end subroutine read_metno_icec_repro

  subroutine read_metno_icec_norepro(fname, data, gr)
    use nfw_mod
    use mod_measurement
    use mod_grid
    implicit none

    character(*), intent(in) :: fname
    type (measurement), allocatable, intent(out) :: data(:)
    type(grid), intent(out) :: gr

    logical :: ex
    integer :: ncid
    integer :: xc_id, yc_id
    integer :: nx, ny
    integer :: lon_id, lat_id, icec_id, cfl_id, flag_id
    real*8, allocatable :: lon(:,:), lat(:,:), icec(:,:), cfl(:, :)
    integer, allocatable :: flag(:,:)

    integer :: i, j, nobs

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
    allocate(icec(nx, ny))
    allocate(cfl(nx, ny))
    allocate(flag(nx, ny))
    call nfw_inq_varid(fname, ncid, 'lon', lon_id)
    call nfw_inq_varid(fname, ncid, 'lat', lat_id)
    call nfw_inq_varid(fname, ncid, 'ice_conc', icec_id)
    call nfw_inq_varid(fname, ncid, 'confidence_level', cfl_id)
    call nfw_inq_varid(fname, ncid, 'status_flag', flag_id)
    call nfw_get_var_double(fname, ncid, lon_id, lon)
    call nfw_get_var_double(fname, ncid, lat_id, lat)
    call nfw_get_var_double(fname, ncid, icec_id, icec)
    call nfw_get_var_double(fname, ncid, cfl_id, cfl)
    call nfw_get_var_int(fname, ncid, flag_id, flag)
    call nfw_close(fname, ncid)

    print *, 'filling the measurements array...'

    allocate(data(nx * ny))

    ! 0.995 is the max allowed by the model
    where (9950.0d0 <= icec .and. icec <= 10000.d0)
       icec = 9950.0d0
    end where   


    nobs = 0
    do j = 1, ny
       do i = 1, nx
          nobs = nobs + 1
          if (flag(i, j) /= 0 .OR. cfl(i,j) < 4) then
             data(nobs) % status = .false.
             cycle
          end if
          data(nobs) % id = 'ICEC'
          data(nobs) % d = icec(i, j) * 1d-4
          if (cfl(i,j)==4) then
             data(nobs) % var = 0.25d0
         else if (cfl(i,j) == 5 ) then 
             data(nobs) % var = 0.01d0
          end if
!          data(nobs) % var =  0.01d0 + (0.5d0 - abs(0.5d0 - data(nobs) % d)) ** 2
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
    print *, '  ', nobs, 'primary ICEC observations'
    print *, '  ', minval(data % d), ' <= icec <= ', maxval(data % d)

    gr = default_grid
    gr % nx = nx
    gr % ny = ny
    gr%reg = .true.
    gr % order = 2
    gr%ux = '10 km'
    gr%uy = '10 km'
    gr%set = .true.

    deallocate(lat, lon, icec, cfl, flag)
  end subroutine read_metno_icec_norepro


end module m_read_metno_icec
