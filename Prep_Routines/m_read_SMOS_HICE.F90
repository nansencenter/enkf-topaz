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

#if defined (CYSMOSV2)
    integer, allocatable :: tmp_int(:,:)
    real*8, dimension(1) :: scalefac
    integer,dimension(1) :: fillval
#endif
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

#if defined (CYSMOSV2)
    allocate(tmp_int(nx,ny))
#endif
    !call nfw_inq_varid(fname, ncid, 'longitude', lon_id)
    !call nfw_inq_varid(fname, ncid, 'latitude', lat_id)
    call nfw_inq_varid(fname, ncid, 'lon', lon_id)
    call nfw_inq_varid(fname, ncid, 'lat', lat_id)
    call nfw_inq_varid(fname, ncid, 'analysis_sea_ice_thickness', hice_id)
    call nfw_inq_varid(fname, ncid, 'analysis_sea_ice_thickness_unc', hvar_id)
    print *, lon_id,lat_id,hice_id,hvar_id,fice_id
#if defined (CYSMOSV2)
    call nfw_inq_varid(fname, ncid, 'sea_ice_concentration', fice_id)

    call nfw_get_var_double(fname, ncid, lon_id, lon)
    call nfw_get_var_double(fname, ncid, lat_id, lat)

    call nfw_get_var_int(fname, ncid, hice_id, tmp_int)
    call nfw_get_att_double(fname, ncid, hice_id, 'scale_factor', scalefac)
    call nfw_get_att_int(fname, ncid, hice_id, '_FillValue', fillval)
    hice=real(tmp_int)*scalefac(1)
    where(tmp_int==fillval(1))  hice=0
 

    call nfw_get_var_int(fname, ncid, hvar_id, tmp_int)
    call nfw_get_att_double(fname, ncid, hvar_id,  'scale_factor', scalefac)
    call nfw_get_att_int(fname, ncid, hvar_id, '_FillValue', fillval)
    hvar=real(tmp_int)*scalefac(1)
    where(tmp_int==fillval(1))  hvar=0

    call nfw_get_var_int(fname, ncid, fice_id, tmp_int)
    call nfw_get_att_double(fname, ncid, fice_id,  'scale_factor', scalefac)
    call nfw_get_att_int(fname, ncid, fice_id, '_FillValue', fillval)
    fice=real(tmp_int)*scalefac(1)
    where(tmp_int==fillval(1))  fice=0

    deallocate(tmp_int)
#else
    call nfw_inq_varid(fname, ncid, 'ice_concentration', fice_id)

    call nfw_get_var_double(fname, ncid, lon_id, lon)
    call nfw_get_var_double(fname, ncid, lat_id, lat)
    call nfw_get_var_real(fname, ncid, hice_id, hice)
    call nfw_get_var_real(fname, ncid, hvar_id, hvar)
    call nfw_get_var_real(fname, ncid, fice_id, fice)

#endif
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
#if defined (CYSMOS_Error)
!ifdef CYSMOS_Error
        !data(nobs) % var =  (hvar(i,j)+min(0.5,0.1+0.15*hice(i,j))) ** 2 !  Offset 0.1-0.5 
        ! tuning at 20th April 2017
!        data(nobs) % var =  (hvar(i,j)+min(0.25,0.1+0.075*hice(i,j))) ** 2 ! Offset 0.1-0.25 
!elseif CYSMOS_Error2
        ! tuning at 20th May 2019 relative to the new version of CS2SMOS
        if (hice(i,j)<3) then
          data(nobs) % var =  (hvar(i,j)+max(0.02,0.1*exp(-hice(i,j)*1.5))) ** 2 ! Offset 0.1-0.02 
        else
          data(nobs) % var =  (hvar(i,j)+min(0.2,0.02*exp(1.8*(hice(i,j)-3)))) ** 2 ! Offset 0.02-0.2 
        endif
#else
        data(nobs) % var =  (hvar(i,j) * 1.0) ** 2 ! Exaggerate, factor 2 
#endif
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


  subroutine read_mltp4_hice(fname,data, gr)
    use mod_measurement
    use mod_grid
    use netcdf
    use nfw_mod
    implicit none

    character(*), intent(in) :: fname
    type (measurement), allocatable, intent(out) :: data(:)
    type(grid), intent(out) :: gr
    real, parameter :: Tkmax=10.0,w0=2.0
    logical :: ex
    integer :: ncid
    integer :: xc_id, yc_id
    integer :: nx, ny
    integer :: lon_id, lat_id, hice_id, hvar_id
    real*4, allocatable :: lon(:,:), lat(:,:), hice(:,:)
    real*4, allocatable :: hvar(:,:),fice(:,:)
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
    allocate(fice(nx, ny))
    call nfw_inq_varid(fname, ncid, 'longitude', lon_id)
    call nfw_inq_varid(fname, ncid, 'latitude', lat_id)
    call nfw_inq_varid(fname, ncid, 'sit', hice_id)
    call nfw_inq_varid(fname, ncid, 'sic', hvar_id)
    call nfw_get_var_real(fname, ncid, lon_id, lon)
    call nfw_get_var_real(fname, ncid, lat_id, lat)
    call nfw_get_var_real(fname, ncid, hice_id, hice)
    call nfw_get_var_real(fname, ncid, hvar_id,fice)
    call nfw_close(fname, ncid)

    print *, 'filling the measurements array...'
    allocate(data(nx * ny))

    nobs=0
    do j=1, ny,2
      do i=1,nx,2
        nobs = nobs + 1
        if (fice(i, j) <= 0.1 .or.hice(i,j)>Tkmax.or.hice(i,j)<= 0.05) then
          data(nobs) % status = .false.
          cycle
        end if
        data(nobs) % id = 'HICE'
        data(nobs) % d = hice(i, j) 
        if (hice(i,j)<1.) then
          data(nobs) % var =  (max(0.15,0.3*exp(-hice(i,j)*2.0))) ** 2     !  

        else if (hice(i,j)<2.0) then
          data(nobs) % var =  (w0*max(0.15,0.3*exp(-hice(i,j)*2.0))) ** 2     !  

        else
          data(nobs) % var =  (w0*0.15*exp(.15*(hice(i,j)-2.0))) ** 2  !  
        endif
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

    print *, '  ', nobs, 'primary ML-based HICE from TP4'
    print *, '  ', minval(data % d), ' <= hice <= ', maxval(data % d)

     gr = default_grid
    gr % nx = nx
    gr % ny = ny
    gr%reg = .true.
    gr % order = 2
    gr%ux = '12.5 km'
    gr%uy = '12.5 km'
    gr%set = .true.

    deallocate(lat, lon, hice, hvar)
 

  end subroutine read_mltp4_hice
 

end module m_read_SMOS_HICE
