module m_read_OSISAF_data

contains

  ! Reads OSISAF ice drift data
  ! 2012-11-13: Geir Arne Waagbø (met.no)
  subroutine read_OSISAF_data(driftfile, gr, data, numdata, var, offset)
    use nfw_mod
    use mod_measurement
    use mod_grid
    implicit none

    character(*), intent(in) :: driftfile
    integer, intent(in) :: numdata
    type(measurement), dimension(numdata),intent(inout) :: data
    type(grid), intent(in) :: gr
    real, intent(in) :: var
    character(len=1), intent(in) :: offset

    integer :: dimids(2)
    integer , dimension(2) :: dimsizes

    integer :: lon_id, lat_id, dx_id, dy_id, qflag_id
    real*8, dimension(:,:), allocatable :: drlon, drlat, drdX, drdY
    integer, dimension(:,:), allocatable :: qflag

    integer :: ncid
    real*8, dimension(1) :: fillval

    integer :: i,j,k,icomp
    integer :: drnx, drny
    logical :: valid

    ! Get dimensions of drift file
    call nfw_open(driftfile, nf_nowrite, ncid)
    call nfw_inq_varid(driftfile, ncid, 'dX', dx_id)
    call nfw_inq_vardimid(driftfile, ncid, dx_id, dimids)
    do i = 1, 2
       call nfw_inq_dimlen(driftfile, ncid, dimids(i), dimsizes(i))
    end do
    drnx=dimsizes(1)
    drny=dimsizes(2)

    if (gr % reg) then
       print *,'NB: OSISAF data should be specified as irregular !'
       print *,'    Currently it is set as regular..'
       print *,'(read_OSISAF_data)'
       call exit(1) 
    end if

    ! Which should match numdata dimension 
    ! NB !!! Mult by 2 for two vector components
    if (2*drnx*drny /= numdata .or. 2*drnx*drny /= gr%nx) then
       print *,'Different dimensions - data file and specified'
       print *,'dimsizes(1)=',drnx
       print *,'dimsizes(2)=',drny
       print *,'nx         =',gr%nx
       print *,'(read_OSISAF_data)'
       call exit(1) 
    end if

    ! Read data from drift file
    allocate(drlon(drnx,drny))
    allocate(drlat(drnx,drny))
    allocate(drdX(drnx,drny))
    allocate(drdY(drnx,drny))
    allocate(qflag(drnx,drny))

    call nfw_inq_varid(driftfile, ncid, 'lon', lon_id)
    call nfw_get_var_double(driftfile, ncid, lon_id, drlon)
    call nfw_inq_varid(driftfile, ncid, 'lat', lat_id)
    call nfw_get_var_double(driftfile, ncid, lat_id, drlat)
    call nfw_inq_varid(driftfile, ncid, 'dX', dx_id)
    call nfw_get_var_double(driftfile, ncid, dx_id, drdX)

    ! Change dY_v1p4 to dY from version 1.4 of OSISAF product file
    call nfw_inq_varid(driftfile, ncid, 'dY_v1p4', dy_id)
    call nfw_get_var_double(driftfile, ncid, dy_id, drdY)

    call nfw_get_att_double(driftfile, ncid, dx_id, '_FillValue', fillval)

    where (abs(drdX - fillval(1)) < 1e-4 * fillval(1))
       drdX = gr % undef
    end where

    call nfw_get_att_double(driftfile, ncid, dy_id, '_FillValue', fillval)

    where (abs(drdY - fillval(1)) < 1e-4 * fillval(1))
       drdY = gr % undef
    end where

    ! Change status_flag_v1p4 to status_flag from version 1.4 of OSISAF product file
    call nfw_inq_varid(driftfile, ncid, 'status_flag_v1p4', qflag_id)
    call nfw_get_var_int(driftfile, ncid, qflag_id, qflag)

    call nfw_close(driftfile, ncid)

    k = 0
    do icomp = 1, 2
       do j = 1, drny
          do i = 1, drnx
             k = k + 1

             valid = qflag(i,j) == 30 ! Only use observations with quality_flag==30
             if (icomp==1) then
                data(k)%id = 'DX'//offset
                data(k)%d = drdX(i,j)
                valid = valid .and.  abs((drdX(i,j)-gr%undef) / gr%undef) > 1e-4
             else
                data(k)%id = 'DY'//offset
                data(k)%d = drdY(i,j)
                valid =  valid .and. abs((drdY(i,j)-gr%undef) / gr%undef) > 1e-4
             end if

             if (.not. valid .or. mod(i,2)==0 .or. mod(j,2)==0) then
                ! Skip invalid observations, or observations on grid points with
                ! even i- or j-indices (to avoid over assimilation)
                data(k)%d = gr%undef
             end if

             data(k)%ipiv = i  ! Not really used for ice drift
             data(k)%jpiv = j  ! Not really used for ice drift
             data(k)%i_orig_grid = i ! Used for ice drift
             data(k)%j_orig_grid = j ! Used for ice drift
             data(k)%lat=drlat(i,j)
             data(k)%lon=ang180(real(drlon(i,j)))
             ! Each vector represents the average drift of a 120kmx120km area of sea ice
             ! The a1 value should be in meters, although other values are in km
             data(k)%a1 = 1.4*60000 ! 1.4 represents square root of 2
             data(k)%ns = 1
             data(k)%var = var ! fom idrft.hdr specification
             data(k)%depth = 0.0
             data(k)%status = valid
          enddo
       enddo
    enddo
    print *, 'Number of data read:', k, gridpoints(gr)
  end subroutine


  ! Reads OSISAF ice drift data
  ! 2012-11-13: Geir Arne Waagbø (met.no)
  ! 2014-9-18:  in order to consistently read the dY from 2010 to 2013
  ! considering before 2011/11/13 no dY_v1p4 
  ! Adding the uncertainty of DX/DY from provider in June 2023 
  subroutine read_OSISAF_data_lv4(driftfile, gr, data, numdata, var, offset)
    use nfw_mod
    use mod_measurement
    use mod_grid
    implicit none

    character(*), intent(in) :: driftfile
    integer, intent(in) :: numdata
    type(measurement), dimension(numdata),intent(inout) :: data
    type(grid), intent(in) :: gr
    real, intent(in) :: var
    character(len=1), intent(in) :: offset

    integer :: dimids(2)
    integer , dimension(2) :: dimsizes

    integer :: lon_id, lat_id, dx_id, dy_id, qflag_id
    real*8, dimension(:,:), allocatable :: drlon, drlat, drdX, drdY
    integer, dimension(:,:), allocatable :: qflag

    integer :: uncertdxy_id
    real*8, dimension(:,:), allocatable :: undxy 

    integer :: ncid
    real*8, dimension(1) :: fillval

    integer :: i,j,k,icomp
    integer :: drnx, drny
    logical :: valid

    ! Get dimensions of drift file
    call nfw_open(driftfile, nf_nowrite, ncid)
    call nfw_inq_varid(driftfile, ncid, 'dX', dx_id)
    call nfw_inq_vardimid(driftfile, ncid, dx_id, dimids)
    do i = 1, 2
       call nfw_inq_dimlen(driftfile, ncid, dimids(i), dimsizes(i))
    end do
    drnx=dimsizes(1)
    drny=dimsizes(2)

    if (gr % reg) then
       print *,'NB: OSISAF data should be specified as irregular !'
       print *,'    Currently it is set as regular..'
       print *,'(read_OSISAF_data)'
       call exit(1) 
    end if

    ! Which should match numdata dimension 
    ! NB !!! Mult by 2 for two vector components
    if (2*drnx*drny /= numdata .or. 2*drnx*drny /= gr%nx) then
       print *,'Different dimensions - data file and specified'
       print *,'dimsizes(1)=',drnx
       print *,'dimsizes(2)=',drny
       print *,'nx         =',gr%nx
       print *,'(read_OSISAF_data)'
       call exit(1) 
    end if

    ! Read data from drift file
    allocate(drlon(drnx,drny))
    allocate(drlat(drnx,drny))
    allocate(drdX(drnx,drny))
    allocate(drdY(drnx,drny))
    allocate(qflag(drnx,drny))
    allocate(undxy(drnx,drny))

    call nfw_inq_varid(driftfile, ncid, 'lon', lon_id)
    call nfw_get_var_double(driftfile, ncid, lon_id, drlon)
    call nfw_inq_varid(driftfile, ncid, 'lat', lat_id)
    call nfw_get_var_double(driftfile, ncid, lat_id, drlat)
    call nfw_inq_varid(driftfile, ncid, 'dX', dx_id)
    call nfw_get_var_double(driftfile, ncid, dx_id, drdX)

    ! Change dY_v1p4 into dY OSISAF product file
    call nfw_inq_varid(driftfile, ncid, 'dY', dy_id)
    call nfw_get_var_double(driftfile, ncid, dy_id, drdY)

    ! Adding the uncertainty of DX/DY (one standard deviation)
    call nfw_inq_varid(driftfile, ncid, 'uncert_dX_and_dY', uncertdxy_id)
    call nfw_get_var_double(driftfile, ncid, uncertdxy_id, undxy)

    call nfw_get_att_double(driftfile, ncid, dx_id, '_FillValue', fillval)

    where (abs(drdX - fillval(1)) < 1e-4 * fillval(1))
       drdX = gr % undef
    end where

    call nfw_get_att_double(driftfile, ncid, dy_id, '_FillValue', fillval)

    where (abs(drdY - fillval(1)) < 1e-4 * fillval(1))
       drdY = gr % undef
    end where

    where (abs(undxy - fillval(1)) < 1e-4 * fillval(1))
       drdY = gr % undef
       drdX = gr % undef
    end where

    !print *,'fillval:', fillval(1),gr%undef
    !print *,        abs((fillval(1)-gr%undef) / gr%undef)
    ! Change status_flag_v1p4 to status_flag from version 1.4 of OSISAF product file
    call nfw_inq_varid(driftfile, ncid, 'status_flag', qflag_id)
    call nfw_get_var_int(driftfile, ncid, qflag_id, qflag)

    call nfw_close(driftfile, ncid)

    k = 0
    do icomp = 1, 2
       do j = 1, drny
          do i = 1, drnx
             k = k + 1

             valid = qflag(i,j) == 30 ! Only use observations with quality_flag==30
             if (icomp==1) then
                data(k)%id = 'DX'//offset
                data(k)%d = drdX(i,j)
                valid = valid .and.  abs((drdX(i,j)-gr%undef) / gr%undef) > 1e-4
                valid = valid .and.  abs(drdX(i,j))<gr%undef
             else
                data(k)%id = 'DY'//offset
                data(k)%d = drdY(i,j)
                valid =  valid .and. abs((drdY(i,j)-gr%undef) / gr%undef) > 1e-4
                valid = valid .and.  abs(drdY(i,j))<gr%undef
             end if

             if (.not. valid .or. mod(i,2)==0 .or. mod(j,2)==0) then
                ! Skip invalid observations, or observations on grid points with
                ! even i- or j-indices (to avoid over assimilation)
                data(k)%d = gr%undef
             end if

             data(k)%ipiv = i  ! Not really used for ice drift
             data(k)%jpiv = j  ! Not really used for ice drift
             data(k)%i_orig_grid = i ! Used for ice drift
             data(k)%j_orig_grid = j ! Used for ice drift
             data(k)%lat=drlat(i,j)
             data(k)%lon=ang180(real(drlon(i,j)))
             ! Each vector represents the average drift of a 120kmx120km area of sea ice
             ! The a1 value should be in meters, although other values are in km
             data(k)%a1 = 1.4*60000 ! 1.4 represents square root of 2
             data(k)%ns = 1
             !data(k)%var = var ! fom idrft.hdr specification
             data(k)%var = max(var,undxy(i,j)**2) ! fom idrft.hdr specification
             data(k)%depth = 0.0
             data(k)%status = valid
          enddo
       enddo
    enddo
    print *, 'Number of data read:', k, gridpoints(gr)
  end subroutine

end module


