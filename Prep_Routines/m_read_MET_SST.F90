module m_read_MET_SST
  ! Reads CLS SLA data after having read the grid in read_CLS_SST_grid

  integer, parameter, private :: STRLEN = 512

contains

  subroutine read_MET_SST(filename,gr,data,Dflg)
    use mod_measurement
    use mod_grid
    use m_spherdist
    use netcdf
    use nfw_mod
    implicit none

    type (measurement),  intent(inout) :: data(:)
    type (grid),         intent(inout) :: gr       ! measurement grid
    character(len=80),      intent(in) :: filename
    integer,                intent(in) :: Dflg     ! different product formats

    ! Variable ids
    integer :: lon_ID, lat_ID,vsst_ID, vstd_ID, vmask_ID
    ! Data arrays
    real*8, allocatable :: sst(:,:), lon(:), lat(:), std(:,:)
    integer, allocatable :: mask(:,:)
    integer :: ncid ! observations
    real*8, dimension(1) :: undef_sst
    integer :: i, j, count1
    real, parameter :: eps = 0.01  ! test for undefined values

    ! enlarge times for the error variance changed in Sep 2025
    real            :: Lscale

    ! filen name
    logical         :: ex

    print *, 'read_MET_SST:'

    inquire(file=trim(filename),exist=ex)
    if (ex) then
       ! Reading the observation file 
       call nfw_open(filename, nf_nowrite, ncid)
       ! Get dimension id in netcdf file ...
       !nb total of data
       allocate(lon(gr%nx), lat(gr%ny), sst(gr%nx,gr%ny), std(gr%nx, gr%ny), mask(gr%nx, gr%ny))

       ! Variable ids in netcdf file
       call nfw_inq_varid(filename, ncid, 'lat', lat_ID)
       call nfw_inq_varid(filename, ncid,'lon', lon_ID)
       call nfw_inq_varid(filename, ncid,'analysed_sst' ,vsst_ID)
       if (Dflg==1) then   ! to avoid the discontinuety in the real operational run
          ! nrt product OSTIA SST
          call nfw_inq_varid(filename, ncid,'analysis_error' ,vstd_ID)
          ! enlarge times default used in TOPAZ4 
          Lscale=6.0
      
       elseif (Dflg==2) then
          ! ESACCI SST <=2016
          call nfw_inq_varid(filename, ncid,'analysed_sst_uncertainty' ,vstd_ID)
          Lscale=4.0
       else
          ! C3S2 SST <=2024
          call nfw_inq_varid(filename, ncid,'analysis_error_sst' ,vstd_ID)
          Lscale=6.0
       endif
       call nfw_inq_varid(filename, ncid,'mask' ,vmask_ID)
       
       ! Variable _FillValue attributes
       call nfw_get_att_double(filename, ncid, vsst_ID, '_FillValue', undef_sst(1))
       gr % undef = undef_sst(1)
       
       ! actual variable values (for dimensions of var -- see ncdump, or improve this program)
       ! NB: note that index dimensions are different between fortran and C internals. 
       ! "ncdump" gives C internals.
       call nfw_get_var_double(filename, ncid, lon_ID, lon)
       call nfw_get_var_double(filename, ncid, lat_ID, lat)
       call nfw_get_var_double(filename, ncid, vsst_ID, sst)
       call nfw_get_var_double(filename, ncid, vstd_ID, std)
       call nfw_get_var_int(filename, ncid, vmask_ID, mask)
       print '(1x, a, 2f10.2)', '    range Lon = ', minval(lon), maxval(lon)
       print '(1x, a, 2f10.2)', '    range Lat = ', minval(lat), maxval(lat)
       print '(1x, a, 2f10.2)', '    range sst (K) = ', minval(sst), maxval(sst)
       print '(1x, a, 2i10)',   '    range mask = ', minval(mask), maxval(mask)
       call nfw_close(filename, ncid)
       count1=1
       do i=1,gr%nx
          do j=1,gr%ny
            !here we only consider:
            !data above -30 of lat; valid, within reasonable range,
            ! and only open ocean (mask == 1)
            !                 
            if (lat(j) > -30 .and.&
                     abs(sst(i,j)-undef_sst(1)) > eps .and.&
                     mask(i,j) == 1 .and. & 
                     sst(i,j) > -190 .and. &
                     sst(i,j) < 4500 .and. &
                     std(i,j) > 0.0) then
                   data(count1)%id = 'SST'
                   data(count1)%d = sst(i,j)*0.01  
                   data(count1)%ipiv = count1 !whatever it is filled afterwards
                   data(count1)%jpiv = 1   !whaterver it is filled afterwards
                   data(count1)%lat=lat(j)
                   data(count1)%lon=lon(i)
                   data(count1)%a1 = spherdist(real(lon(i))-.5*gr%dx,real(lat(j))-.5*gr%dy,real(lon(i))+.5*gr%dx,real(lat(j))+.5*gr%dy)
                   data(count1)%ns = 1 ! 1 for data with a spatial extent
                   data(count1)%var = (std(i,j) * 0.01 * Lscale) ** 2 ! Exaggerate, factor 8
                   data(count1)%date = 0
                   data(count1)%depth = 0.0
                   data(count1)%status = .true.
                   count1=count1+1
            endif
          enddo   !i
        enddo    !j
       print*, '    # of obs read = ', count1
       deallocate(lat, lon, sst, mask)
    end if ! ex
    print *, 'MAX var(SST) = ', maxval(data(1 : count1) % var)
  end subroutine read_MET_SST
end module m_read_MET_SST
