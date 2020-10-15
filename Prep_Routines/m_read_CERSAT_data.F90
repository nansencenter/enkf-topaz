module m_read_CERSAT_data

contains



  subroutine read_CERSAT_data_rep(driftfile, gr, data, numdata, var)
    use nfw_mod
    use mod_measurement
    use mod_grid
    use m_spherdist
    implicit none

    character(*), intent(in) :: driftfile
    integer, intent(in) :: numdata
    type(measurement), dimension(numdata) :: data
    type(grid), intent(in) :: gr
    real, intent(in) :: var

    integer :: dimids(4)
    integer , dimension(4) :: dimsizes

    integer :: lon_id, lat_id, zon_id, mer_id, qua_id
    real*8, dimension(:,:), allocatable :: drlon, drlat, drmer, drzon
    real*8, dimension(:),   allocatable :: tmplon, tmplat
    real*8, dimension(:,:,:), allocatable :: tmpmer, tmpzon

    integer :: ncid, varid
    real*8, dimension(1) :: scalefac, fillval, addoffset

    integer :: i,j,k,icomp,idate
    integer :: drnx, drny
    logical :: valid
    real    :: Sunit

    ! Get dimensions of drift file
    call nfw_open(driftfile, nf_nowrite, ncid)
    call nfw_inq_varid(driftfile, ncid,'eastward_sea_ice_velocity', &
         varid)
    call nfw_inq_vardimid(driftfile, ncid, varid, dimids)
    do i = 1, 4
       call nfw_inq_dimlen(driftfile, ncid, dimids(i), dimsizes(i))
    end do

    print *,'numdata=',numdata
    print *,'dimsizes(1)=',dimsizes(1)
    print *,'dimsizes(2)=',dimsizes(2)
    print *,'dimsizes(3)=',dimsizes(3)
    print *,'dimsizes(4)=',dimsizes(4)

    ! Which should match numdata dimension 
    ! NB !!! Mult by 2 for two vector components
    if (4 * dimsizes(1) * dimsizes(2) /= numdata) then 
       print *,'Different dimensions - data file and specified'
       print *,'nx         =',gr%nx
       print *,'(read_CERSAT_data_rep)'
       call exit(1) 
    end if

    ! Read data from drift file
    drnx=dimsizes(1)
    drny=dimsizes(2)
    print *,'dimsizes(1)=',dimsizes(1)
    print *,'dimsizes(2)=',dimsizes(2)
    allocate(drlon(drnx,drny))
    allocate(drlat(drnx,drny))
    allocate(tmplon(drnx), tmplat(drny))
    allocate(tmpmer(drnx,drny,2),tmpzon(drnx,drny,2))
    allocate(drmer(drnx,drny))
    allocate(drzon(drnx,drny))
    call nfw_inq_varid(driftfile, ncid, 'longitude', lon_id)
    call cersat_readfield(driftfile, ncid, lon_id, tmplon, drnx,0)
    call nfw_inq_varid(driftfile, ncid, 'latitude', lat_id)
    call cersat_readfield(driftfile, ncid, lat_id, tmplat, drny,0)

    do i=1, dimsizes(1)
       drlat(i,:)=tmplat
    end do
    do i=1, dimsizes(2)
       drlon(:,i)=tmplon
    end do

    call nfw_inq_varid(driftfile, ncid, 'eastward_sea_ice_velocity',zon_id)
    call cersat_readfield(driftfile, ncid, zon_id, tmpzon, 2*drnx * drny,1)
    call nfw_inq_varid(driftfile, ncid, 'northward_sea_ice_velocity', mer_id)
    call cersat_readfield(driftfile, ncid, mer_id, tmpmer, 2*drnx * drny,1)

    call nfw_get_att_double(driftfile,ncid, zon_id,'_FillValue', fillval)
    call nfw_get_att_double(driftfile,ncid, zon_id,'scale_factor',scalefac)
    call nfw_get_att_double(driftfile,ncid, zon_id, 'add_offset', addoffset)

    !where (abs(tmpzon - fillval(1)) <abs(1e-4 * fillval(1)) )
    where (abs(tmpzon - (fillval(1) * scalefac(1) + addoffset(1))) < &
         1e-4 * abs(fillval(1) * scalefac(1) + addoffset(1)))
       tmpzon = gr % undef
    end where

    call nfw_get_att_double(driftfile, ncid, mer_id, '_FillValue', fillval)
    call nfw_get_att_double(driftfile, ncid, mer_id, 'scale_factor', scalefac)
    call nfw_get_att_double(driftfile, ncid, mer_id, 'add_offset', addoffset)

    ! Flag zonal motion for fill values
    !where (abs(tmpmer - fillval(1)) <abs(1e-4 * fillval(1)) )
    where (abs(tmpmer - (fillval(1) * scalefac(1) + addoffset(1))) < &
         1e-4 * abs(fillval(1) * scalefac(1) + addoffset(1)))
       tmpmer = gr % undef
    end where
    call nfw_close(driftfile, ncid)


    Sunit=24*3600*0.001
    k = 0
    do icomp = 1, 2
       do idate=1,2       ! there are two times observations
          drmer=tmpmer(:,:,idate)
          drzon=tmpzon(:,:,idate)
      
       do j = 1, drny ! gr%ny
          do i = 1, drnx ! gr%nx
             k = k + 1
             valid=.true.
             if (icomp==1) then
                data(k)%id = 'VICE'
                data(k)%d = drmer(i,j)*Sunit ! Convert to km/day
                valid = valid .and.  abs( (drmer(i,j)-gr%undef)   / gr%undef)   > 1e-4
             else
                data(k)%id = 'UICE'
                data(k)%d = drzon(i,j)*Sunit ! Convert to km/day
                valid =  valid .and. abs( (drzon(i,j)-gr%undef)   / gr%undef)   > 1e-4
             end if

             if (.not.valid) then
                data(k)%d = gr%undef
             end if

             data(k)%ipiv = i  ! Not really used for ice drift
             data(k)%jpiv = j  ! Not really used for ice drift
             data(k)%i_orig_grid = i ! Used for ice drift
             data(k)%j_orig_grid = j ! Used for ice drift
             data(k)%lat=drlat(i,j)
             data(k)%lon=ang180(real(drlon(i,j)))
             !LB: Data support is assumed = a square grid cell
             !support diameter in meters stored in %a1 (tricky, isn't it ?)
             ! KAL -- hardcoded from data
             data(k)%a1 = 1.4 * 60000.0
             data(k)%ns = 1
             ! To be decided - obs units are meters O(1e4)
             ! CERSAT grid cells are of ~62.5 km - We assume the errors are
             ! roughly ~10 km
             !KAL data(k)%var = (10)**2
             data(k)%var = var ! fom idrft.hdr specification
             data(k)%depth = 0.0
             data(k)%date = 1 - idate 
             data(k)%status = valid
          enddo
       enddo
    enddo

    enddo
    deallocate(drlon,drlat,tmplon,tmplat,drmer,drzon)
    deallocate(tmpmer,tmpzon)
    print*, 'Number of data read:', k, gridpoints(gr)

  end subroutine read_CERSAT_data_rep








  subroutine read_CERSAT_data(driftfile, gr, data, numdata, var)
    use nfw_mod
    use mod_measurement
    use mod_grid
    use m_spherdist
    implicit none

    character(*), intent(in) :: driftfile
    integer, intent(in) :: numdata
    type(measurement), dimension(numdata) :: data
    type(grid), intent(in) :: gr
    real, intent(in) :: var

    integer :: dimids(2)
    integer , dimension(2) :: dimsizes

    integer :: lon_id, lat_id, zon_id, mer_id, qua_id
    real*8, dimension(:,:), allocatable :: drlon, drlat, drmer, drzon
    integer, dimension(:,:), allocatable :: qflag

    integer :: ncid, varid
    real*8, dimension(1) :: scalefac, fillval, addoffset

    integer :: i,j,k,icomp
    integer :: drnx, drny
    logical :: valid
    integer :: tmpint, bit(0:8)
#if defined (CERSAT_QF)
    real    :: tmpf
    integer :: ik
#endif

    ! Get dimensions of drift file
    call nfw_open(driftfile, nf_nowrite, ncid)
    call nfw_inq_varid(driftfile, ncid, 'zonal_motion', varid)
    call nfw_inq_vardimid(driftfile, ncid, varid, dimids)
    do i = 1, 2
       call nfw_inq_dimlen(driftfile, ncid, dimids(i), dimsizes(i))
    end do

    if (gr % reg) then
       print *,'NB: CERSAT data should be specified as irregular !'
       print *,'    Currently it is set as regular..'
       print *,'(read_CERSAT_data)'
       call exit(1) 
    end if

    ! Which should match numdata dimension 
    ! NB !!! Mult by 2 for two vector components
    if (2 * dimsizes(1) * dimsizes(2) /= numdata .or. &
         gr % nx /= dimsizes(1) * dimsizes(2) * 2) then
       print *,'Different dimensions - data file and specified'
       print *,'dimsizes(1)=',dimsizes(1)
       print *,'dimsizes(2)=',dimsizes(2)
       print *,'nx         =',gr%nx
       print *,'(read_CERSAT_data)'
       call exit(1) 
    end if

    ! Read data from drift file
    drnx=dimsizes(1)
    drny=dimsizes(2)
    allocate(drlon(drnx,drny))
    allocate(drlat(drnx,drny))
    allocate(drmer(drnx,drny))
    allocate(drzon(drnx,drny))
    allocate(qflag(drnx,drny))
    call nfw_inq_varid(driftfile, ncid, 'longitude', lon_id)
    !call nfw_get_var_double(driftfile, ncid, lon_id, drlon)
    call cersat_readfield(driftfile, ncid, lon_id, drlon, drnx * drny,1)
    call nfw_inq_varid(driftfile, ncid, 'latitude', lat_id)
    !call nfw_get_var_double(driftfile, ncid, lat_id, drlat)
    call cersat_readfield(driftfile, ncid, lat_id, drlat, drnx * drny,1)
    call nfw_inq_varid(driftfile, ncid, 'zonal_motion', zon_id)
    !call nfw_get_var_double(driftfile, ncid, zon_id, drzon)
    call cersat_readfield(driftfile, ncid, zon_id, drzon, drnx * drny,1)
    call nfw_inq_varid(driftfile, ncid, 'meridional_motion', mer_id)
    !call nfw_get_var_double(driftfile, ncid, mer_id, drmer)
    call cersat_readfield(driftfile, ncid, mer_id, drmer, drnx * drny,1)

    call nfw_get_att_double(driftfile, ncid, zon_id, '_FillValue', fillval)
    call nfw_get_att_double(driftfile, ncid, zon_id, 'scale_factor', scalefac)
    call nfw_get_att_double(driftfile, ncid, zon_id, 'add_offset', addoffset)

    where (abs(drzon - (fillval(1) * scalefac(1) + addoffset(1))) <&
         1e-4 * fillval(1) * scalefac(1) + addoffset(1))
       drzon = gr % undef
    end where

    call nfw_get_att_double(driftfile, ncid, mer_id, '_FillValue', fillval)
    call nfw_get_att_double(driftfile, ncid, mer_id, 'scale_factor', scalefac)
    call nfw_get_att_double(driftfile, ncid, mer_id, 'add_offset', addoffset)

    ! Flag zonal motion for fill values
    where (abs(drmer - (fillval(1) * scalefac(1) + addoffset(1))) <&
         1e-4 * fillval(1) * scalefac(1) + addoffset(1))
       drmer = gr % undef
    end where

    call nfw_inq_varid(driftfile, ncid, 'quality_flag', qua_id)
    call nfw_get_var_int(driftfile, ncid, qua_id, qflag)

    call nfw_close(driftfile, ncid)

    k = 0
    do icomp = 1, 2
       do j = 1, drny ! gr%ny
          do i = 1, drnx ! gr%nx
             k = k + 1

             ! Qualit flag bits - may be signed 
             tmpint = qflag(i,j)
             bit(7) = tmpint/128;  tmpint = tmpint - bit(7)*128 ! Not used
             bit(6) = tmpint/ 64;  tmpint = tmpint - bit(6)* 64 ! Validated using all available info
             bit(5) = tmpint/ 32;  tmpint = tmpint - bit(5)* 32 ! Validated using local consistency
             bit(4) = tmpint/ 16;  tmpint = tmpint - bit(4)* 16 ! Cost function used
             bit(3) = tmpint/  8;  tmpint = tmpint - bit(3)*  8 ! Two identical drift vectors
             bit(2) = tmpint/  4;  tmpint = tmpint - bit(2)*  4 ! SSMI V selected
             bit(1) = tmpint/  2;  tmpint = tmpint - bit(1)*  2 ! SSMI H used
             bit(0) = tmpint/  1;  tmpint = tmpint - bit(0)*  1 ! Quickscat used
             
             valid = qflag(i,j) < 127 ! Intermediate solution until I figure out the byte stuff
             if (icomp==1) then
                data(k)%id = 'VICE'
                data(k)%d = drmer(i,j)*.001 ! Convert to km
                valid = valid .and.  abs( (drmer(i,j)-gr%undef)   / gr%undef)   > 1e-4
             else
                data(k)%id = 'UICE'
                data(k)%d = drzon(i,j)*.001 ! Convert to km
                valid =  valid .and. abs( (drzon(i,j)-gr%undef)   / gr%undef)   > 1e-4
             end if

             if (.not.valid) then
                data(k)%d = gr%undef
             end if

             data(k)%ipiv = i  ! Not really used for ice drift
             data(k)%jpiv = j  ! Not really used for ice drift
             data(k)%i_orig_grid = i ! Used for ice drift
             data(k)%j_orig_grid = j ! Used for ice drift
             data(k)%lat=drlat(i,j)
             data(k)%lon=ang180(real(drlon(i,j)))
             !LB: Data support is assumed = a square grid cell
             !support diameter in meters stored in %a1 (tricky, isn't it ?)
             ! KAL -- hardcoded from data
             data(k)%a1 = 1.4 * 16000.0
             data(k)%ns = 1
             ! To be decided - obs units are meters O(1e4)
             ! CERSAT grid cells are of ~30 km - We assume the errors are
             ! roughly ~15 km
             !KAL data(k)%var = (15)**2
             data(k)%var = var ! fom idrft.hdr specification
             data(k)%depth = 0.0
             data(k)%status = valid
          enddo
       enddo
    enddo
    print*, 'Number of data read:', k, gridpoints(gr)
    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end subroutine read_CERSAT_data


  subroutine cersat_readfield(fname, ncid, varid, v, vlen,Kscale)
    use nfw_mod
    implicit none
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: Kscale
    integer :: vlen
    real(8), intent(out) :: v(vlen)
    
    real(8) :: scale_factor(1)
    real(8) :: offset(1)

    if (Kscale==1) then
       print *, trim(fname)//' scale_factor'
       call nfw_get_att_double(fname, ncid, varid, 'scale_factor', scale_factor)
       call nfw_get_att_double(fname, ncid, varid, 'add_offset', offset)
       call nfw_get_var_double(fname, ncid, varid, v)
       v = v * scale_factor(1) + offset(1)
    else
       call nfw_get_var_double(fname, ncid, varid, v)
    endif
  end subroutine cersat_readfield

end module m_read_CERSAT_data


