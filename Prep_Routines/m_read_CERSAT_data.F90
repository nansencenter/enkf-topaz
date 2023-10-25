module m_read_CERSAT_data

contains

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
    call cersat_readfield(driftfile, ncid, lon_id, drlon, drnx * drny)
    call nfw_inq_varid(driftfile, ncid, 'latitude', lat_id)
    !call nfw_get_var_double(driftfile, ncid, lat_id, drlat)
    call cersat_readfield(driftfile, ncid, lat_id, drlat, drnx * drny)
    call nfw_inq_varid(driftfile, ncid, 'zonal_motion', zon_id)
    !call nfw_get_var_double(driftfile, ncid, zon_id, drzon)
    call cersat_readfield(driftfile, ncid, zon_id, drzon, drnx * drny)
    call nfw_inq_varid(driftfile, ncid, 'meridional_motion', mer_id)
    !call nfw_get_var_double(driftfile, ncid, mer_id, drmer)
    call cersat_readfield(driftfile, ncid, mer_id, drmer, drnx * drny)

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
#if defined (CERSAT_QF)
             tmpf=exp((qflag(i,j)-20)/60)
             data(k)%var = var*(tmpf/2)**2 ! relax by the QF specification
#else
             data(k)%var = var ! fom idrft.hdr specification
#endif
             data(k)%depth = 0.0
             data(k)%status = valid
          enddo
       enddo
    enddo
    print*, 'Number of data read:', k, gridpoints(gr)

    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print *,'!!!!!!!!! Adjust obs errors  !!!!!!!!!!!!!!!!!!!'
    print *,'!!!!!!!Use qflag in valid as well!!!!!!!!!!!!!!!'
    print *,'!!!!!!!!!!CHECK use of qflag !!!!!!!!!!!!!!!!!!!'
    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end subroutine read_CERSAT_data


  subroutine cersat_readfield(fname, ncid, varid, v, vlen)
    use nfw_mod
    implicit none

    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer :: vlen
    real(8), intent(out) :: v(vlen)
    
    real(8) :: scale_factor(1)
    real(8) :: offset(1)

    call nfw_get_att_double(fname, ncid, varid, 'scale_factor', scale_factor)
    call nfw_get_att_double(fname, ncid, varid, 'add_offset', offset)
    call nfw_get_var_double(fname, ncid, varid, v)
    v = v * scale_factor(1) + offset(1)
  end subroutine cersat_readfield

end module m_read_CERSAT_data


