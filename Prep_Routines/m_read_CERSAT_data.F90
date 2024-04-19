module m_read_CERSAT_data

contains



  subroutine read_CERSAT_data_rep(driftfile, gr, data, numdata, var,offset)
    use nfw_mod
    use mod_measurement
    use mod_grid
    use m_spherdist
    use m_rk2
    implicit none

    real,   parameter :: undef=-1e14
    character(*), intent(in) :: driftfile
    integer, intent(in) :: numdata
    type(measurement), dimension(numdata) :: data
    type(grid), intent(in) :: gr
    real, intent(in) :: var
    character(len=1), intent(in) :: offset

    integer :: dimids(4)
    integer , dimension(4) :: dimsizes

    integer :: lon_id, lat_id, zon_id, mer_id, qua_id
    real, dimension(:,:), allocatable :: drlon, drlat, drmer, drzon
    real, dimension(:,:), allocatable :: rdx,rdy      ! Resolution in observatoin grid
    real, dimension(:,:), allocatable :: DDX,DDY      ! Resolution in observatoin grid
    real, dimension(:),   allocatable :: tmplon, tmplat
    real, dimension(:,:), allocatable :: tmpmer, tmpzon
    integer, dimension(:,:), allocatable :: Imask

    integer :: ncid, varid
    real, dimension(1) :: scalefac, fillval, addoffset

    integer :: i,j,k,icomp,idate
    integer :: drnx, drny
    logical :: valid
    
    !  positions
    real, dimension(:,:), allocatable :: x0,y0,x,y
    ! RK velocity estimates
    real, dimension(:,:,:), allocatable :: urk,vrk
    real                              :: delt
    integer                           :: nstep,istep
    !  new positions
    real                              :: lon0,lat0,lon1,lat1,lon2,lat2
    real                              :: w0,w1
    integer                           :: tmpi0,tmpi1,tmpj0,tmpj1
    real                              :: lon_1,lat_1,lon_2,lat_2

    read(offset,'(I1)') idate 
    ! Get dimensions of drift file
    call nfw_open(driftfile, nf_nowrite, ncid)
    call nfw_inq_varid(driftfile, ncid,'eastward_sea_ice_velocity', &
         varid)
    call nfw_inq_vardimid(driftfile, ncid, varid, dimids)
    do i = 1, 4
       call nfw_inq_dimlen(driftfile, ncid, dimids(i), dimsizes(i))
    end do

    print *,'numdata=',numdata
    ! Which should match numdata dimension 
    ! NB !!! Mult by 2 for two vector components
    if (2 * dimsizes(1) * dimsizes(2) /= numdata) then 
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
    allocate(tmpmer(drnx,drny),tmpzon(drnx,drny),Imask(drnx,drny))
    allocate(drmer(drnx,drny))
    allocate(drzon(drnx,drny))
    allocate(rdx(drnx,drny),rdy(drnx,drny))
    allocate(DDX(drnx,drny),DDY(drnx,drny))
    call nfw_inq_varid(driftfile, ncid, 'longitude', lon_id)
    call cersat_readfield(driftfile, ncid, lon_id, tmplon, drnx,0)
    call nfw_inq_varid(driftfile, ncid, 'latitude', lat_id)
    call cersat_readfield(driftfile, ncid, lat_id, tmplat, drny,0)

    do i=1, dimsizes(1)
       drlat(i,:)=tmplat
    end do
    do i=1, dimsizes(2)
       do j=1,dimsizes(1)
          !drlon(j,i)=ang180(real(tmplon(j)))
          drlon(j,i)=tmplon(j)
       end do
    end do
    ! derive the grid resolutions
    do j=1,drny
       do i=2,drnx-1
          rdx(i,j)=(spherdist(drlon(i-1,j),drlat(i-1,j),drlon(i,j),drlat(i,j))+  &
                   spherdist(drlon(i,j),drlat(i,j),drlon(i+1,j),drlat(i+1,j)))/2
       end do
    end do
    do j=2,drny-1
       do i=1,drnx
          rdy(i,j)=(spherdist(drlon(i,j-1),drlat(i,j-1),drlon(i,j),drlat(i,j))+  &
                   spherdist(drlon(i,j),drlat(i,j),drlon(i,j+1),drlat(i,j+1)))/2
       end do
    end do
    do j=1,drny
       rdx(1,j)=spherdist(drlon(1,j),drlat(1,j),drlon(2,j),drlat(2,j))
       rdx(drnx,j)=spherdist(drlon(drnx,j),drlat(drnx,j),drlon(drnx-1,j),drlat(drnx-1,j))
    end do
    do i=1,drnx
       rdy(i,1)=spherdist(drlon(i,1),drlat(i,1),drlon(i,2),drlat(i,2))
       rdy(i,drny)=spherdist(drlon(i,drny),drlat(i,drny),drlon(i,drny-1),drlat(i,drny-1))
    end do

    call nfw_inq_varid(driftfile, ncid, 'eastward_sea_ice_velocity',zon_id)
    call cersat_readfield(driftfile, ncid, zon_id, tmpzon, drnx * drny,1)
    call nfw_inq_varid(driftfile, ncid, 'northward_sea_ice_velocity', mer_id)
    call cersat_readfield(driftfile, ncid, mer_id, tmpmer, drnx * drny,1)

    call nfw_get_att_double(driftfile,ncid, zon_id,'_FillValue', fillval)
    call nfw_get_att_double(driftfile,ncid, zon_id,'scale_factor',scalefac)
    call nfw_get_att_double(driftfile,ncid, zon_id, 'add_offset', addoffset)

    Imask=1
    where (abs(tmpzon - (fillval(1) * scalefac(1) + addoffset(1))) < &
         1e-4 * abs(fillval(1) * scalefac(1) + addoffset(1)))
       tmpzon = 0.
       tmpmer =0. 
       Imask=0
    end where

    call nfw_get_att_double(driftfile, ncid, mer_id, '_FillValue', fillval)
    call nfw_get_att_double(driftfile, ncid, mer_id, 'scale_factor', scalefac)
    call nfw_get_att_double(driftfile, ncid, mer_id, 'add_offset', addoffset)

    call nfw_close(driftfile, ncid)

    where (abs(tmpmer - (fillval(1) * scalefac(1) + addoffset(1))) < &
         1e-4 * abs(fillval(1) * scalefac(1) + addoffset(1)))
       tmpzon = 0.
       tmpmer =0. 
       Imask=0
    end where
     
    ! skip the observations at north of 85N
    where (drlat>=85)
       tmpzon = 0.
       tmpmer =0. 
       Imask=0
    end where

    !print *, 'tmpzon: ', tmpzon(2,:,1)
    !print *, 'tmpmer: ', tmpmer(2,:,1)
    allocate(x0   (drnx,drny))
    allocate(y0   (drnx,drny))
    allocate(x    (drnx,drny))
    allocate(y    (drnx,drny))
    allocate(urk  (drnx,drny,2))
    allocate(vrk  (drnx,drny,2))
    x0=undef
    y0=undef
    do i=2,drnx-1
       do j=2,drny-1
          if(Imask(i,j)==1) then
             x0(i,j)=i
             y0(i,j)=j
          endif
       end do
    end do
   
    ! accumulate the drift DX/DY using 
    ! Subroutine advances drift over obe day using second order RK
    !subroutine rk2(u,v,scpx,scpy,nx,ny,x,y,drnx,drny,delt,undef)

    ! Start the time loop from start to end time. We use 6 hours
    ! as the RK(2) time step
    x=x0; y=y0
    delt=86400/4.
    !nstep=2/(delt/86400.)    ! drifting for two days
    nstep=3
    do istep=1,nstep
       urk(:,:,1)=tmpzon(:,:)
       urk(:,:,2)=tmpzon(:,:)
       vrk(:,:,1)=tmpmer(:,:)
       vrk(:,:,2)=tmpmer(:,:)
       !call rk2(urk,vrk,rdx,rdy,drnx,drny,x,y,drnx,drny,delt,undef);
       ! add the limit for the maximal drift where near the boundary?
       call rk2(urk,vrk,rdx,rdy,drnx,drny,x,y,drnx,drny,delt,undef);
       ! filter some odd results
       do i=1,drnx
          do j=1,drny
             if (x(i,j)>=drnx) then
                x(i,j)=undef; 
                y(i,j)=undef; 
                urk(i,j,:)=0; 
                vrk(i,j,:)=0; 
             elseif (x(i,j)<=1) then
                x(i,j)=undef;
                y(i,j)=undef;
                urk(i,j,:)=0; 
                vrk(i,j,:)=0; 
             endif
             if (y(i,j)>=drny) then
                y(i,j)=undef;
                x(i,j)=undef;
                urk(i,j,:)=0; 
                vrk(i,j,:)=0; 
             elseif (y(i,j)<=1) then
                y(i,j)=undef;
                x(i,j)=undef;
                urk(i,j,:)=0; 
                vrk(i,j,:)=0; 
             endif
             if (abs(x(i,j)-undef)<0.001.or.abs(y(i,j)-undef)<0.001.or. &
                 abs(x0(i,j)-undef)<0.001.or. abs(y0(i,j)-undef)<0.001) then
                y(i,j)=undef;
                x(i,j)=undef;
             endif
          end do
       end do

    end do
    print *
    print *,'Displacement max in grid coordinates:'
    print *,'x:',maxval(abs(x-x0))
    print *,'y:',maxval(abs(y-y0))

    do i=1,drnx
       do j=1,drny
          if (abs(x(i,j)-x0(i,j))<0.001.and.abs(y(i,j)-y0(i,j))<0.001) then
             DDX(i,j)=gr%undef
             DDY(i,j)=gr%undef
          elseif (Imask(i,j)==0.or.abs(x(i,j)-undef)<0.001 &
             .or.abs(y(i,j)-undef)<0.001) then
             DDX(i,j)=gr%undef
             DDY(i,j)=gr%undef
          else
             lon0=drlon(i,j)
             lat0=drlat(i,j)
             tmpi0=max(1,floor(x(i,j)))
             tmpi1=min(drnx,ceiling(x(i,j)))
             w0=x(i,j)-real(tmpi0)
             tmpj0=max(1,floor(y(i,j)))
             tmpj1=min(drny,ceiling(y(i,j)))
             w1=y(i,j)-real(tmpj0)  
             lat1=(1-w1)*((1-w0)*drlat(tmpi0,tmpj0)+w0*drlat(tmpi1,tmpj0)) + &
                   w1*(drlat(tmpi0,tmpj1)*(1-w0)+drlat(tmpi1,tmpj1)*w0)

             lon_1=drlon(tmpi0,tmpj0)
             lon_2=drlon(tmpi0,tmpj1)
             if (lon_2>lon_1+180) lon_2=lon_2-360
             if (lon_2<lon_1-180) lon_2=lon_2+360
             lon1=lon_1*(1-w1)+w1*lon_2
             lon_1=drlon(tmpi1,tmpj0)
             lon_2=drlon(tmpi1,tmpj1)
             if (lon_2>lon_1+180) lon_2=lon_2-360
             if (lon_2<lon_1-180) lon_2=lon_2+360
             lon2=lon_1*(1-w1)+w1*lon_2
             if (lon2>lon1+180) lon2=lon2-360
             if (lon2<lon1-180) lon2=lon2+360
             lon_1=(1-w0)*lon1+w0*lon2
             lon1=lon_1

             DDX(i,j)=spherdist(lon0,lat0,lon1,lat0)*sign(1.,lon1-lon0)
             DDY(i,j)=spherdist(lon0,lat0,lon0,lat1)*sign(1.,lat1-lat0)

          endif
          
       end do
    end do

    k = 0
    do icomp = 1, 2
       do j = 1, drny ! gr%ny
          !do i = 1, drnx    ! gr%nx
          ! thinning observations
          do i = 1, drnx,2 ! gr%nx
             k = k + 1
             valid=.true.
             if (icomp==1) then
                data(k)%id = 'DX'//offset
                data(k)%d = DDX(i,j)*0.001 ! Convert to km
                valid = valid .and.  abs( (DDX(i,j)-gr%undef)   / gr%undef)   > 1e-4
             else
                data(k)%id = 'DY'//offset
                data(k)%d = DDY(i,j)*0.001 ! Convert to km
                valid =  valid .and. abs( (DDY(i,j)-gr%undef)   / gr%undef)   > 1e-4
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

    deallocate(drlon,drlat,tmplon,tmplat,drmer,drzon)
    deallocate(tmpmer,tmpzon)
    deallocate(rdx,rdy)
    deallocate(DDX,DDY)
    deallocate(x0,y0,x,y)
    deallocate(urk,vrk)
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


