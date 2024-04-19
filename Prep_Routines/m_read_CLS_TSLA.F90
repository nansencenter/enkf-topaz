module m_read_CLS_TSLA

  integer, parameter, private :: STRLEN = 512
  real(8), parameter, private :: RE_MULTIPLE = 0.7d0
  character(*), parameter, private :: RE_FNAME = "re_sla.nc"

contains

  subroutine read_CLS_TSLA(filename, gr, data)
    use mod_measurement
    use mod_grid
    use nfw_mod
    implicit none

    character(*), intent(in) :: filename
    type(grid), intent(inout) :: gr ! CLS measurement grid
    type(measurement), intent(inout) :: data(:)

    integer :: data_ID, track_ID, cycl_ID
    integer :: vNbpoint_ID, vLongitude_ID, vLatitude_ID, vBegindate_ID, vSLA_ID

    ! array dimensions
    integer :: nb, ntracks, ncycles

    ! data arrays
    real(8), allocatable :: vsla(:,:), vlon(:), vlat(:), vbegindate(:,:)
    integer, allocatable :: vnbpoint(:)
    logical, allocatable :: isgood(:,:)

    integer :: ncid
    real(8), dimension(1) :: undef_sla, undef_lat, undef_lon, undef_begindate
    real(8) :: varsat
    integer, dimension(1) :: undef_nbpoint
    integer :: i, j, k, nobs, obsid, sid, age
    real(8), parameter :: EPS = 0.01  ! test for undefined values
    character(STRLEN) :: ftemplate
    character(STRLEN) :: fname
    character(STRLEN) :: fpath
    logical :: ex

    print *, 'read_CLS_TSLA():'

    fpath='./'
    read(filename,'(i7)') age
    nobs = 0
    do sid = 1, 7 ! loop over satellite ID
       select case(sid)
       case(1)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_en*.nc'
          varsat = 0.0009 ! 3 cm for ENVISAT 
          print *, '  ENVISSAT:'
       case(2)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_j1*.nc'
          varsat = 0.0009 ! 3 cm for ENVISAT Jason1
          print *, '  Jason1:'
       case(3)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_j2*.nc'
          varsat = 0.0009 ! 3 cm for ENVISAT Jason2
          print *, '  Jason2:'
       case(4)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_e1*.nc'
          varsat = 0.0075 ! 8.5 cm for e1
          print *, '  ERS1:'
       case(5)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_e2*.nc'
          varsat = 0.0075 ! 8.5 cm for e2
          print *, '  ERS2:'
       case(6)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_tp*.nc'
          varsat = 0.0030 ! 5.5 cm for TOPEX 
          print *, '  TOPEX:'
       case(7)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_g2*.nc'
          varsat = 0.0030 ! GEOSAT
          print *, '  GEOSAT2:'
       end select

       call fname_fromtemplate(ftemplate, fname)
       inquire(file = trim(fname), exist = ex)
       if (.not. ex) then
          cycle
       end if

       ! Reading the observation file of satellite
       call nfw_open(fname, nf_nowrite, ncid)
       call nfw_inq_dimid(fname, ncid, 'Data', data_ID)
       call nfw_inq_dimid(fname, ncid, 'Tracks', track_ID)
       call nfw_inq_dimid(fname, ncid, 'Cycles', cycl_ID)
          
       call nfw_inq_dimlen(fname, ncid, data_ID, nb)
       call nfw_inq_dimlen(fname, ncid, track_ID, ntracks)
       call nfw_inq_dimlen(fname, ncid, cycl_ID, ncycles)
       print '(1x, a, 3i8)', '    dimensions (# obs, # tracks, # cycles):', nb, ntracks, ncycles
       
       allocate(vlon(nb), vlat(nb), vsla(ncycles, nb))
       allocate(vnbpoint(ntracks), vbegindate(ncycles, ntracks))
       allocate(isgood(ncycles, ntracks))
          
       ! Variable ids in netcdf file
       call nfw_inq_varid(fname, ncid, 'Latitudes', vLatitude_ID)
       call nfw_inq_varid(fname, ncid,'Longitudes', vLongitude_ID)
       call nfw_inq_varid(fname, ncid,'BeginDates', vBegindate_ID)
       call nfw_inq_varid(fname, ncid,'NbPoints', vNbpoint_ID)
       call nfw_inq_varid(fname, ncid,'SLA', vSLA_ID)
       
       ! Variable _FillValue attributes
       call nfw_get_att_double(fname, ncid, vLatitude_ID , '_FillValue', undef_lat(1))
       call nfw_get_att_double(fname, ncid, vLongitude_ID , '_FillValue', undef_lon(1))
       call nfw_get_att_double(fname, ncid, vSLA_ID, '_FillValue', undef_sla(1))
       call nfw_get_att_int(fname, ncid, vNbpoint_ID, '_FillValue', undef_nbpoint(1))
       call nfw_get_att_double(fname, ncid,vBegindate_ID, '_FillValue', undef_begindate(1))
       gr % undef = undef_sla(1)
          
       call nfw_get_var_double(fname, ncid, vLongitude_ID, vlon)
       call nfw_get_var_double(fname, ncid, vLatitude_ID, vlat)
       call nfw_get_var_double(fname, ncid, vSLA_ID, vsla)
       !lon = ang180(lon)
       vlon = vlon * 1.e-06
       vlat = vlat * 1.e-06
       print '(1x, a, 2f10.2)', '    range Lon = ', minval(vlon), maxval(vlon)
       print '(1x, a, 2f10.2)', '    range Lat = ', minval(vlat), maxval(vlat)
       print '(1x, a, 2f10.2)', '    range SLA = ', minval(vsla), maxval(vsla)
          
       call nfw_get_var_int(fname, ncid, vNbpoint_ID, vnbpoint)
       call nfw_get_var_double(fname, ncid, vBegindate_ID, vbegindate)
       print '(1x, a, 2i8)', '    range nbpoints = ', minval(vnbpoint), maxval(vnbpoint)
       print *, '    age = ', age
       isgood = .false.
       where (vbegindate /= undef_begindate(1))
          vbegindate = age - floor(vbegindate) - 1
          isgood = .true.
       end where
       print '(3x,a,2G10.3)','  range begin_date (days from assim) = ', &
            minval(pack(vbegindate, isgood)), maxval(pack(vbegindate, isgood))
       call nfw_close(fname, ncid)

       ! Here we set the reference the date with respect to the assimilation
       ! date (0=today, 6=is 6 day old).
       ! Fanf: We assume that the data from the same pass have the same
       ! date=begindate(passnb).
       ! We also assume that envisat, J1 and J2 have similar accuracy, and 
       ! thus use data%var to store the age of the data. Only data that are 
       ! younger than 6 days are retained such that we do not assimilate the
       ! same obs twice.
       do k = 1, ncycles 
          obsid = 0
          do i = 1, ntracks
             do j = 1, vnbpoint(i)
                obsid = obsid + 1
                ! only consider data above -30 of lat 
                if (vlat(obsid) <= -30.0 .or.&
                     vbegindate(k, i) >= 7 .or. vbegindate(k, i) <= -1 .or.&
                     vlon(obsid) == undef_lon(1) .or.&
                     vlat(obsid) == undef_lat(1) .or.&
                     vsla(k, obsid) == undef_sla(1)) then
                   cycle
                end if
                nobs = nobs + 1
                data(nobs) % id = 'TSLA'
                data(nobs) % d = vsla(k, obsid) * 0.001  ! conversion to meters
                data(nobs) % ipiv = -1 ! to be filled
                data(nobs) % jpiv = -1
                data(nobs) % lat = vlat(obsid)
                data(nobs) % lon = ang180(real(vlon(obsid)))
                data(nobs) % a1 = -1.0e10 ! to be filled
                data(nobs) % a2 = -1.0e10
                data(nobs) % a3 = -1.0e10
                data(nobs) % a4 = -1.0e10
                data(nobs) % ns = 0
                data(nobs) % var = varsat
                data(nobs) % date = int(vbegindate(k, i))
                data(nobs) % depth = 0.0
                data(nobs) % status = .true.
             enddo   ! Vnbpoint
          enddo    ! track
       enddo   ! cycle
       print*, '    # of obs read so far = ', nobs
       deallocate(vlat, vlon, vsla, vnbpoint, vbegindate, isgood)
    end do ! satellite id
    gr % nx = nobs
  end subroutine read_CLS_TSLA

 subroutine read_MYO_TSLA(julday,dayinweek, gr, data)
    use mod_measurement
    use mod_grid
    use nfw_mod
    implicit none
!Fanf: this routine assume that you have one seperate file for each day.
!Call prepobs 7 times (for each cycle days with the date and the corrsponding
!day in the cycle

!    character(*), intent(in) :: filename
    character(*), intent(in) :: julday, dayinweek
    type(grid), intent(inout) :: gr ! MYO measurement grid
    type(measurement), intent(inout) :: data(:)

    integer :: time_ID !data_ID, track_ID, cycl_ID
    integer :: vNbpoint_ID, vLongitude_ID, vLatitude_ID, vBegindate_ID, vSLA_ID, vtime_ID

    ! array dimensions
    integer :: nb !, ntracks, ncycles

    ! data arrays
    real(8), allocatable :: vsla(:), vlon(:), vlat(:), vtime(:)!, vbegindate(:,:)
    logical, allocatable :: isgood(:)

    integer :: ncid
    real(8), dimension(1) :: undef_sla, undef_lat, undef_lon, undef_begindate, undef_time
    real(8) :: varsat
    integer, dimension(1) :: undef_nbpoint
    integer :: i, j, k, nobs, obsid, sid, age, idayinweek
    real(8), parameter :: EPS = 0.01  ! test for undefined values
    character(STRLEN) :: ftemplate
    character(STRLEN) :: fname
    character(STRLEN) :: fpath
    logical :: ex

    print *, 'read_MYO_TSLA():'

    fpath='./'
    read(julday,'(i7)') age
    read(dayinweek,'(i2)') idayinweek 
    nobs = 0
    do sid = 1, 17 ! loop over satellite ID
       select case(sid)
       case(1)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_en.nc'
          varsat = 0.0009 ! 3 cm for ENVISAT 
          print *, '  ENVISSAT:'
       case(11)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_enn.nc'
          varsat = 0.0009 ! 3 cm for ENVISAT 
          print *, '  ENVISSAT:'
       case(2)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_j1.nc'
          varsat = 0.0009 ! 3 cm for ENVISAT Jason1
          print *, '  Jason1:'
       case(10)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_j1*.nc'
          varsat = 0.0009 ! 3 cm for ENVISAT Jason1
          print *, '  Jason1:'
       case(13)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_h2*.nc'
          varsat = 0.0016 ! 4 cm for Haiyang2
          print *, '  H2:'
       case(3)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_j2*.nc'
          varsat = 0.0009 ! 3 cm for ENVISAT Jason2
          print *, '  Jason2:'
       case(4)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_e1*.nc'
          varsat = 0.0075 ! 8.5 cm for e1
          print *, '  ERS1:'
       case(5)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_e2*.nc'
          varsat = 0.0075 ! 8.5 cm for e2
          print *, '  ERS2:'
       case(6)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_tp*.nc'
          varsat = 0.0030 ! 5.5 cm for TOPEX 
          print *, '  TOPEX:'
       case(7)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_g2*.nc'
          varsat = 0.0030 ! 5.5 cm for GEOSAT
          print *, '  GEOSAT2:'
       case(8)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_c2*.nc'
          varsat = 0.0030 ! CRYOSAT-2
          print *, '  CRYOSAT2:'
       case(9)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_a*.nc'
          varsat = 0.0009 ! ALTIKA
          print *, '  ALTIKA:'
       case(14)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_j3*.nc'
          varsat = 0.0009 !  J3
          print *, '  J3:'
       case(15)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_s3a.nc'
          varsat = 0.0009 !  S3a/S3b
          print *, '  S3a:'
       case(16)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_s3b.nc'
          varsat = 0.0009 !  S3a/S3b
          print *, '  S3b:'
       case(17)
          ftemplate = trim(fpath)//'sla_'//trim(julday)//'_s6*.nc'
          varsat = 0.0009 !  S6a
          print *, '  S6:'
       end select
       call fname_fromtemplate(ftemplate, fname)

       inquire(file = trim(fname), exist = ex)
       if (.not. ex) then
          cycle
       end if
     

       ! Reading the observation file of satellite
       call nfw_open(fname, nf_nowrite, ncid)
       call nfw_inq_dimid(fname, ncid, 'time', time_ID)
          
       call nfw_inq_dimlen(fname, ncid, time_ID, nb)
       print '(1x, a, i8)', '    dimensions (# obs):', nb !, ntracks, ncycles
       
       allocate(vlon(nb), vlat(nb), vsla(nb), vtime(nb))
       allocate(isgood(nb))
          
       ! Variable ids in netcdf file
       call nfw_inq_varid(fname, ncid, 'latitude', vLatitude_ID)
       call nfw_inq_varid(fname, ncid,'longitude', vLongitude_ID)
       call nfw_inq_varid(fname, ncid, 'time', vtime_ID)
!       call nfw_inq_varid(fname, ncid,'SLA', vSLA_ID)
       ! modified by Jiping on 09-June-2017
       call nfw_inq_varid(fname, ncid,'sla_unfiltered', vSLA_ID)
       ! used by nrt product
       !call nfw_inq_varid(fname, ncid,'sla_filtered', vSLA_ID)

       ! Variable _FillValue attributes
       call nfw_get_att_double(fname, ncid, vSLA_ID, '_FillValue', undef_sla(1))
       gr % undef = undef_sla(1)
          
       call nfw_get_var_double(fname, ncid, vLongitude_ID, vlon)
       call nfw_get_var_double(fname, ncid, vLatitude_ID, vlat)
       call nfw_get_var_double(fname, ncid, vSLA_ID, vsla)
       call nfw_get_var_double(fname, ncid, vtime_ID, vtime)
       !lon = ang180(lon)
       vlon = vlon * 1.e-06
       vlat = vlat * 1.e-06
       print '(1x, a, 2f10.2)', '    range Lon = ', minval(vlon), maxval(vlon)
       print '(1x, a, 2f10.2)', '    range Lat = ', minval(vlat), maxval(vlat)
       print '(1x, a, 2f10.2)', '    range SLA = ', minval(vsla), maxval(vsla)
          
       print *, '    age = ', age
       print '(3x,a,G10.3)','  Days before assim = ', idayinweek
       call nfw_close(fname, ncid)

       ! Here we set the reference the date with respect to the assimilation
       ! date (0=today, 6=is 6 day old).
       ! Fanf: We assume that the data from the same pass have the same
       ! date=begindate(passnb).
       ! We also assume that envisat, J1 and J2 have similar accuracy, and 
       ! thus use data%var to store the age of the data. Only data that are 
       ! younger than 6 days are retained such that we do not assimilate the
       ! same obs twice.
       do k = 1, nb 
          ! only consider data above -30 of lat 
#if defined (Rio22)
          ! add considering data south of 80 of lat 
          if (vlat(k) <= -30.0 .or. vlat(k) > 85.0 .or. & 
               vsla(k) == undef_sla(1)) then
             cycle
          end if
#else
          ! add considering data south of 80 of lat 
          if (vlat(k) <= -30.0 .or. &
              (vlat(k) > 82.0 .and. (vlon(k)<120 .or. vlon(k)>320)) .or. &
              (vlat(k) > 78.0 .and. (vlon(k)>120 .and. vlon(k)<320)) .or. &
               vsla(k) == undef_sla(1)) then
             cycle
          end if
#endif
          nobs = nobs + 1
          data(nobs) % id = 'TSLA'
          data(nobs) % d = real(vsla(k)*0.001)  ! conversion to meters
          data(nobs) % ipiv = -1 ! to be filled
          data(nobs) % jpiv = -1
          data(nobs) % lat = vlat(k)
          data(nobs) % lon = ang180(real(vlon(k)))
          data(nobs) % a1 = -1.0e10 ! to be filled
          data(nobs) % a2 = -1.0e10
          data(nobs) % a3 = -1.0e10
          data(nobs) % a4 = -1.0e10
          data(nobs) % ns = 0
          data(nobs) % var = varsat
          data(nobs) % date = idayinweek
          data(nobs) % depth = 0.0
          data(nobs) % status = .true.
       enddo   ! cycle
       print*, '    # of obs read so far = ', nobs,k 
       deallocate(vlat, vlon, vsla, vtime, isgood)
    end do ! satellite id
    gr % nx = nobs
  end subroutine read_MYO_TSLA

  subroutine set_re_TSLA(nrobs, obs, nx, ny, modlon, modlat)
    use mod_measurement
    use nfw_mod

    integer, intent(in) :: nrobs
    type(measurement), dimension(nrobs), intent(inout) :: obs
    integer, intent(in) :: nx, ny
    real, dimension(nx, ny), intent(in)  ::  modlon, modlat

    integer :: ncid, reid
    real*8, dimension(nx, ny) :: re
    real :: reo
    integer :: o

    print *, '  reading representation error from "', trim(RE_FNAME), '"'

    call nfw_open(RE_FNAME, nf_nowrite, ncid)
    call nfw_inq_varid(RE_FNAME, ncid, 're_sla', reid)
    call nfw_get_var_double(RE_FNAME, ncid, reid, re)
    call nfw_close(RE_FNAME, ncid)
    
    do o = 1, nrobs
       reo = re(obs(o) % ipiv, obs(o) % jpiv)
       if (reo < 0 .or. reo > 1.0d5) then
          cycle
       end if
       ! PS 1.4.2010 Increased the multiple for representation error from
       ! 0.3 to 0.5 - it seems that with 0.3 it wants to do more in the Gulf
       ! Stream region than the model can sustain.
       ! PS June 2010 - further increased the multiple to 0.7.
       obs(o) % var = obs(o) % var + RE_MULTIPLE * reo
    end do
  end subroutine set_re_TSLA


  subroutine fname_fromtemplate(ftemplate, fname)
    character(*), intent(in) :: ftemplate
    character(*), intent(inout) :: fname

    character(STRLEN) :: command ! (there may be a limit of 80 on some systems)
    integer :: ios

    command = 'ls '//trim(ftemplate)//' 2> /dev/null > tsla_files.txt'
    print *, trim(command)
    call system(trim(command));

    open(10, file = 'tsla_files.txt')
    read(10, fmt = '(a)', iostat = ios) fname
    close(10)
    if (ios /= 0) then
       fname = ""
    end if
  end subroutine fname_fromtemplate

end module m_read_CLS_TSLA
