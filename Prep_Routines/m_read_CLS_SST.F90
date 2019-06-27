module m_read_CLS_SST
! Reads CLS SST data after having read the grid in read_CLS_SST_grid
  contains

  subroutine read_CLS_SST(fname,gr,data)
  use mod_measurement
!  use mod_dimensions
  use mod_grid
  use m_spherdist
  use netcdf
  use m_nf90_err
  implicit none

  type (measurement),  intent(inout) :: data(:)
  type (grid),         intent(inout) :: gr ! CLS measurement grid
  character(len=80),   intent(in) :: fname

!dimension ids
  integer :: NbLatitudes_ID, NbLongitudes_ID, LatLon_ID

! Variable ids
  integer :: vNbLatitudes_ID, vNbLongitudes_ID, vGrid0001_ID, vGrid0004_ID

! Array dimensions
  integer :: LatLon, NbLatitudes, NbLongitudes

! Data arrays
  real,allocatable :: sst(:,:), sst_var(:,:), lon(:),lat(:)!, latlon0(:), dlatlon(:)

! utilitary
  integer ncid, ijmax(2)
  real undef,undef_lat, undef_lon, undef_var
  integer i, j,k
  logical valid
  real, parameter :: eps = 0.01  ! test for undefined values

! Open file
!  filename='sst_topaz_19510.nc'
  call nf90_err(NF90_OPEN(trim(fname),NF90_NOCLOBBER,ncid))

! Get dimension id in netcdf file ...
  call nf90_err(nf90_Inq_Dimid(ncid,'LatLon',LatLon_ID))
  call nf90_err(nf90_Inq_Dimid(ncid,'NbLatitudes',NbLatitudes_ID))
  call nf90_err(nf90_Inq_Dimid(ncid,'NbLongitudes',NbLongitudes_ID))

! Get dimension length from id
  call nf90_err(nf90_Inquire_Dimension(ncid,LatLon_ID,len=LatLon))
  call nf90_err(nf90_Inquire_Dimension(ncid,NbLatitudes_ID,len=NbLatitudes))
  call nf90_err(nf90_Inquire_Dimension(ncid,NbLongitudes_ID,len=NbLongitudes))
  print*, 'Dimensions:', NbLatitudes, NbLongitudes, LatLon

! State which variable you want here.. Available vars are shown when you do
! "ncdump -h " on the netcdf file. This is for SST and SDEV of SST
  allocate(lon(NbLongitudes))
  allocate(lat(NbLatitudes))
!  allocate(latlon0(LatLon))
!  allocate(dlatlon(LatLon))
  allocate(sst_var(NbLatitudes,NbLongitudes))
  allocate(sst(NbLatitudes,NbLongitudes))

! Variable ids in netcdf file
  call nf90_err(nf90_inq_varid(ncid,'NbLatitudes' ,vNbLatitudes_ID),'NbLatitudes')
  call nf90_err(nf90_inq_varid(ncid,'NbLongitudes' ,vNbLongitudes_ID),'NbLongitudes')
  call nf90_err(nf90_inq_varid(ncid,'Grid_0001' ,vGrid0001_ID),'Grid_0001')
  call nf90_err(nf90_inq_varid(ncid,'Grid_0004' ,vGrid0004_ID),'Grid_0004')

! Variable _FillValue attributes
  call nf90_err(nf90_get_att(ncid,vNbLatitudes_ID , '_FillValue',undef_lat))
  call nf90_err(nf90_get_att(ncid,vNbLongitudes_ID ,'_FillValue',undef_lon))
  call nf90_err(nf90_get_att(ncid,vGrid0001_ID ,   '_FillValue',undef))
  call nf90_err(nf90_get_att(ncid,vGrid0004_ID ,   '_FillValue',undef_var))
  print*, 'Undefined values are ', undef_lat, undef_lon, undef, undef_var
  gr%undef = undef

! actual variable values (for dimensions of var -- see ncdump, or improve this program)
! NB: note that index dimensions are different between fortran and C internals. 
! "ncdump" gives C internals.
  print *,'test'
  call nf90_err(nf90_get_var(ncid,vNbLongitudes_ID  ,lon))
  print *,'Range Lon', minval(lon), maxval(lon)
  call nf90_err(nf90_get_var(ncid,vNbLatitudes_ID   ,lat))
  print *,'Range Lat', minval(lat), maxval(lat)
  call nf90_err(nf90_get_var(ncid,vGrid0001_ID      ,sst))
  print *,'Range SST', minval(sst), maxval(sst)
  call nf90_err(nf90_get_var(ncid,vGrid0004_ID      ,sst_var))
  print *,'Range Std. Dev.', minval(sst_var), maxval(sst_var)

!  print*, 'Latitudes'
!  print '(12f8.2)',lat
!  print*, 'Longitudes'
!  print '(12f8.2)',lon
  print '(4a10)','Lat','Lon','Temp[C]','Error[C]'
!   print*,lat,lon,temp(i),err_temp(i),saln(i),err_saln(i)
!   write(13,*)lat,lon,saln(i),err_saln(i),depth(i)
!   write(14,*)lat,lon,temp(i),err_temp(i),depth(i)
  ijmax = minloc(sst)
  do i=ijmax(1)-5, ijmax(1)+5
    j = ijmax(2)
    print '(4f10.3)', lat(i), lon(j), sst(i,j), sst_var(i,j)
  enddo

  call nf90_err (nf90_close(ncid))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fill the data(:) vector

  do j=1,NbLongitudes ! gr%ny
  do i=1,NbLatitudes ! gr%nx
     k=(j-1)*gr%nx+i

     data(k)%id = 'SST'
     data(k)%d = sst(i,j)

     data(k)%ipiv = i
     data(k)%jpiv = j

     data(k)%lat=lat(i)
     data(k)%lon=lon(j)

!LB: Data support is assumed = a square grid cell
!support diameter in meters stored in %a1 (tricky, isn't it ?)
     data(k)%a1 = spherdist(lon(j)-.5*gr%dx,lat(i)-.5*gr%dy, &
                            lon(j)+.5*gr%dx,lat(i)+.5*gr%dy)
     data(k)%ns = 1
 
     data(k)%var = sst_var(i,j)  ! corrected: variance is provided

     data(k)%depth = 0.0

     valid =   (abs( (lon(j)-undef_lon) / undef_lon ) > eps  & 
         .and.  abs( (lat(i)-undef_lat) / undef_lat ) > eps  &
         .and.  abs( (sst(i,j)-undef)   / undef )     > eps  &
         .and.  abs( (sst_var(i,j)-undef_var)/undef_var) > eps &
         .and.  sst_var(i,j)> 0 &
         .and.  sst_var(i,j)< 16 )   ! Sdev too high => perturbations too large

     data(k)%status = valid

  enddo
  enddo
  print*, 'Number of data read:', k, gridpoints(gr)

  deallocate(lat,lon,sst,sst_var)

end subroutine read_CLS_SST
end module m_read_CLS_SST
