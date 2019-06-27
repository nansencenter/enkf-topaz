module m_read_CLS_SST_grid
 ! Reads the CLS SST NetCDF dimensions

contains
subroutine read_CLS_SST_grid(filename,gr)
  !use mod_dimensions
  use mod_grid
  use netcdf
  use m_nf90_err
  !use nfw_mod
  implicit none

  character(len=80), intent(in) :: filename
  type(grid),        intent(out) :: gr

!dimension ids
  integer :: NbLatitudes_ID, NbLongitudes_ID, LatLon_ID


! Array dimensions
  integer :: LatLon, NbLatitudes, NbLongitudes
  real, allocatable :: latlon0(:), dlatlon(:)

!variables ids
  integer :: vLatLonMin_ID, vLatLonStep_ID

  integer :: ncid

  gr = default_grid

! Open file
!filename='sst_topaz_19510.nc'
  call nf90_err(NF90_OPEN(trim(filename),NF90_NOCLOBBER,ncid))
  !call nfw_open(trim(filename), nf_nowrite, ncid)

! Get dimension id in netcdf file ...
  call nf90_err(nf90_Inq_Dimid(ncid,'LatLon',LatLon_ID))
  call nf90_err(nf90_Inq_Dimid(ncid,'NbLatitudes',NbLatitudes_ID))
  call nf90_err(nf90_Inq_Dimid(ncid,'NbLongitudes',NbLongitudes_ID))
! Get dimension length from id
  call nf90_err(nf90_Inquire_Dimension(ncid,LatLon_ID,len=LatLon))
  call nf90_err(nf90_Inquire_Dimension(ncid,NbLatitudes_ID,len=NbLatitudes))
  call nf90_err(nf90_Inquire_Dimension(ncid,NbLongitudes_ID,len=NbLongitudes))
!  call nf90_err(nf90_Inquire_Dimension(ncid,GridDepth_ID,len=GridDepth))
  print*, 'Dimensions:', NbLatitudes, NbLongitudes, LatLon

  allocate(latlon0(LatLon))  ! Grid origin coordinates
  allocate(dlatlon(LatLon))  ! dx and dy

! Variable ids in netcdf file
  call nf90_err(nf90_inq_varid(ncid,'LatLonMin' ,vLatLonMin_ID),'LatLonMin')
  call nf90_err(nf90_inq_varid(ncid,'LatLonStep' ,vLatLonStep_ID),'LatLonStep')

! Variables in NetCDF file
  call nf90_err(nf90_get_var(ncid,vLatLonMin_ID     ,latlon0))
  print *, 'Grid Origin ', latlon0
  call nf90_err(nf90_get_var(ncid,vLatLonStep_ID    ,dlatlon))
  print *, 'Grid Size ', dlatlon

  gr%nx=NbLatitudes
  gr%ny=NbLongitudes
  gr%x0=   latlon0(1)
  gr%y0=   latlon0(2)
!  gr%dx=   0.179
  gr%dx=   dlatlon(1)
  gr%dy=   dlatlon(2)
  gr%reg = .true.
  gr%order = 2
  gr%ux = 'deg'
  gr%uy = 'deg'
  gr%set = .true.

  deallocate(latlon0,dlatlon)

  end subroutine read_CLS_SST_grid
end module m_read_CLS_SST_grid

