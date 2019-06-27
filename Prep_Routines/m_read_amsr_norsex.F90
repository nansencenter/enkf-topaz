! File:          m_read_amsr_norsex.F90
!
! Created:       ???
!
! Last modified: 29/06/2010
!
! Purpose:       Reads ICEC data
!
! Description:   
!
! Modifications:
!                29/06/2010 PS:
!                  - set the maximum ICEC to 0.995 to match the model
!                Prior history:
!                  Not documented.

 module m_read_amsr_norsex
! Reads amsr icec from NERSC (norsex)
! This will only work for northern hemisphere data - some minor corrections
! are needed for SH data

  integer, parameter, private :: STRLEN = 512

contains
  subroutine read_amsr_norsex(fname,gr,data,obstype)
  use mod_grid
  use mod_measurement
  implicit none
  type (grid),      intent(out) :: gr
  character(len=*) ,intent(in)  :: fname
  type(measurement), allocatable, intent(out) :: data(:)
  character(len=5) ,intent(in)  :: obstype

  integer :: i, j,k, rlen
  integer*1, allocatable :: iofldi1(:,:)
  integer*4, allocatable :: iofldi4(:,:)
  real   *4, allocatable, dimension(:,:) :: lon,lat,icec
  logical :: ex(3)

  ! The grid stuff should be made more consistent - KAL
  gr = default_grid
  gr%undef=120.
  gr%nx=608
  gr%ny=896
  gr%order=2
  gr%ux='12.5 km' !Roughly
  gr%uy='12.5 km' !Roughly
  gr%set=.true.
  print '(a,3e14.3)','undef          : ', gr%undef
  print *,' No of gridpoints: ', gridpoints(gr)

  ! Test for input files:
  inquire(exist=ex(1),file='psn12lons_v2.dat')
  inquire(exist=ex(2),file='psn12lats_v2.dat')
  inquire(exist=ex(3),file=trim(fname))

  if (any(.not.ex)) then
     print *,'A file is missing:'
     print *,'File flag: ',ex(1),' - name: psn12lons_v2.dat'
     print *,'File flag: ',ex(2),' - name: psn12lats_v2.dat'
     print *,'File flag: ',ex(3),' - name: '//trim(fname)
     print *,'(read_amsr_norsex)'
     call exit(1)
  end if


  ! Allocate fields and read input data
  allocate(icec   (gr%nx,gr%ny))
  allocate(lon    (gr%nx,gr%ny))
  allocate(lat    (gr%nx,gr%ny))
  allocate(iofldi1(gr%nx,gr%ny))
  allocate(iofldi4(gr%nx,gr%ny))
  allocate(data (gr%nx*gr%ny))

  inquire(iolength=rlen) iofldi4
  open(10,file='psn12lons_v2.dat',status='old',form='unformatted',access='direct',recl=rlen)
  read(10,rec=1) iofldi4 
  close(10)
#if defined (LITTLE_ENDIAN) /* Lon/lat input files are big endian */
  do j=1,gr%ny
  do i=1,gr%nx
     call swapendian2(iofldi4(i,j),4)
  end do
  end do
#endif
  lon = real(iofldi4,4) / 100000.0_4


  inquire(iolength=rlen) iofldi4
  open(10,file='psn12lats_v2.dat',status='old',form='unformatted',access='direct',recl=rlen)
  read(10,rec=1) iofldi4 
  close(10)
#if defined (LITTLE_ENDIAN) /* Lon/lat input files are big endian */
  do j=1,gr%ny
  do i=1,gr%nx
     call swapendian2(iofldi4(i,j),4)
  end do
  end do
#endif
  lat = real(iofldi4, 4) / 100000.0_4

  inquire(iolength=rlen) iofldi1
  open(10,file=trim(fname),status='old',form='unformatted',access='direct',recl=rlen)
  read(10,rec=1) iofldi1
  close(10)

  icec=iofldi1
  where(icec>100)  
     icec = real(gr % undef, 4)
  elsewhere
     icec = icec / 100.0_4
     !LB tighten observed pack ice 
     !where (icec>0.9) icec = 1.0
  end where
  ! PS 25/06/2010 0.995 is the max allowed by the model
  where (0.995 <= icec .and. icec <= 1.0)
     icec = 0.995
  end where


  do j=1,gr%ny
  do i=1,gr%nx
 
     k=(j-1)*gr%nx +i
   
     data(k)%id = obstype
     data(k)%d = icec(i,j)
     data(k)%jpiv = j
     data(k)%ipiv = i
     data(k)%lon=lon(i,j)
     data(k)%lat=lat(i,j)

!LB: Data support is assumed = a square grid cell
!support diameter stored in %a1 (tricky, isn't it ?)
     data(k)%a1 = 12500. *sqrt(2.) ! AMSR-E grid diagonal
     data(k)%ns = 1 ! 1 for obs with a spatial extent

     data(k)%status = .not. undefined(data(k)%d,gr) ! active data
     ! PS 17.06.2010 - increased obs error at the ice edge
     ! data(k)%var = 0.01 ! KAL 10%
     data(k) % var = 0.01d0 + (0.5d0 - abs(0.5d0 - icec(i,j))) ** 2
     data(k) % depth = 0.0
  end do
  end do

  call icec2nc(gr % nx, gr % ny, icec, lon, lat)
end subroutine read_amsr_norsex

subroutine icec2nc(ni, nj, icec, lon, lat)
  use nfw_mod

  integer, intent(in) :: ni
  integer, intent(in) :: nj
  real*4, intent(in) :: icec(ni, nj), lon (ni, nj), lat(ni, nj)

  character(STRLEN) :: fname
  integer :: ncid
  integer :: nij_id(2), icec_id, lon_id, lat_id

  fname = 'icec.nc';
  call nfw_create(fname, nf_clobber, ncid)
  call nfw_def_dim(fname, ncid, 'ni', ni, nij_id(1));
  call nfw_def_dim(fname, ncid, 'nj', nj, nij_id(2));
  call nfw_def_var(fname, ncid,  'icec', nf_float, 2, nij_id, icec_id)
  call nfw_def_var(fname, ncid,  'lon', nf_float, 2, nij_id, lon_id)
  call nfw_def_var(fname, ncid,  'lat', nf_float, 2, nij_id, lat_id)
  call nfw_enddef(fname, ncid)

  call nfw_put_var_real(fname, ncid, icec_id, icec)
  call nfw_put_var_real(fname, ncid, lon_id, lon)
  call nfw_put_var_real(fname, ncid, lat_id, lat)
  call nfw_close(fname, ncid)
end subroutine icec2nc
 
end module m_read_amsr_norsex
