 module m_read_CLS_data
! Reads SLA and SST data from CLS, Toulouse, France
! Files are given as .asc (lat,lon,data) 
! The data points are surface data and therefore the data(k)%depths=0
! This subroutine also prepares the gridd which the data
! Is provided on.

contains

  subroutine read_CLS_data(fname,obstype,dformat,gr,form,data,factor,var)
  use mod_measurement
  use mod_grid
  use m_spherdist
  implicit none

  type (measurement), intent(inout)  :: data(:)
  type (grid), intent(in)           :: gr ! CLS measurement grid
  real, intent(in) :: factor, var

  character(len=80), intent(in) :: fname,dformat
  character(len=3), intent(in)::form
  character(len=*), intent(in)::obstype
  integer :: k, telli, tellj, nrdat, dum
  logical :: ex, found, fleeting
  real :: lon, lat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read  observation file
  if (trim(form) == '0') stop 'read_CLS_data: illegal format '
  inquire (file=fname, exist=ex)
  if (.not. ex) then
     print *, 'read_CLS_data: file ', fname, ' not found.'
     stop
  end if

!std = 0.0; lat = 0.0; lon = 0.0

!!! Find out if data column is type integer or not

  found = .false.
  found = ((scan(dformat,'i') > 0) .or. (scan(dformat,'I') > 0))
  fleeting = .not. found
 
  open (10, file=fname, form='formatted', status='old')

  telli=1  
  tellj=1

  do k = 1, gridpoints(gr) 
     data(k)%id = obstype

     if (fleeting) then  ! Data column floating point
       read (10,dformat,end=999,err=999) lat, lon, data(k)%d
     else  ! Data column integer valued
       read (10,dformat,end=999,err=999) lat, lon, dum
       data(k)%d = real(dum)
     end if
!      print*,'lat',lat,'lon', lon,'data',data(k)%d

!NBNBN Avoid sla data in region 3S to 3N (due to strange Ifremer mean ssh in this region):

!     if (trim(data(k)%id) == 'ssh' .or. trim(data(k)%id) == 'SSH') then
!        if ((lat.ge.-3.0).and.(lat.le.3.0)) then
!           data(k)%d = 999.9
!        endif
!     endif

     if (.not. undefined(data(k)%d,gr)) then
         data(k)%d = data(k)%d*factor ! Convert to proper units
     end if

     data(k)%jpiv = telli
     data(k)%ipiv = tellj
!        iloop(k) = telli
!        jloop(k) = tellj
        

     telli = telli + 1

     if (telli > gr%ny) then
        tellj=tellj+1
        telli = 1
     endif

     data(k)%lon=lon
     data(k)%lat=lat

!LB: Data support is assumed = a square grid cell
!support diameter stored in %a1 (tricky, isn't it ?)
     data(k)%a1 = spherdist(lon-.5*gr%dx,lat-.5*gr%dy,lon+.5*gr%dx,lat+.5*gr%dy)
     data(k)%ns = 1
     data(k)%status = .not. undefined(data(k)%d,gr) ! active data

     if (trim(obstype) == 'SST') then
        data(k) % status = data(k) % status .and.&
             abs(data(k) % d + 1.8d0) > 1.0d-6
     end if
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! In the case of SSH data the var parameter is given for each data point !!!!
!
     if (trim(data(k)%id) == 'ssh' .or. trim(data(k)%id) == 'SSH') then
       data(k)%var = var !!!NBNBNB + std**2 
     else 
       data(k)%var = var
     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     data(k)%depth = 0.0
  enddo ! k = 1, gridpoints(gr)
999 continue
  nrdat =k-1
  print*, 'Number of data read:', nrdat
  close(10)

end subroutine read_CLS_data

end module m_read_CLS_data
