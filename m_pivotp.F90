module m_pivotp
  use m_confmap
  implicit none

contains

  ! This subroutine computes the pivot point of each of the observations
  ! in the temporary array tmpobs of type observation. The pivot point
  ! is the biggest i and the biggest j, (i,j) is the computation points/
  ! the grid, that is less than the position of the observation.
  !
  subroutine pivotp(lon, lat, ipiv, jpiv)
   real, intent(in) ::  lon, lat
   integer, intent(out) :: ipiv, jpiv

   real :: tmptan
   real :: lontmp
   
    if (.not. confmap_initialised) then
       print *, 'ERROR: oldtonew(): confmap not initialised'
       stop
    end if

   ! fix for wrap-around
   ! Knut: For some exotic grids the wrap-around
   ! is not needed. By exotic grid I mean Conman,
   ! where the poles are on the other side of the earth,
   ! and the eastern limit is actually WEST of the western 
   ! limit.... (di < 0)
   !if (lon < wlim) then
   if (lon < wlim .and. di > 0. ) then
      lontmp = lon + 360.0
   else
      lontmp = lon
   endif

   ipiv = int((lontmp - wlim) / di) + 1

   if (mercator) then
      if (abs(lat) < 89.999) then
         tmptan = tan(0.5 * rad * lat + 0.25 * pi_1)
         jpiv = int((log(tmptan) - slim * rad) / (rad * dj)) + 1
      else
         jpiv= - 999
      endif
   else
      jpiv = int((lat - slim) / dj) + 1
   endif
 end subroutine pivotp

end module m_pivotp
