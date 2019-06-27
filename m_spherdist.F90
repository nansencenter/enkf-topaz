module m_spherdist

contains

  ! Computes the distance between geo. pos. lon1, lat1 and lon2, lat2.
  ! http://en.wikipedia.org/wiki/Great-circle_distance
  !
  ! Input is in degrees, output in meters
  !
  !
  !FC: 29/02/12 add min max to avoid NaN from acos
real function spherdist(lon1, lat1, lon2, lat2)
  implicit none

  real, intent(in) :: lon1, lat1, lon2, lat2 ! pos. in degrees

  real, parameter :: INVRAD = 3.14159265358979323846d0 / 180.0d0
  real, parameter :: REARTH = 6371000.0d0
  real  :: rlon1, rlat1, rlon2, rlat2 ! pos. in radians

  rlon1 = lon1 * INVRAD !lon1 in rad
  rlat1 = lat1 * INVRAD !90-lat1 in rad 
  rlon2 = lon2 * INVRAD ! lon2 in rad
  rlat2 = lat2 * INVRAD !90 - lat2 in rad 

  spherdist = REARTH * acos(min(max(sin(rlat1) * sin(rlat2)&
         + cos(rlat1) * cos(rlat2) * cos(rlon1 - rlon2),-1.),1.))
end function spherdist

end module m_spherdist
