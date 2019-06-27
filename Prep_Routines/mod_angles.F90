! Subprograms for the conversion of angles. 
module mod_angles

contains
   
   function ang360(ang)
   ! Maps arbitrary angle to [0, 360) degrees.

      real ang360

      real, intent(in) :: ang

      ang360 = mod(ang, 360.0) - (sign(1.0,ang)-1.0)*180.0

   end function ang360

   function ang180(ang)
   ! Maps arbitrary angle to [-180, 180) degrees.
   ! Use this whenever two angles are subtracted.
   ! Requires ang360.

      real ang180

      real, intent(in) :: ang

      ang180 = ang360(ang)
      ang180 = ang180 - 180.0*(sign(1.0,ang180-180.0)+1.0)

   end function ang180

end module mod_angles
