module mod_measurement_oldnew

   
   type measurement_old
      real d                       ! Measurement value
      real var                     ! Error variance of measurement
      character(len=5) id          ! Type of measurement ('SST', 'SLA', 'COL', 'HICE', 'FICE', 'SAL', 'TEM'
      real lon                     ! Longitude position
      real lat                     ! Latitude position
      real depths                  ! depths of position 
      integer ipiv                 ! i-pivot point in grid
      integer jpiv                 ! j-pivot point in grid
      integer ns                   ! representativity in mod cells (meas. support)
                                   ! ns=0 means: point measurements
                                   ! used in m_Generate_element_Sij.F90
      real a1                      ! bilinear coeffisients (if ni=0)
      real a2                      ! bilinear coeffisients
      real a3                      ! bilinear coeffisients
      real a4                      ! bilinear coeffisients
      logical status               ! active or not
   end type measurement_old

   
   type measurement_new
      real d                       ! Measurement value
      real var                     ! Error variance of measurement
      character(len=5) id          ! Type of measurement ('SST', 'SLA', 'COL', 'HICE', 'FICE', 'SAL', 'TEM'
      real lon                     ! Longitude position
      real lat                     ! Latitude position
      real depths                  ! depths of position 
      integer ipiv                 ! i-pivot point in grid
      integer jpiv                 ! j-pivot point in grid
      integer ns                   ! representativity in mod cells (meas. support)
                                   ! ns=0 means: point measurements
                                   ! used in m_Generate_element_Sij.F90
      real a1                      ! bilinear coeffisients (if ni=0)
      real a2                      ! bilinear coeffisients
      real a3                      ! bilinear coeffisients
      real a4                      ! bilinear coeffisients
      logical status               ! active or not
      integer i_orig_grid          ! KAL - ice drift needs orig grid index as well   !!! NEW
      integer j_orig_grid          ! KAL - ice drift needs orig grid index as well   !!! NEW
   end type measurement_new

contains

   subroutine oldtonew(old,new)
   implicit none

      type (measurement_old), intent(in)  :: old
      type (measurement_new), intent(out) :: new

      new%d     =old%d
      new%var   =old%var
      new%id    =old%id
      new%lon   =old%lon
      new%lat   =old%lat
      new%depths=old%depths
      new%ipiv  =old%ipiv
      new%jpiv  =old%jpiv
      new%ns    =old%ns
      new%a1    =old%a1
      new%a2    =old%a2
      new%a3    =old%a3
      new%a4    =old%a4
      new%status=old%status

      new%i_orig_grid=0
      new%j_orig_grid=0

   end subroutine





end module mod_measurement_oldnew

