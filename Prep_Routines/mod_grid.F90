module mod_grid
! Contains the type definition for regular (or irregular grids) together
! with a selection of subprograms for extracting information about the
! grid. 
! 
! 28.1.99, Oyvind.Breivik@nrsc.no.
!
! Future extensions: function checkgrid returns zero if grid contains errors
! or is not set. Include function overloading so that checkgrid may return both
! integer and real.

!!! Module

use mod_angles

!!! Type definition

   ! Type grid contains information for constructing a 1D, 2D, or 3D grid. The 
   ! grid may be periodic and physical units may be added to keep track of
   ! the physical dimensions of start points and resolution of the grid.
   !
   ! Oyvind Breivik, 30.12.98.

   type grid
      integer :: nx, ny, nz ! No of grid points 
      real    :: dx, dy, dz ! Resolution 
      real    :: x0, y0, z0 ! Start point (lower left)
      real    :: undef      ! Undefined value, typically 999.0
      integer :: order      ! 1D, 2D or 3D grid? Default is 2.
      logical :: px, py, pz ! Periodic grid in x, y, z? Default is .false.
      logical :: reg        ! Regular grid? Default is .true.
                            ! If not, order should be 1, indicating an
                            ! array of unevenly spaced data rather than a
                            ! proper grid. In this case, resolution and
                            ! start point become meaningless.
      logical :: set        ! Grid initialized or containing default settings?
      character(len=10) :: ux, uy, uz ! Physical units, 'deg' denotes degrees,
                                      ! default is '1', nondimensional.
   end type grid

   type (grid), parameter :: default_grid = grid(0, 0, 0, 0.0, 0.0, &
      0.0, 0.0, 0.0, 0.0, 999.0, 0, .false., .false., .false., &
      .true., .false., '1', '1', '1')

contains

!!! Subprograms

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function gridpoints(gr)
   ! Calculates the total number of grid points N in a regular grid of
   ! type grid or an irregular array of type grid. Returns zero if grid is not
   ! initialized.

   ! Oyvind Breivik, 30.12.98.

   !!! Interface

      integer gridpoints

      type (grid), intent (in) :: gr

      select case (gr%order)
         case (1)
            gridpoints = gr%nx
         case (2)
            gridpoints = gr%nx*gr%ny
         case (3)
            gridpoints = gr%nx*gr%ny*gr%nz
      end select

      if (.not. gr%reg) then ! Irregular grid?
         gridpoints = gr%nx
      end if

      if (.not. gr%set) then ! Grid initialized or containing default values?
         gridpoints = 0      ! If not initialized, return zero.
      end if

   end function gridpoints


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function gridindex(x,dimid,gr)
   ! Finds corresponding grid index for coordinate x for grid dimension dimid, 
   ! where dimid = 1 denotes x, dimid = 2
   ! denotes y, and dimid = 3 denotes z. If dimid < 0, the grid index is
   ! rounded down using INT, so that the corresponding grid point is ``to
   ! the left'' of x. Otherwise NINT is used and the nearest grid point is
   ! found. 
   !
   ! A return value of zero indicates that x is out of range or grid not 
   ! initialized.
   ! Note that (x-x0) is mapped to [-180, 180] degrees if and only if the
   ! variable ux, uy, or uz (depending again on dimid) equals 'deg'. This is
   ! to ensure that crossing the zero longitude branch cut is handled correctly.
   ! A return value of -1 indicates that dimid is illegal (greater than the
   ! order of the grid).
   !
   ! Requires module mod_angles.

   !!! Interface

      integer gridindex

      real, intent (in)        :: x
      integer, intent (in)     :: dimid
      type (grid), intent (in) :: gr

   !!! Locals

      real    :: x0, x1, dx, e
      integer :: nx
      logical :: closest, deg

   !!! Initialize

      closest = (dimid > 0)

      select case (abs(dimid))  ! Choose correct grid dimension
         case (1) 
            x0 = gr%x0
            dx = gr%dx
            nx = gr%nx
            deg = (gr%ux == 'deg')
         case (2)
            x0 = gr%y0
            dx = gr%dy
            nx = gr%ny
            deg = (gr%uy == 'deg')
         case (3)
            x0 = gr%z0
            dx = gr%dz
            nx = gr%nz
            deg = (gr%uz == 'deg')
      end select

      x1 = x - x0

      if (closest) then
         e = dx/2          ! Small value epsilon
      else
         e = 0.0
      end if

      if (deg) then
         x1 = ang360(x1+e) ! Adding dx/2 is a trick to avoid the branch cut
         x1 = x1-e         ! when finding the closest grid point.
      end if

      if (.not. closest) then
         x1 = x1 - dx/2       ! Round down
      end if

      gridindex = nint(x1/dx) + 1

      if (gridindex < 1 .or. gridindex > nx) then
         gridindex = 0
      end if

      if (abs(dimid) > gr%order) then
         gridindex = -1
      end if

      if (.not. gr%set) then
         gridindex = 0
      end if

   end function gridindex


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function gridpos(i,dimid,gr)
   ! Returns the position of grid node i along dimension dimid in grid.
   !
   ! If dimid < 0 and the physical unit of the grid is degrees, 
   ! -180 <= gridpos < 180 [deg]. Otherwise, 0 <= gridpos < 360 [deg].
   !
   ! Requires module mod_angles.

   !!! Interface 

      real gridpos

      integer,     intent (in) :: i, dimid
      type (grid), intent (in) :: gr

   !!! Locals

      real x0, dx
      logical deg

      select case (abs(dimid))
         case (1)
            x0 = gr%x0
            dx = gr%dx
            deg = (gr%ux == 'deg')
         case (2)
            x0 = gr%y0
            dx = gr%dy
            deg = (gr%uy == 'deg')
         case (3)   
            x0 = gr%z0
            dx = gr%dz
            deg = (gr%uz == 'deg')
      end select

      gridpos = x0 + real(i-1)*dx

      if (deg) then
         if (dimid < 0) then
            gridpos = ang180(gridpos)
         else
            gridpos = ang360(gridpos)
         end if
      end if

   end function gridpos


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function ingrid(x,dimid,gr)
   ! Is x within [x0, x1]? Here x0 and x1 denote the physical 
   ! bounds of the grid along dimension dimid. If dimid < 0 then ingrid
   ! checks the interval [-dx/2+x0, x1+dx/2] instead.
   !
   ! Requires module mod_angles.

   !!! Interface

      logical ingrid

      real,        intent (in) :: x
      integer,     intent (in) :: dimid
      type (grid), intent (in) :: gr

   !!! Locals

      real x0, x1, dx
      integer nx
      logical deg

      select case (abs(dimid))
         case (1)
            dx = gr%dx
            x0 = gr%x0
            nx = gr%nx
            deg = (gr%ux == 'deg')
         case (2)
            dx = gr%dy
            x0 = gr%y0
            nx = gr%ny
            deg = (gr%uy == 'deg')
         case (3)
            dx = gr%dz
            x0 = gr%z0
            nx = gr%nz
            deg = (gr%uz == 'deg')
      end select

      x1 = gridpos(nx,dimid,gr)

      if (dimid < 0) then
         x0 = x0 - dx/2
         x1 = x1 + dx/2
      end if

      ingrid = (x0 <= x) .and. (x <= x1)

      if (deg) then
         ingrid = ang360(x1-x0) >= ang360(x-x0)
      end if

      ingrid = ingrid .and. gr%set
         
   end function ingrid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function undefined(d,gr)
   ! True if d == gr%undef.

      logical undefined

      real, intent (in) :: d
      type (grid), intent (in) :: gr

      undefined = abs(d-gr%undef) < 0.01

   end function undefined


   
end module mod_grid
