 module m_read_CLS_header
! Reads the CLS header stored as sla.hdr or sst.hdr

contains
  subroutine read_CLS_header(fnamehdr,gr,dataformat,form,factor,var)
  use mod_grid
  implicit none

  type (grid), intent(out)       :: gr

  character(len=80),intent(in) :: fnamehdr
  character(len=80),intent(out) :: dataformat
  character(len=80) :: title
  character(len=3), intent (out)::form
  integer :: lastchar
  real, intent(out) :: factor, var

  logical :: ex

  gr = default_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read .hdr file information
   inquire (file=fnamehdr,exist=ex)
   if (ex) then
     open (10,file=fnamehdr)
     read (10,'(a)') title
     read (10,'(a3)') form
     read (10,'(a80)') dataformat
     read (10,*) gr%undef, factor, var
     lastchar = scan(dataformat,')')
     dataformat(lastchar+1:80) = ' '
     print '(2a)','title      : ', trim(title)
     print '(2a)','Form       : ', form
     print '(2a)','data format: ', trim(dataformat)
     print '(a,3e14.3)','undef factor var: ', gr%undef,factor,var

!Reads the observation gridd
      if (form == 'reg') then
         read (10,*) gr%nx,gr%ny,gr%x0,gr%y0,gr%dx,gr%dy
         gr%reg = .true.
         gr%order = 2
         gr%ux = 'deg'
         gr%uy = 'deg'
      elseif (form == 'irr') then
         read (10,*) gr%nx
         gr%reg = .false.
         gr%order = 1
      else
         stop 'readhdr: Header error, format should be reg or irr'
      end if
   gr%set = .true.
   close (10)

   else
      form = '0' ! File not found.
      gr%set = .false.
      print*, title
      stop 'read_CLS_header: no data header'
   end if

      print *,' No of gridpoints: ', gridpoints(gr)
!      print '(a,a3,a,f8.4)','Error variance of dataset ',obstype, ' is', var
end subroutine read_CLS_header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_read_CLS_header
