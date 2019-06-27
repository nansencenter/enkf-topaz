module m_get_mod_grid

logical, parameter, private :: VERBOSE = .false.

contains 
subroutine get_mod_grid(modlon,modlat,depths,mindx,meandx,nx,ny)
   use mod_raw_io
#if defined (QMPI)
   use qmpi
#else
   use qmpi_fake
#endif
   implicit none
   integer, intent(in) :: nx,ny
   real, dimension(nx,ny), intent(out) :: modlon,modlat,depths
   real,intent(out)   :: mindx,meandx

   character(len=7) tag7    
   character(len=80) fname    
   logical ex

   real*4 :: amin,amax,scpx(nx,ny), scpy(nx,ny), spval=0
   real*4 :: infld(nx,ny) ! For reading from .a files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read position from model files
   inquire(file='newpos.uf',exist=ex)
   if (.not.ex) then
      if (master) then
         print *,'newpos.uf file does not exist'
         print *,'(get_mod_grid)'
      end if
      !call stop_mpi()
   else
      if (master .and. VERBOSE) print *,'Reading lon/lat from newpos.uf'
      open(10,file='newpos.uf',form='unformatted',status='old')
        read(10)modlat,modlon
      close(10)
   end if


! Try again -- from regional.grid.a
   if (.not.ex) then
      inquire(exist=ex,file='regional.grid.a')
      if (.not. ex) then
         if (master)  then
            print *,'ERROR: regional.grid.a is not present '
            print *,'(get_mod_grid)'
         end if
         call stop_mpi() ! No more options
      else
         if (master .and. VERBOSE) print *,'Reading lon/lat from regional.grid.a'
         call READRAW(infld,amin,amax,nx,ny,.false.,spval,'regional.grid.a',1)
         modlon=infld
         call READRAW(infld,amin,amax,nx,ny,.false.,spval,'regional.grid.a',2)
         modlat=infld
      end if
   end if



          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read depths from model files
   write (tag7,'(i3.3,a,i3.3)') nx,'x',ny
   fname='depths'//trim(tag7)//'.uf'
   inquire(file=fname,exist=ex)
   if (.not.ex) then
      if (master) then
         print *,'ERROR: depths???x???.uf file does not exist'
         print *,'(get_mod_grid)'
      end if
      call stop_mpi()
   end if
   open(10,file="depths"//trim(tag7)//".uf",form='unformatted',status='old')
      read(10)depths
   close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KAL --  Read scpx, scpy from regional.grid.a

   inquire(exist=ex,file='regional.grid.a')
   if (.not. ex) then
      if (master)  then
         print *,'regional.grid.a is not present '
         print *,'(get_mod_grid)'
      end if
      call stop_mpi()
   end if

   call READRAW(scpx,amin,amax,nx,ny,.false.,spval,'regional.grid.a',10)
   call READRAW(scpy,amin,amax,nx,ny,.false.,spval,'regional.grid.a',11)

   mindx = min(real(minval(scpx)), real(minval(scpy)))
   if (master .and. VERBOSE) then
      print *,'MINIMUM grid size from scpx/scpy : ',mindx
   end if

   meandx = sum(scpx, mask = depths > 1.0d0 .and. depths < 1.0d25) / &
          real(count(depths > 1.0d0 .and. depths < 1.0d25))
   if (master .and. VERBOSE) then
      print *,'MEAN grid size from scpx/scpy : ',meandx
   end if

   ! Safety check ..
   if (mindx<2000.) then
      if (master) then
         print *,'min grid size lower than safety threshold - fix if you want'
         print *,'(get_mod_grid)'
      end if
      call stop_mpi()
   end if

   ! Safety check .. This one is not that critical so the value is set high
   if (mindx>500000.) then
      if (master) then
         print *,'min grid size higher than safety threshold - fix if you want'
         print *,'(get_mod_grid)'
      end if
      call stop_mpi()
   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine  get_mod_grid
end module  m_get_mod_grid
