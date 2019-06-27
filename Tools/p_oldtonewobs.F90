program oldtonewobs
   use mod_measurement_oldnew
   implicit none

   type (measurement_old) :: oldobs
   type (measurement_new) :: newobs

   integer :: reclold,reclnew
   integer :: iobs
   integer :: iosin, iosout
   logical :: ex


   ! copy old obs 
   inquire(exist=ex,file='observations.uf')
   if (ex) then
      call system("cp observations.uf old_observations.uf")
   else
      print *,'observations.uf does not exist'
      call exit(1)
   end if

   ! Open old and new obs files
   inquire(iolength=reclold) oldobs
   inquire(iolength=reclnew) newobs
   open(10, file='old_observations.uf',status='old',recl=reclold,access='direct')
   open(11, file='new_observations.uf',status='replace',recl=reclnew,access='direct')

   iosin=0
   iosout=0
   iobs=1
   do while (iosin==0 .and. iosout==0)

      read(10,rec=iobs,iostat=iosin) oldobs

      if (iosin==0) then

         call oldtonew(oldobs,newobs)

         write(11,rec=iobs,iostat=iosout) newobs

         !print *,newobs%ipiv,newobs%jpiv

         iobs=iobs+1

         if (iosout/=0) then
            print *,'Error when writing to new obs'
            print *,'(oldtonewobs)'
            call exit(1)
         end if
      end if

   end do
   close(10)
   close(11)

   print *,'Processed ',iobs-1,' observations'




   







end program  oldtonewobs

