










! File:          p_fixhycom.F90
!
! Created:       ???
!
! Last modified: 29/06/2010
!
! Purpose:       Fixes EnKF output.
!
! Description:   
!
! Modifications:
!                25/10/2011 FC:
!                  - set the two time levels equal
!                29/06/2010 PS:
!                  - set the maximum ICEC to 0.995 to match the model
!                ?/?/? KAL:
!                  - Modification of the "fixhycom" subroutine, into separate
!                    program, working on a file-by-file basis
!                Prior history:
!                  Not documented.

! KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
! KAL -- Input arguments:
! KAL --     template restart file
! KAL --     ensemble member
! KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

program fixhycom




   use mod_raw_io
   use m_parse_blkdat
   use m_put_mod_fld
   use m_get_mod_fld
   use m_get_mod_grid
   implicit none

   integer*4, external :: iargc
   real, parameter :: onem=9806.
   real, parameter :: PSTARB0 = 1000;

   integer imem                  ! ensemble member
   character(len=80) :: restart,icerestart ! restart template

   character(len=80) :: afile,newfile, char80
   integer          :: fnd, rstind, tmpindx, iafile
   logical          :: ex, allok, nomatch
   character(len=8) :: cfld, ctmp
   character(len=3) :: cproc,cmem
   integer          :: tlevel, vlevel, nproc
   real             :: bmin, bmax, rdummy
   integer          :: idm,jdm,kdm
   real, allocatable:: fld(:,:)
   real*8, allocatable, dimension(:,:) :: &
      ficem,hicem,hsnwm,ticem,tsrfm
   real, allocatable, dimension(:,:) :: depths,modlon,modlat,saln
   real*4, allocatable:: fldr4(:,:)
   real*4 :: spval,amin,amax
   real, allocatable :: press(:)

   real, allocatable, dimension(:,:)   :: dpsum
   real, allocatable, dimension(:,:,:)   :: dp, dpold

   integer,parameter :: numfields=2
   integer :: ios,ios2, reclICE,ifld
   integer :: i,j,k

   real :: mindx,meandx


   icerestart=''
   if (iargc()==2 .or. iargc()==3) then
      call getarg(1,restart)
      call getarg(2,ctmp)
      read(ctmp,*) imem
      write(cmem,'(i3.3)') imem
      if (iargc()==3) call getarg(3,icerestart)
   else
      print *,'"fixhycom" -- A crude routine to correct restart files for obvious errors'
      print *
      print *,'usage: '
      print *,'   fixhycom restart_file ensemble_member <ice_file>'
      print *,'   "restart_file"    the restart file you want to fix (.a-file)"'
      print *,'   "ensemble_member" is the ensemble member - should corr. to that of restart file'
      print *,'   "ice_file"        is optional - it is the restart file for ice fields'
      call exit(1)
   endif

   ! Get dimensions from blkdat
   call parse_blkdat('idm   ','integer',rdummy,idm)
   call parse_blkdat('jdm   ','integer',rdummy,jdm)
   call parse_blkdat('kdm   ','integer',rdummy,kdm)

   if (idm>0 .and. idm < 1e4 .and. jdm>0 .and. jdm<1e4) then
      allocate(fld  (idm,jdm))
      allocate(fldr4(idm,jdm))
      allocate(saln (idm,jdm))
      allocate(ficem(idm,jdm))
      allocate(hicem(idm,jdm))
      allocate(hsnwm(idm,jdm))
      allocate(ticem(idm,jdm))
      allocate(tsrfm(idm,jdm))
      allocate(dpsum(idm,jdm))
      allocate(depths(idm,jdm))
      allocate(modlon(idm,jdm))
      allocate(modlat(idm,jdm))
      allocate(dpold(idm,jdm,kdm))
      allocate(dp   (idm,jdm,kdm))
      allocate(press(kdm+1))
   else
      print *,'fld allocate error'
      stop '(EnKF_postprocess)'
   end if


   ! Remove postfix of restart file
   fnd=max(index(restart,'.a'),index(restart,'.b'))


   ! Inquire for existence
   inquire(exist=ex,file=restart(1:fnd-1)//'.b')
   if (.not.ex) then
      write(*,*) 'Can not find '//restart(1:fnd-1)//'.b'
      stop '(EnKF_postprocess)'
   end if

   print *,restart(1:fnd-1)//'.b'
   newfile='fix'//restart(1:fnd-1)


   ! Get model grid
   call get_mod_grid(modlon,modlat,depths,mindx,meandx,idm,jdm)

   !loop over the two time level
   ! Get layer thickness
   dpsum=0.
   do k=1,kdm
      call get_mod_fld_new(restart(1:fnd-1),dp(:,:,k),imem,'dp      ',k,1,idm,jdm)
      dpsum=dpsum+dp(:,:,k)
   end do
   dpold=dp(:,:,:)


   ! DP correction
   do j=1,jdm
   do i=1,idm

      !!! Move negative layers to neighbouring layers.
      do k = 1, kdm-1
         dp(i,j,k+1) = dp(i,j,k+1) + min(0.0,dp(i,j,k))
         dp(i,j,k  ) = max(dp(i,j,k),0.0)
      end do

      !!! Go backwards to fix lowermost layer.
      do k = kdm, 3, -1
         dp(i,j,k-1) = dp(i,j,k-1) + min(0.0,dp(i,j,k))
         dp(i,j,k  ) = max(dp(i,j,k),0.0)
      end do
      !!! No layers below the sea bed.
      press(  1) = 0.0         
      do k = 1, kdm-1
         press(k+1) = press(k) + dp(i,j,k)
         press(k+1) = min(depths(i,j)*onem,press(k+1))
      end do
      press(kdm+1) = depths(i,j)*onem

      do k = 1, kdm
         dp(i,j,k) = press(k+1) - press(k)
      end do
      if (depths(i,j)>100000. .or. depths(i,j) < 1. ) then
        dp(i,j,:)=dpold(i,j,:)
      endif
   end do
   end do
   do k = 1, kdm
      print *,'max diff is:',maxval(dpold(:,:,k)-dp(:,:,k))/onem,maxloc(dpold(:,:,k)-dp(:,:,k))
   end do



   ! Loop over restart file
   rstind=1 ! Restart index
   allok=.true.
   do while ( allok)

      ! Get header info from restart
      call rst_header_from_index(restart(1:fnd-1)//'.b', &
            cfld,vlevel,tlevel,rstind,bmin,bmax,.true.)

      allok=tlevel/=-1 ! test to see if read was ok

      print *,cfld



      if (allok ) then 
!         Here reading the time record 1 whatever tlevel
         call get_mod_fld_new(restart(1:fnd-1),fld(:,:),imem,cfld,vlevel,1,idm,jdm)


         if (trim(cfld)=='temp') then

            ! need salinity as well
            ! reading the time record 1 whatever tlevel
            call get_mod_fld_new(restart(1:fnd-1),saln(:,:),imem,'saln    ',vlevel,1,idm,jdm)
            ! keep water warmer than freezing point
            do j=1,jdm
            do i=1,idm
               fld(i,j)=max(-.057*saln(i,j),fld(i,j))
               fld(i,j)=min(fld(i,j),35.0) !cut off values that are too warm :FC
            end do
            end do
         else if (trim(cfld)=='saln') then
           do j=1,jdm
           do i=1,idm
              fld(i,j)=max(5.,fld(i,j)) ! LB :no water fresher than 5 psu (Baltic)
              fld(i,j)=min(41.,fld(i,j)) ! FC :no water saltier than 40 psu 
           end do
           end do
         else if (trim(cfld)=='pstarb') then
              fld(:,:) = (sqrt(PSTARB0 ** 2 + fld(:,:) ** 2) + fld(:,:)) / 2.0d0;
         else if (trim(cfld)=='dp') then
            !set it equal to the time level 1 that has been corrected
            fld = dp(:,:,vlevel) 
! in west of Mediterranean
          else if (trim(cfld)=='pbavg' .or. trim(cfld)=='ubavg' .or. trim(cfld)=='vbavg'&
                 .or. trim(cfld)=='u' .or. trim(cfld)=='v') then
            do i = 701, 704
              do j = 482, 484
                  fld(i, j) = 0.0d0
              end do
            end do
         else if (trim(cfld)=='ficem') then
           do j=1,jdm
             do i=1,idm
               fld(i,j)=min(max(0.,fld(i,j)),0.995) 
             end do
           end do
         else if (trim(cfld)=='hicem') then
           do j=1,jdm
             do i=1,idm
               fld(i,j)=min(max(0.,fld(i,j)),15.) 
             end do
           end do
         else if (trim(cfld)=='hsnwm') then
           do j=1,jdm
             do i=1,idm
               fld(i,j)=min(max(0.,fld(i,j)),0.4) 
             end do
           end do
         end if ! No correction for other fields in the hycom restart file
         !dump the field from the first time level into tlevel (1 or 2)
         call put_mod_fld(trim(newfile),fld,imem,cfld,vlevel,tlevel,rstind,idm,jdm)
         rstind=rstind+1
      end if ! read of template was ok
   end do


   ! put_mod_fld does not include the header from the original file
   open(10,file=trim(newfile)//'.b',status='old')
   open(20,file=restart(1:fnd-1)//'.b',status='old')
   ! Supports parallel execution with different members
   open(11,file='tmp'//cmem//'.b',status='replace')

   ! Header from original file
   read(20,'(a80)') char80 ; write(11,'(a80)') char80
   read(20,'(a80)') char80 ; write(11,'(a80)') char80
   close(20)

   ! The rest from the newly created file
   ios=0
   do while(ios==0)
      read(10,'(a80)',iostat=ios) char80  
      if (ios==0) write(11,'(a80)') char80
   end do

   close(10)
   close(11)
   close(20)

   ! Move the new .b file to "newfile"
   !SERIAL call system('mv tmp'//cmem//'.b '//trim(newfile)//'.b')



!###################################################################
!####################### FIX      MODEL #########################
!###################################################################
end program
