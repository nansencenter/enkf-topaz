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
!                16/07/2019 JX:
!                  - Mofication as the reqirements for HYCOM_CICE 
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
   use m_get_mod_fld_nc
   use m_get_mod_grid
   use m_put_mod_fld_nc
   implicit none

   integer*4, external :: iargc
   real, parameter :: onem=9806.
   real, parameter :: PSTARB0 = 1000;
   real                             :: Hflag  !whether assimilate HICE or not

#if defined(HYCOM_CICE)
   integer,                  parameter :: ncat=5
   real,                     parameter :: Thkmin=0.01
   ! 2D variables
   real*8, allocatable, dimension(:,:) :: &
                      ficem,hicem,hsnwm,ticem,tsrfm,sicem,tsnom
   ! 3D variables
   real*8,allocatable,dimension(:,:,:) :: &
                      aicen,vicen,vsnon 

#endif

   integer imem                  ! ensemble member
   integer iyday                 ! date in year 
   character(len=80) :: f_restart,f_icerestart   ! restart nc files 
   character(len=80) :: a_restart,a_icerestart   ! analyzed restart files 

   character(len=80) :: afile,newfile, char80
   integer          :: fnd, rstind, tmpindx, iafile
   logical          :: ex, allok, nomatch
   character(len=8) :: cfld, ctmp
   character(len=3) :: cproc,cmem
   integer          :: tlevel, vlevel, nproc
   real             :: bmin, bmax, rdummy
   integer          :: idm,jdm,kdm
   real, allocatable:: fld(:,:)
   real, allocatable, dimension(:,:) :: depths,modlon,modlat,saln,sss
   real*4, allocatable:: fldr4(:,:)
   real*4 :: spval,amin,amax
   real, allocatable :: press(:)

   real, allocatable, dimension(:,:)   :: dpsum
   real, allocatable, dimension(:,:,:)   :: dp, dpold

   integer,parameter :: numfields=2
   integer :: ios
   integer :: i,j,k

   real :: mindx,meandx

   logical :: Touchice
   integer :: Iceind

   f_restart=''
   f_icerestart=''
   if (iargc()==2 .or. iargc()==6) then
      call getarg(1,a_restart)
      call getarg(2,ctmp)
      read(ctmp,*) imem
      write(cmem,'(i3.3)') imem
      if (iargc()==6)  then
         call getarg(3,f_restart)
         call getarg(4,f_icerestart)
         call getarg(5,ctmp); read(ctmp,*) iyday 
         call getarg(6,ctmp); read(ctmp,*) Hflag 
         print *, Hflag
      end if
   else
      print *,'"fixhycom" -- A crude routine to correct restart files for obvious errors'
      print *
      print *,'usage: '
      print *,'   fixhycom restart_file ensemble_member fore_restart.nc iced_file.nc iyday '
      print *,'   "restart_file"    the restart file you want to fix (.a-file)"'
      print *,'   "ensemble_member" is the ensemble member - should corr. to that of restart file'
      print *,'   "fore_restart"        is optional - it is the restart file for ice fields'
      print *,'   "ice_file"        is optional - it is the restart file for ice fields'
      print *,'   "iyday"           is optional - date in this year'
      call exit(1)
   endif

   ! Get dimensions from blkdat
   call parse_blkdat('idm   ','integer',rdummy,idm)
   call parse_blkdat('jdm   ','integer',rdummy,jdm)
   call parse_blkdat('kdm   ','integer',rdummy,kdm)

   if (idm>0 .and. idm < 1e4 .and. jdm>0 .and. jdm<1e4) then
      allocate(fld  (idm,jdm))
      allocate(fldr4(idm,jdm))
      allocate(sss  (idm,jdm))
      allocate(saln (idm,jdm))
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
   fnd=max(index(a_restart,'.a'),index(a_restart,'.b'))


   ! Inquire for existence
   inquire(exist=ex,file=a_restart(1:fnd-1)//'.b')
   if (.not.ex) then
      write(*,*) 'Can not find '//a_restart(1:fnd-1)//'.b'
      stop '(EnKF_postprocess)'
   end if

   print *,a_restart(1:fnd-1)//'.b'
   newfile='fix'//a_restart(1:fnd-1)


   ! Get model grid
   call get_mod_grid(modlon,modlat,depths,mindx,meandx,idm,jdm)

   !loop over the two time level
   ! Get layer thickness
   dpsum=0.
   do k=1,kdm
      call get_mod_fld_new(a_restart(1:fnd-1),dp(:,:,k),imem,'dp      ',k,1,idm,jdm,0)
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
   !do k = 1, kdm
   !   print *,'max diff is:',maxval(dpold(:,:,k)-dp(:,:,k))/onem,maxloc(dpold(:,:,k)-dp(:,:,k))
   !end do

   ! Loop over restart file
   rstind=1 ! Restart index
   allok=.true.
   do while ( allok)

      ! Get header info from restart
      call rst_header_from_index(a_restart(1:fnd-1)//'.b', &
            cfld,vlevel,tlevel,rstind,bmin,bmax,.true.)

#if defined(HYCOM_CICE)
      Touchice=trim(cfld) /= 'icec' .and. trim(cfld) /= 'hice' .and. &
        trim(cfld) /= 'hsnwm'.and.trim(cfld) /= 'aicen'.and. &
        trim(cfld) /= 'vicen'.and.trim(cfld) /= 'vsnon'.and. &
        trim(cfld) /= 'ticem'.and.trim(cfld) /= 'tsrfm'.and. & 
        trim(cfld) /= 'sicem'.and.trim(cfld) /= 'tsnom'.and. & 
        trim(cfld) /= 'ficem'.and.trim(cfld) /= 'hicem'

      allok=tlevel/=-1.and.Touchice    ! test to see if read was ok
      if (.not.Touchice) Iceind=rstind
#else
      allok=tlevel/=-1 ! test to see if read was ok

#endif

      if (allok ) then 
         print '(a8,i3)',cfld, vlevel
!         Here reading the time record 1 whatever tlevel
         call get_mod_fld_new(a_restart(1:fnd-1),fld(:,:),imem,cfld,vlevel,1,idm,jdm,0)


         if (trim(cfld)=='temp') then

            ! need salinity as well
            ! reading the time record 1 whatever tlevel
            call get_mod_fld_new(a_restart(1:fnd-1),saln(:,:),imem,'saln    ',vlevel,1,idm,jdm,0)
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
           if (vlevel==1) then
              sss=fld
              where(depths<0.1.or.depths>20000)
                 sss=-999.9
              end where
           endif
         else if (trim(cfld)=='pstarb') then
              fld(:,:) = (sqrt(PSTARB0 ** 2 + fld(:,:) ** 2) + fld(:,:)) / 2.0d0;
         else if (trim(cfld)=='dp') then
            !set it equal to the time level 1 that has been corrected
            fld = dp(:,:,vlevel) 
!#if defined (TOPAZ)
!! in west of Mediterranean
!          else if (trim(cfld)=='pbavg' .or. trim(cfld)=='ubavg' .or. trim(cfld)=='vbavg'&
!                 .or. trim(cfld)=='u' .or. trim(cfld)=='v') then
!            do i = 701, 704
!              do j = 482, 484
!                  fld(i, j) = 0.0d0
!              end do
!            end do
!#endif
         end if ! No correction for other fields in the hycom restart file
         !dump the field from the first time level into tlevel (1 or 2)
         call put_mod_fld(trim(newfile),fld,imem,cfld,vlevel,tlevel,rstind,idm,jdm)
         rstind=rstind+1
      end if ! read of template was ok
   end do

   ! put_mod_fld does not include the header from the original file
   ! Supports parallel execution with different members
   open(11,file='tmp'//cmem//'.b',status='replace')
     print *,trim(a_restart(1:fnd-1))//'.b'
     open(20,file=a_restart(1:fnd-1)//'.b',status='old')
       ! Header from original file
       !rewind(11)
       read(20,'(a80)') char80 ; 
       write(11,'(a)',ADVANCE='NO') trim(char80)//achar(13)//achar(10)
       read(20,'(a80)') char80 ; write(11,'(a80)') trim(char80)
     close(20)
  !   open(10,file=trim(newfile)//'.b',status='old')
  !   ! The rest from the newly created file
  !     ios=0
  !     do while(ios==0)
  !       read(10,'(a80)',iostat=ios) char80  
  !       if (ios==0) write(11,'(a80)') trim(char80)
  !     end do
  !   close(10)
   close(11)

   ! Move the new .b file to "newfile"
   !SERIAL call system('mv tmp'//cmem//'.b '//trim(newfile)//'.b')
   ! note to keep mind: tmp???.b to be used replacing the fixrestart or
   ! transport to overwrite the forecasted restart file after postprocessing
   ! for CICE


!###################################################################
!####################### FIX   ICE   MODEL #########################
!###################################################################
#if defined(HYCOM_CICE)
    allocate(ficem(idm,jdm))
    allocate(hicem(idm,jdm))
    allocate(hsnwm(idm,jdm))
    allocate(ticem(idm,jdm))
    allocate(tsrfm(idm,jdm))
    allocate(sicem(idm,jdm))
    allocate(tsnom(idm,jdm))
    allocate(aicen(idm,jdm,ncat))
    allocate(vicen(idm,jdm,ncat))
    allocate(vsnon(idm,jdm,ncat))

   ficem=0; hicem=0; hsnwm=0;
   ticem=0; tsrfm=0;
   sicem=0; tsnom=0;
   aicen=0; vicen=0; vsnon=0;

   ! Loop over restart file
   print *, 'dealing with the ice variables ... '
   rstind=max(1,Iceind) ! Restart index
   allok=.true.

   do while ( allok)

      print '(a8,i3,i3)',cfld, vlevel,tlevel
      ! Get header info from restart
      call rst_header_from_index(a_restart(1:fnd-1)//'.b', &
            cfld,vlevel,tlevel,rstind,bmin,bmax,.true.)

      allok=tlevel/=-1 ! test to see if read was ok

      if (rstind < 0) then
         print *, 'ERROR: fixhycom to read ', trim(a_restart(1:fnd-1)), '.b: "',&
              trim(cfld), '" not found for ',rstind
         stop
      end if
   
      if (allok) then
         ! KAL -- les datafelt vi fann fraa header fila (indx)
         spval=2**100
         call READRAW(fldr4              ,& ! Midlertidig felt som skal lesast
                   amin, amax            ,& ! max/min fraa data (.a) fila 
                   idm,jdm               ,& ! dimensjonar
                   .true.,spval         ,& ! dette brukast for  sette "no value" verdiar
                   trim(a_restart)         ,& ! fil som skal lesast fraa
                   rstind)                  ! index funne over
         ! Sjekk p at vi har lest rett - samanlign max/min fr filene
         if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. &
                 abs(bmax-amax).gt.abs(bmax)*1.e-4     ) then
             print *,'Inconsistency between .a and .b files'
             print *,'.a : ',amin,amax
             print *,'.b : ',bmin,bmax
             print *,cfld,vlevel,tlevel,rstind
             call exit(1)
         end if
        
         if (vlevel>0) then
            select case (cfld)
               case ('aicen')
                 !print *,'aicen : ',amin,amax
                 aicen(:,:,vlevel)=real(fldr4,8);
               case ('vicen')
                 !print *,'vicen : ',amin,amax
                 vicen(:,:,vlevel)=real(fldr4,8);
               case ('vsnon')
                 !print *,'vsnon : ',amin,amax
                 vsnon(:,:,vlevel)=real(fldr4,8);
            end select
         else
            select case (cfld)
               case ('ficem')
                 !print *,'ficem : ',amin,amax
                 ficem=real(fldr4,8);
               case ('hicem')
                 !print *,'hicem : ',amin,amax
                 hicem=real(fldr4,8);
               case ('hsnwm')
                 !print *,'hsnwm : ',amin,amax
                 hsnwm=real(fldr4,8);
            end select
         endif
      endif
      rstind=rstind+1
   end do  ! read the template cycle is end

!  preliminary filter grid mask for sea ice  
   ficem=min(max(0.,ficem), 1.)
   hicem=min(max(0.,hicem),15.)
   hsnwm=min(max(0.,hsnwm), .8)

   where(hicem<=0)
      hicem=0;  hsnwm=0; ficem=0;  
   endwhere
   where(ficem<=0)
      hicem=0;  hsnwm=0
   endwhere

   if (iargc()==6)  then
!      using the same ice restat file as before
!      call fix_cice(ficem,hicem,sss,idm,jdm,ncat, &
!           f_restart,f_icerestart,iyday,0) 

!     Only takeoff the grid if the sea ice melted fully 
      call fix_cice(ficem,hicem,sss,idm,jdm,ncat, &
           f_restart,f_icerestart,iyday,Hflag,1) 

   end if
#endif

print *, 'fixhycom is done!'



end program
