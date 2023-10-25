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
!                1/03/2011 Ehouarn:
!                  -modification of the fixhycom subroutine: interpolation of 
!                    biogeochemical tracers on the analysis grid according to
!                    hycom remapping in order to insure conservation.  
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
! KAL -- TODO: Fix for ice fields

program fixhycom_eco




   use mod_raw_io
   use m_parse_blkdat
   use m_put_mod_fld
   use m_get_mod_fld
   use m_get_mod_grid
#if defined (ECO)
   use m_fixhycom_eco_metno
#endif   
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
   integer :: i,j,k,ktr

   real :: mindx,meandx
   
#if defined (ECO)   
    real,dimension(:,:,:), allocatable::dpfor,cfi 
    character(len=80) :: restfor
    real,dimension(:,:,:,:), allocatable::tracerf
    real, dimension(:,:), allocatable::trcraij
    real, dimension(:), allocatable::prsf,dpf
    integer::ntracr,ktrcr
    real::dpthin
    character(2)::ctrcr
    logical, dimension(:), allocatable::lcm
    
    real, dimension(:,:,:), allocatable::temp,sal
    integer::kisop
#endif

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
#if defined (ECO)
   call parse_blkdat('ntracr','integer',rdummy,ntracr) 	  
#endif
      
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
#if defined (ECO)
      allocate(dpfor(idm,jdm,kdm))
      allocate(cfi(kdm,ntracr,2))
      allocate(prsf(kdm+1))
      allocate(tracerf(idm,jdm,kdm,ntracr))
      allocate(trcraij(kdm,ntracr))
      allocate(lcm(kdm))
      allocate(dpf(kdm))
      allocate(sal(idm,jdm,kdm))
      allocate(temp(idm,jdm,kdm))
#endif
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

#if defined (ECO)
   !files where are stored the forecast fields!
   restfor='forecast'//cmem
   dpthin = onem*0.001
#endif

   ! Get layer thickness
   dpsum=0.
   do k=1,kdm
      !call get_mod_fld(dp(:,:,k),1,'dp      ',k,1)
      call get_mod_fld_new(restart(1:fnd-1),dp(:,:,k),imem,'dp      ',k,1,idm,jdm,0)
      dpsum=dpsum+dp(:,:,k)
#if defined (ECO)
      !reading of forecast fields: tracers, T and S
      call get_mod_fld_new(trim(restfor),dpfor(:,:,k),imem,'dp      ',k,1,idm,jdm,0)
      do ktrcr=1,ntracr
         write(ctrcr,'(i2.2)') ktrcr
         cfld='tracer'//ctrcr
         call get_mod_fld_new(trim(restfor),tracerf(:,:,k,ktrcr),imem,cfld,k,1,idm,jdm,0)
      enddo
      call get_mod_fld_new(trim(restfor),temp(:,:,k),imem,'temp    ',k,1,idm,jdm,0)
      call get_mod_fld_new(trim(restfor),sal(:,:,k),imem,'saln     ',k,1,idm,jdm,0)
#endif      
   end do
   print *,maxval(dpsum-depths*onem)
   dpold=dp



   ! DP correction
   do j=1,jdm
   do i=1,idm

      !!! Move negative layers to neighbouring layers.
      do k = 1, kdm-1
         dp(i,j,k+1) = dp(i,j,k+1) + min(0.0,dp(i,j,k))
         dp(i,j,k  ) = max(dp(i,j,k),0.0)
         !dp(i,j,k  ) = max(dp(i,j,k),1.e-3*onem)
      end do

      !!! Go backwards to fix lowermost layer.
      do k = kdm, 3, -1
         dp(i,j,k-1) = dp(i,j,k-1) + min(0.0,dp(i,j,k))
         dp(i,j,k)   =   max(dp(i,j,k),0.0)
	 !dp(i,j,k)   =   max(dp(i,j,k),1.e-3*onem)
      end do

      !!! No layers below the sea bed.
      press(  1) = 0.0
#if defined (ECO)
      !computation of the forecast layer interfaces (prsf)!
      prsf(1)=0.
#endif               
      do k = 1, kdm-1
         press(k+1) = press(k) + dp(i,j,k)
         press(k+1) = min(depths(i,j)*onem,press(k+1))
#if defined (ECO)
         prsf(k+1) = prsf(k) + dpfor(i,j,k)
#endif 	 
      end do
      press(kdm+1) = depths(i,j)*onem
#if defined (ECO)
      prsf(kdm+1)=depths(i,j)*onem
#endif     

      do k = 1, kdm
         dp(i,j,k) = press(k+1) - press(k)
#if defined (ECO)
         !definition of the isopycnal layers!
         dpf(k)=max(dpfor(i,j,k),dpthin) 
         kisop=compute_kisop(temp(i,j,:),sal(i,j,:),kdm)
         if     (k.le.max(2,kisop)) then
          lcm(k) = .false.  !fixed layers are never PCM
         else
! ---       thin and isopycnal layers remapped with PCM.
            lcm(k) = dpfor(i,j,k).le.dpthin
         endif 
#endif  	       
      end do

!eho 2/11/11      
      if(depths(i,j)>10000. .or. depths(i,j)<1.)then
        dp(i,j,:)=dpold(i,j,:)
      endif
      
#if defined (ECO)     
!      if(depths(i,j)==0.)then
      if(depths(i,j)>10000. .or. depths(i,j)<1.)then
        cycle
      else	
	call hybgen_weno_coefs(tracerf(i,j,:,:),dpf,lcm,cfi,kdm,ntracr,dpthin)
        call hybgen_weno_remap(tracerf(i,j,:,:),prsf,dpfor(i,j,:),cfi,trcraij,&
                        press,kdm,kdm,ntracr,dpthin)
        tracerf(i,j,1:kdm,1:ntracr)=trcraij(1:kdm,1:ntracr)  
      endif
      
#endif
      
   end do
   end do

   do k = 1, kdm
      print *,maxval(dpold(:,:,k)-dp(:,:,k))/onem
   end do

#if defined (ECO)     
   deallocate(trcraij,dpf,lcm,cfi,prsf,dpfor)
#endif





   ! Loop over restart file
   rstind=1 ! Restart index
   allok=.true.
   do while ( allok)

      ! Get header info from restart
      call rst_header_from_index(restart(1:fnd-1)//'.b', &
            cfld,vlevel,tlevel,rstind,bmin,bmax,.true.)

      allok=tlevel/=-1 ! test to see if read was ok

      print *,cfld



      if (allok) then 

         ! Get actual field  - for now we use the READRAW routine (later on we
         ! should switch to m_get_mod_fld
!         call READRAW(fldr4,amin,amax,idm,jdm,.false.,spval,restart(1:fnd-1)//'.a',rstind)
!         fld=fldr4

         ! Sjekk p at vi har lest rett - samanlign max/min fr filene
!         if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. &
!                 abs(bmax-amax).gt.abs(bmax)*1.e-4     ) then
!            print *,'Inconsistency between .a and .b files'
!            print *,'.a : ',amin,amax
!            print *,'.b : ',bmin,bmax
!            print *,cfld,vlevel,tlevel
!            call exit(1)
!         end if

	 call get_mod_fld_new(restart(1:fnd-1),fld(:,:),imem,cfld,vlevel,1,idm,jdm,0)

         if (trim(cfld)=='temp') then

            ! need salinity as well
            call get_mod_fld_new(restart(1:fnd-1),saln(:,:),imem,'saln    ',vlevel,1,idm,jdm,0)

            !if (tlevel==-1) then
            !   print *,'Could not get salinity field'
            !   call exit(1)
            !end if

            ! keep water warmer than freezing point
            do j=1,jdm
            do i=1,idm
               fld(i,j)=max(-.057*saln(i,j),fld(i,j))
               fld(i,j)=min(35.,fld(i,j)) !FC: cut off values that are too warm
            end do
            end do
         else if (trim(cfld)=='saln') then
           do j=1,jdm
           do i=1,idm
              fld(i,j)=max(5.,fld(i,j)) ! LB :no water fresher than 5 psu (Baltic)
              fld(i,j)=min(41.,fld(i,j))! FC: no water saltier than 40 psu
           end do
           end do
         else if (trim(cfld)=='pstarb') then
              fld(:,:) = (sqrt(PSTARB0 ** 2 + fld(:,:) ** 2) + fld(:,:)) / 2.0d0;
         else if (trim(cfld)=='dp') then
            fld = dp(:,:,vlevel) ! NB, one time level 
#if defined (ECO)  
         else if (cfld(1:6)=='tracer') then  	    
	    !updating the file!
	    ktrcr=tracr_get_incr(cfld(7:8))
	    if (ktrcr==-1)then
	      print*,'alert tracer unknow'
	      exit
	    endif
            fld(:,:)= tracerf(:,:,vlevel,ktrcr)
	    
#endif	 
	 end if ! No correction for other fields in the hycom restart file

         call put_mod_fld(trim(newfile),fld,imem,cfld,vlevel,tlevel,rstind,idm,jdm)


            
         rstind=rstind+1
      end if ! read of template was ok
   end do

#if defined (ECO)
   deallocate(tracerf)
#endif


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
!####################### FIX   ICE   MODEL #########################
!###################################################################





   if (iargc()==3) then

      inquire(exist=ex,file=trim(icerestart))
      if (.not.ex) then
         print *,icerestart//' does not exist!'
         print *,'(fixhycom)'
         stop
      end if

      inquire(iolength=reclICE)ficem,hicem,hsnwm,ticem,tsrfm  !,iceU,iceV
      open(10,file=icerestart,form='unformatted',access='direct',recl=reclICE)
      read(10,rec=imem)ficem,hicem,hsnwm,ticem,tsrfm !,iceU,iceV
      close(10)

      ! PS 25/06/2010 max(ficem) = 0.995 - max the model allows
      ficem=min(max(0.,ficem), 0.995)
      hicem=min(max(0.,hicem),15.)
      hsnwm=min(max(0.,hsnwm), 4.)

      open (10,file='fix'//trim(icerestart),form='unformatted',access='direct',recl=reclICE)
      write(10,rec=imem)ficem,hicem,hsnwm,ticem,tsrfm !,iceU,iceV
      close(10)
   end if
      




!KAL
!KAL   ! This is to make it easier for supporting programs
!KAL   open(10,file='file.EnKF_postprocess',status='replace')
!KAL   write(10,'(a)') 'analysis'//cmem
!KAL   close(10)
!KAL
!KAL
!KAL
!KAL   ! ice processing Loop over restart file
!KAL   print *
!KAL   print *
!KAL   print *,'processing ice restart file'
!KAL   rstind=1 ! Restart index
!KAL   allok=.true.
!KAL   call system("cp ensemble_TMP_ICE.uf ensemble_TMP_ICE_final.uf")
!KAL   do ifld=1,numfields
!KAL      cfld=  fieldnames (ifld)
!KAL      vlevel=fieldlevels(ifld)
!KAL      if (trim(cfld)=='icec' .or. trim(cfld)=='hice') then
!KAL         nomatch=.true.
!KAL         do iafile=1,nproc ! List of procs used in analysis
!KAL            write(cproc,'(i3.3)') iafile-1
!KAL
!KAL            ! NB - time level=1
!KAL            ! NB2 - the files dumped in the analysis lack a header (last argument
!KAL            ! is false)
!KAL            call rst_index_from_header(trim(afile)//'.b',cfld,vlevel,1, &
!KAL                                             tmpindx,bmin,bmax,.false.) 
!KAL
!KAL
!KAL            if (tmpindx/=-1) then
!KAL               print '(a8," -- layer:",i4,"  match : record, file",i4," ",a)', cfld, vlevel,tmpindx,trim(afile)
!KAL               nomatch=.false.
!KAL               exit
!KAL            end if
!KAL
!KAL            ! Read field from analysed file 
!KAL            call READRAW(fldr4,amin,amax,idm,jdm,.false.,spval,trim(afile)//'.a',tmpindx)
!KAL            fld=fldr4
!KAL
!KAL           ! Sjekk p at vi har lest rett - samanlign max/min fr filene
!KAL           if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. &
!KAL                   abs(bmax-amax).gt.abs(bmax)*1.e-4     ) then
!KAL              print *,'Inconsistency between .a and .b files'
!KAL              print *,'.a : ',amin,amax
!KAL              print *,'.b : ',bmin,bmax
!KAL              print *,cfld,vlevel,tlevel
!KAL              call exit(1)
!KAL           end if
!KAL
!KAL         end do
!KAL
!KAL
!KAL         ! Check if we got ice concentration or ice thickness
!KAL         if (.not.nomatch) then
!KAL
!KAL
!KAL            inquire(iolength=reclICE)ficem,hicem,hsnwm,ticem,tsrfm  !,iceU,iceV
!KAL            open(10,file=trim(icetemplate),form='unformatted',access='direct',recl=reclICE,status='old')
!KAL            open(11,file='ensemble_TMP_ICE_final.uf',form='unformatted',access='direct',recl=reclICE,status='unknown')
!KAL            read(10,rec=imem,iostat=ios)ficem,hicem,hsnwm,ticem,tsrfm !,iceU,iceV
!KAL            if (trim(cfld)=='icec') ficem=fld
!KAL            if (trim(cfld)=='hice') hicem=fld
!KAL            write(11,rec=imem,iostat=ios2)ficem,hicem,hsnwm,ticem,tsrfm !,iceU,iceV
!KAL            close(10)
!KAL            close(11)
!KAL
!KAL            if (ios/=0 .or.ios2/=0) then
!KAL               print *,ios
!KAL               print *,ios2
!KAL               print *,'Error when writing to ice ens file'
!KAL               call exit(1)
!KAL            end if
!KAL         end if
!KAL      end if
!KAL   end do
!KAL
!KAL
!KAL   print *,'Normal exit of EnKF_postprocess'
!KAL   print *,'TODO: Process dp inconsistencies'
!KAL   call exit(0)
!KAL
!KAL
!KAL
end program
