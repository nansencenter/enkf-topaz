










program EnKF_postprocess
! KAL -- The new EnKF is 1-parallelized, and each thread will dump
! KAL -- to its own file named "analysisXXX_procXXX.[ab]
! KAL --
! KAL -- This routine will gather the
! KAL -- analysis from the separate analyzed files into one complete restart file.
! KAL -- To do this, a "template" restart file must be specified. This file
! KAL -- copies non-existing variables in the analyzed fields from the template,
! KAL -- and into the final analysis. The final files produced by this routine 
! KAL -- are named "analysisXXX.[ab]".
! KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
! KAL -- Input arguments:
! KAL --     template restart file
! KAL --     template ice restart file
! KAL --     ensemble member
! KAL --     number of 1 threads used in analysis
! KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
! KAL --
! KAL -- NB: This routine will not check or modify the fields. To do that, 
! KAL --     use "consistency" (check) and "fixhycom" (fix)



   use mod_raw_io
   use m_parse_blkdat
   use m_put_mod_fld
   implicit none

   integer*4, external :: iargc

   integer imem                  ! ensemble member
   character(len=80) :: template,icetemplate ! restart template

   character(len=80) :: afile
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
   real*4, allocatable:: fldr4(:,:)
   real*4 :: spval,amin,amax

   real, allocatable, dimension(:,:)   :: dpsum

   integer,parameter :: numfields=2
   integer :: ios,ios2, reclICE,ifld
   character(len=8) :: fieldnames(numfields)
   integer :: fieldlevels(numfields)


   if (iargc()==4) then
      call getarg(1,template)
      call getarg(2,icetemplate)
      call getarg(3,ctmp)
      read(ctmp,*) imem
      write(cmem,'(i3.3)') imem
      call getarg(4,ctmp)
      read(ctmp,*) nproc
   else
      print *,'usage: EnKF_postprocess restart_template ice_template ensemble_member nproc'
      call exit(1)
   endif

   ! Get dimensions from blkdat
   call parse_blkdat('idm   ','integer',rdummy,idm)
   call parse_blkdat('jdm   ','integer',rdummy,jdm)
   call parse_blkdat('kdm   ','integer',rdummy,kdm)

   if (idm>0 .and. idm < 1e4 .and. jdm>0 .and. jdm<1e4) then
      allocate(fld  (idm,jdm))
      allocate(fldr4(idm,jdm))
      allocate(ficem(idm,jdm))
      allocate(hicem(idm,jdm))
      allocate(hsnwm(idm,jdm))
      allocate(ticem(idm,jdm))
      allocate(tsrfm(idm,jdm))
      allocate(dpsum(idm,jdm))
   else
      print *,'fld allocate error'
      stop '(EnKF_postprocess)'
   end if


   ! Remove postfix of template file
   fnd=max(index(template,'.a'),index(template,'.b'))


   ! Inquire for existence
   inquire(exist=ex,file=template(1:fnd-1)//'.b')
   if (.not.ex) then
      write(*,*) 'Can not find '//template(1:fnd-1)//'.b'
      stop '(EnKF_postprocess)'
   end if

   print *,template(1:fnd-1)//'.b'

   ! Loop over restart file
   dpsum=0.
   rstind=1 ! Restart index
   allok=.true.
   do while ( allok)

      ! Get header info from template
      call rst_header_from_index(template(1:fnd-1)//'.b', &
            cfld,vlevel,tlevel,rstind,bmin,bmax,.true.)

      allok=tlevel/=-1 ! test to see if read was ok


      if (allok) then 

         ! Get actual field  - for now we use the READRAW routine (later on we
         ! should switch to m_get_mod_fld
         call READRAW(fldr4,amin,amax,idm,jdm,.false.,spval,template(1:fnd-1)//'.a',rstind)
         fld=fldr4


         !print *,cfld,tlevel, vlevel


         ! From the template, we have the stuff we need, now we go looking
         ! in the EnKF analysis files output for the correct input. Note that
         ! all time levels are set equal in the restart. This may introduce 
         ! imbalances, but the alternative is probably worse.
         !
         ! KAL -- need list of analysis files here, for now we hardcode
         nomatch=.true.
         do iafile=1,nproc ! List of procs used in analysis

            write(cproc,'(i3.3)') iafile-1
            
            ! Temporary name, will change
            afile='analysis'//cmem//'_proc'//cproc


            ! NB - time level=1
            ! NB2 - the files dumped in the analysis lack a header (last argument
            ! is false)
            if (trim(cfld)=='ficem') then
              ctmp='icec'
            elseif (trim(cfld)=='hicem') then
              ctmp='hice'
            else
              ctmp=trim(cfld)
            endif
            call rst_index_from_header(trim(afile)//'.b',ctmp,vlevel,1, &
                                          tmpindx,bmin,bmax,.false.) 

            if (tmpindx/=-1) then
               !if (tlevel/=1) print *,'--> replacing time level with 1'
               print '(a8," -- time, layer:",2i4,"  match : record, file",i4," ",a)', cfld,tlevel, vlevel,tmpindx,trim(afile)
               nomatch=.false.

               ! Read field from analysed file 
               call READRAW(fldr4,amin,amax,idm,jdm,.false.,spval,trim(afile)//'.a',tmpindx)
               fld=fldr4

              ! Sjekk p at vi har lest rett - samanlign max/min fr filene
              if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. &
                      abs(bmax-amax).gt.abs(bmax)*1.e-4     ) then
                 print *,'Inconsistency between .a and .b files'
                 print *,'.a : ',amin,amax
                 print *,'.b : ',bmin,bmax
                 print *,cfld,vlevel,tlevel
                 call exit(1)
              end if



               ! put into final, processed file -- imem is not used, actually
               call put_mod_fld('analysis'//cmem,fld,imem,cfld,vlevel,tlevel,rstind,idm,jdm)
               

               exit
            end if
         end do


         if (nomatch) then
            print '(a8," -- time, layer:",2i4," - no match - replace with template")',cfld,tlevel,vlevel

            ! put template values into final, processed file -- imem is not used, actually
            call put_mod_fld('analysis'//cmem,fld,imem,cfld,vlevel,tlevel,rstind,idm,jdm)
         end if
            
         rstind=rstind+1
      end if ! read of template was ok
   end do


   print *,'Normal exit of EnKF_postprocess'
   print *,'TODO: Process dp inconsistencies'
   call exit(0)



end program
