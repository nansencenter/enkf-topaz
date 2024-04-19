!KAL -- A change of Laurents cosistency program to fit with the new stuff
!KAL -- it is called on a file-by-file (member - by - member) basis
!KAL --
!KAL -- Two argument, the file to check, and the corresponding ensemble member. 
!KAL -- the latter is needed when checking ice files as well.
!KAL -- Output file is  consistency_"input_file"
!KAL

program consistency
  use mod_raw_io
  use m_parse_blkdat
  use m_get_mod_grid
  use m_get_mod_fld
  use nfw_mod
  use mod_testinfo
  implicit none

  integer*4, external :: iargc
  integer, parameter :: maxweird = 100 ! give up reporting after 100 weird values 
  real, parameter :: onem=9806.0, undef=-1e14

  logical isweird                      ! error locally at certain i,j,k
  integer iens,i,j,k, count_weird ! counting weird values in A

  character(len=80) :: rstfile, rstbase, outfile, icerstfile
  character(len=8) :: cfld
  character(len=80) :: cmem

  real, dimension(:,:), allocatable :: readfld,modlon,modlat,depths, dpsum
  real, dimension(:,:), allocatable :: countprobs, dpprobs, tempprobs
  real*4, dimension(:,:),  allocatable :: fldr4
  logical :: process_ice
  integer :: iter
  integer :: ncid
  integer :: dimids(2) 
  integer :: lon_id, lat_id, tot_id,tem_id, dp_id

  real    :: bmin, bmax
  real*4  :: amin, amax, spval
  integer :: vlevel,tlevel,rstind
  logical :: readok, ex
  integer :: testindex
  integer :: idm,jdm,kdm
  integer :: ios, fnd
  integer :: imem, itime
  real    :: rdummy

  real*8, allocatable, dimension(:,:) :: ficem,hicem,hsnwm,ticem,tsrfm
  integer :: reclice

  type(testinfo), dimension(:), allocatable :: tests
  integer :: numtest

  character(len=80) :: ncfile
  integer :: dimx,dimy,var_id,var2d(2),ierr

   real :: mindx,meandx

  spval=2**100
   process_ice=.false.
   imem=1 ! Only really needed when reading ice restart file
   if (iargc()==1) then
      call getarg(1,rstfile)
   elseif (iargc()==3) then
      call getarg(1,rstfile)
      call getarg(2,icerstfile)
      call getarg(3,cmem)
      read(cmem,*) imem
      process_ice=.true.
   else
      print *,'Usage: consistency restartfile icerestartfile ensemble_member'
      call exit(1)
   end if
   fnd=max(index(rstfile,'.a'),index(rstfile,'.b'))
   rstbase=rstfile(1:fnd-1)

   ! Set up test info
   numtest=9
   allocate(tests(9))

   call tests_init(tests, numtest)

   !Get model dimensions - 
   call parse_blkdat('idm   ','integer',rdummy,idm)
   call parse_blkdat('jdm   ','integer',rdummy,jdm)
   call parse_blkdat('kdm   ','integer',rdummy,kdm)

   ! Allocate fields
   allocate(readfld(idm,jdm))
   allocate(fldr4  (idm,jdm))
   allocate(modlon (idm,jdm))
   allocate(modlat (idm,jdm))
   allocate(depths (idm,jdm))
   allocate(dpsum  (idm,jdm))
   allocate(countprobs (idm,jdm))
   allocate(dpprobs    (idm,jdm))
   allocate(tempprobs  (idm,jdm))

   ! Get model grid
   call get_mod_grid(modlon,modlat,depths,mindx,meandx,idm,jdm)

   ! Loop through the file header, extract field information, then:
   !   |--> extract header info
   !   |--> match header info with test cases
   !     |--> If match, read field from .a - file
   !     |--> Check for inconcistencies
   outfile='consistency_'//trim(rstbase)
   open(10, file=trim(outfile), access='sequential',status='replace')
   write(10,*) '**************************************************************'
   write(10,*) 'THIS FILE CONTAINS ERRORS DETECTED IN THE HYCOM STATE VARIABLE'
   write(10,*) 'Go to the end of this file for a summary of all errors '
   write(10,*) '**************************************************************'
   write(10,*) ''
   countprobs=0.
   dpprobs=0.
   tempprobs=0.
   rstind=1
   ios=0
   count_weird=0
   readok=.true.
   do while (readok)
      ! Get header info
      call rst_header_from_index(trim(rstbase)//'.b', &
           cfld,vlevel,tlevel,rstind,bmin,bmax,.true.)

      readok=tlevel/=-1 ! test to see if read was ok
      if (readok) then
         call matchtest(cfld,tests,numtest,testindex)

         if (testindex/=-1) then
            print *,'Checking : ',cfld,vlevel,tlevel
            call READRAW(fldr4,amin,amax,idm,jdm,.true.,spval,&
                 trim(rstbase)//'.a',rstind)
            readfld=fldr4

            write(10,'(a,2i5)') 'Testing '//cfld//' at time and layer :',&
                 tlevel,vlevel

            do j=1,jdm
               do i=1,idm
                  if (depths(i,j)>.1) then

                     isweird=.false.
                     if (readfld(i,j)>tests(testindex)%max) then
                        tests(testindex)%toolarge=tests(testindex)%toolarge+1
                        isweird=.true.
                     else if (readfld(i,j)<tests(testindex)%min) then
                        tests(testindex)%toosmall=tests(testindex)%toosmall+1
                        isweird=.true.
                     end if

                     if (tests(testindex)%toosmall + tests(testindex)%toolarge&
                          <maxweird .and. isweird) then
                        write(10,'(a,4i5,e14.4)') '   '//cfld//'&
                             Error at i,j,z,t:',i,j,vlevel,tlevel,readfld(i,j)
                     end if

                     if (isweird) countprobs(i,j)=countprobs(i,j)+1
                     if (isweird.and.trim(cfld)=='dp')&
                          dpprobs(i,j)=dpprobs(i,j)+1
                     if (isweird.and.trim(cfld)=='temp')&
                          tempprobs(i,j)=tempprobs(i,j)+1
                  end if
               end do
            end do

            if ( tests(testindex)%toosmall + tests(testindex)%toolarge&
                 >=maxweird) then
               write(10,*) 'Found ', tests(testindex)%toosmall +&
                    tests(testindex)%toolarge,' errors for ', cfld
            end if
         else
            print *,'Skipping : ',cfld,vlevel,tlevel
         end if
      end if

      rstind=rstind+1
   end do

   ! KAL -- test 2, see if layer thicknesses sum up to depths*onem
   do itime=1,2
      dpsum=0.
      do k=1,kdm
         call get_mod_fld_new(trim(rstbase),readfld,imem,'dp      ',k,itime,&
              idm,jdm,0)
         dpsum=dpsum+readfld
      end do
      print '(a,i3)','Max difference dpsum / depths at time index ',itime
      print *,maxval(dpsum-depths*onem)/onem
   end do

   ! KAL -- test 3, see if ice thickness 
   if (process_ice) then
      allocate(ficem(idm,jdm))
      allocate(hicem(idm,jdm))
      allocate(hsnwm(idm,jdm))
      allocate(ticem(idm,jdm))
      allocate(tsrfm(idm,jdm))

      print *,'TODO -- check ice fields'
      inquire(iolength=reclICE)ficem,hicem,hsnwm,ticem,tsrfm  !,iceU,iceV
      open(20,file=trim(icerstfile),form='unformatted',access='direct',recl=reclICE,status='old')
      read(20,rec=imem,iostat=ios)ficem,hicem,hsnwm,ticem,tsrfm !,iceU,iceV
      close(20)

      do iter=1,2
         if (iter==1) then
            cfld='icec'
            call matchtest(cfld,tests,numtest,testindex)
            readfld=ficem
         else if (iter==2) then
            cfld='hice'
            call matchtest(cfld,tests,numtest,testindex)
            readfld=hicem
         end if

         if (testindex/=-1) then
            write(10,'(a,2i5)') 'Testing '//cfld//' at time and layer :',&
                 tlevel,vlevel

            do j=1,jdm
               do i=1,idm
                  if (depths(i,j)>.1) then
                     isweird=.false.
                     if (readfld(i,j)>tests(testindex)%max) then
                        tests(testindex)%toolarge=tests(testindex)%toolarge+1
                        isweird=.true.
                     else if (readfld(i,j)<tests(testindex)%min) then
                        tests(testindex)%toosmall=tests(testindex)%toosmall+1
                        isweird=.true.
                     end if
                     
                     if (tests(testindex)%toosmall + tests(testindex)%toolarge&
                          <maxweird .and. isweird) then
                        write(10,'(a,4i5,e14.4)') '   '//cfld//&
                             ' Error at i,j,z,t:',i,j,vlevel,tlevel,readfld(i,j)
                     end if
                  end if

                  if (isweird) countprobs(i,j)=countprobs(i,j)+1
               end do
            end do

            if (tests(testindex)%toosmall + tests(testindex)%toolarge&
                 >=maxweird) then
               write(10,*) 'Found ', tests(testindex)%toosmall +  tests(testindex)%toolarge,' errors for ',cfld
            end if

         end if
      end do
   end if
   close(10)

   print *,minval(countprobs),maxval(countprobs)
   where (depths<.1) 
      countprobs=undef
      dpprobs=undef
      tempprobs=undef
   end where
   print *,minval(countprobs),maxval(countprobs)

   ! Netcdf - distribution of "problematic" areas

   ncfile=trim(outfile)//'.nc'
   call nfw_create(ncfile, nf_clobber, ncid)
   call nfw_def_dim(ncfile, ncid, 'idm', idm, dimids(1))
   call nfw_def_dim(ncfile, ncid, 'jdm', jdm, dimids(2))
   call nfw_def_var(ncfile, ncid, 'lon', nf_float, 2, dimids, lon_id)
   call nfw_def_var(ncfile, ncid, 'lat', nf_float, 2, dimids, lat_id)
   call nfw_def_var(ncfile, ncid, 'tempprob', nf_float, 2, dimids, tem_id)
   call nfw_def_var(ncfile, ncid, 'totprob', nf_float, 2, dimids, tot_id)
   call nfw_def_var(ncfile, ncid, 'dpprob', nf_float, 2, dimids, dp_id)
   call nfw_enddef(ncfile, ncid)

   call nfw_put_var_double(ncfile, ncid, lon_id, modlon)
   call nfw_put_var_double(ncfile, ncid, lat_id, modlat)
   call nfw_put_var_double(ncfile, ncid, tot_id, countprobs)
   call nfw_put_var_double(ncfile, ncid, tem_id, tempprobs)
   call nfw_put_var_double(ncfile, ncid, dp_id, dpprobs)
   call nfw_close(ncfile, ncid)



end program consistency
