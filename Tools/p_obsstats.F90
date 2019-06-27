program p_obsstats
! Computes the EnKF analysis

! For parallelization
#if defined (QMPI)
   use qmpi
   use distribute
#else
   use qmpi_fake
#endif
   use mod_measurement
   use mod_sphere_tools
   use m_get_mod_grid
   use m_get_mod_nrens
   use m_uobs
   use m_obs
   use m_prep_4_EnKF
   use m_set_random_seed2
   use m_parse_blkdat
   use netcdf
   implicit none
   
   integer nrens                               ! Size of ensemble
   real, dimension(:,:), allocatable :: modlon,modlat,depths

!Generated variables for EnKF analysis
!Global Analysis
   real, allocatable, dimension(:,:) :: S
   real, allocatable, dimension(:,:) :: E
   real, allocatable, dimension(:)   :: d, meanD, meanS, RMSD, RMSE, RMSS, &
      mask_obs

!Local Analysis
   real radius
   integer n_obs_local
   character(len=12) clocal

!Local variables in main program
   integer iargc
   integer i,j

!Parallab: Given random seed
   integer seedsze
   integer, allocatable ::  putseed(:)

!KAL   real, allocatable, dimension(:,:,:,:) :: subS, X3
   real(8) rtc, old, time0, time1

!KAL -- just for testing parse_blkdat
   real :: rdummy
   integer :: idm, jdm, kdm

   real :: mindx,meandx


#if defined(MATLAB)
!#include </export/fimm/local/Matlab-R14sp3/extern/include/fintrf.h> fimm
#include </usr/local/Matlab-6.5/extern/include//fintrf.h> 
   MWPOINTER :: mxCreateNumericMatrix, mxGetPr, mxClassIDFromClassName, matopen,  &
      mxCreateDoubleMatrix, matPutVariableAsGlobal, mp, pa1
   integer matputvariable, matclose
   real*8, dimension(:,:), allocatable :: matio
   integer :: status
#endif

   ! Netcdf output
   integer :: obsdim, var_id,ncid,ierr2
   character(len=80) :: ncfile

   logical :: ex
   character(len=20) :: regname
   integer,parameter :: maxcorners=10
   real, dimension(maxcorners) :: crnlon, crnlat
   integer           :: numcorners,nrobsbox
   integer, dimension(:), allocatable :: pointinbox


   integer :: iuobs

   integer :: ios

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read the characteristics of the assimilation to be carried out.

   if (iargc()/=0) then
      stop 'usage: obstats [ no argument ] '
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set a variable random seed
!//lab   call set_random_seed2

! Remove any randomness in the results for Parallab
   call random_seed(size=seedsze)
   allocate(putseed(seedsze))
   putseed(:)=13
   call random_seed(put=putseed)
   deallocate(putseed)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Get model dimensions - first step on the path to a model-independent EnKF
   call parse_blkdat('idm   ','integer',rdummy,idm)
   call parse_blkdat('jdm   ','integer',rdummy,jdm)
   call parse_blkdat('kdm   ','integer',rdummy,kdm)

! Allocate model grid
   allocate(modlon(idm,jdm))
   allocate(modlat(idm,jdm))
   allocate(depths(idm,jdm))
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read measurements and store in d
   if (master) then
      print*,' '
      print*,'Start reading files, calling obs_readobs()'
   end if
   call obs_readobs
   if (master) then
      print '(2a,i6)', obs(1)%id,' Number of obs =',nobs
      print '(2a,2f6.2)', obs(1)%id,'first obs and var= ',obs(1)%d, obs(1)%var
      print '(2a,2f6.2)', obs(nobs)%id,'last  obs and var= ',obs(nobs)%d, obs(nobs)%var
   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Get model grid
   if (master) then
      print*,' '
      print*,'EnKF: Start reading files, get_mod_grid'
   end if
   call get_mod_grid(modlon,modlat,depths,mindx,meandx,idm,jdm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read ensemble and store in A

   nrens=get_mod_nrens(idm,jdm)
   if (master) print *,'NR ENS (available) = ',nrens,' !!!!'
   !call stop_mpi()
   if (master) then
      print*,'EnKF: Start calculations of input to the analysis'
   end if

   allocate(S(nobs, nrens), d(nobs))
   allocate(meanD(nobs), meanS(nobs))

time0 = rtc()

! KAL -- move ssh stuff in here, as well as reading "observed" variable
   call uobs_get(obs%id, nobs, .true.)
   call obs_prepareobs
   print*, 'Unique obs ', unique_obs(:)
   call prep_4_EnKF(nrens,d,S,depths,meandx*0.001,idm,jdm,kdm)

   print *,'prep_4_EnKF ferdig'



! mean innovations 

   ! old code (allegedly wrong - PS)
   !
   !meanD(:)=0.
   !do j=1,nrens
   !   do i=1,nobs
   !      meanD(i)=meanD(i)+D(i,j)
   !   enddo
   !enddo
   !meanD(:)=meanD(:)/real(nrens)

   meanD(1 : nobs) = obs(1 : nobs) % d - meanS(1 : nobs)

! Innovation RMS without obs. perturbations
   allocate(RMSD(nobs))
   RMSD(:)=0.
   do i=1,nobs
      RMSD(i)=RMSD(i)+meanD(i)**2
   enddo
   RMSD(:)=sqrt(RMSD(:))

! observation std. deviation
   allocate(RMSE(nobs))
   RMSE(:)=0.
   do j=1,nrens
      do i=1,nobs
         RMSE(i)=RMSE(i)+E(i,j)**2
      enddo
   enddo
   RMSE(:)=sqrt(RMSE(:)/real(nrens-1))

! model std. deviation
   allocate(RMSS(nobs))
   RMSS(:)=0.
   do j=1,nrens
      do i=1,nobs
         RMSS(i)=RMSS(i)+S(i,j)**2
      enddo
   enddo
   RMSS(:)=sqrt(RMSS(:)/real(nrens-1))

   print *,'obs stats ferdig'

! In case there are more than 1 obs -- order them in this mask
   allocate(mask_obs(nobs))
!   nuobs=1
!   unique_obs(nuobs)=obs(1)%id
!   do i=1,nobs
!      if (all(unique_obs(1:nuobs)/=obs(i)%id)) then
!         nuobs=nuobs+1
!
!         if (nuobs>maxnuobs) then
!            print *, '(obsstats: too many unique obs ids)'
!            call exit(1)
!         end if
!
!         unique_obs(nuobs)=obs(i)%id
!      end if
!   end do
      
   do iuobs=1,nuobs
      where ( obs%id == unique_obs(iuobs) )  mask_obs=iuobs
   end do

   print *,'obs mask done'
      

   ! Produce tec fields - should be split into several files/zones if observations
   ! are of different types
   open(10,file='innovationstats.tec',status='replace')
   write(10,*)'TITLE = "innovation statistics"'
   write(10,*)'VARIABLES = "i-index" "j-index" "lon" "lat" "meaninnov"' // &
              '"RMSinnov" "stdobs" "stdmodel" "maskobs" '
   write(10,'(a,i7)')' ZONE  F=BLOCK, I=',nobs
   write(10,900)(obs(i)%ipiv + obs(i)%a1 + obs(i)%a4,i=1,nobs)
   write(10,900)(obs(i)%jpiv + obs(i)%a2 + obs(i)%a3,i=1,nobs)
   write(10,900)(obs(i)%lon                         ,i=1,nobs)
   write(10,900)(obs(i)%lat                         ,i=1,nobs)
   write(10,900)(meanD(i)                           ,i=1,nobs)
   write(10,900)(RMSD (i)                           ,i=1,nobs)
   write(10,900)(RMSE (i)                           ,i=1,nobs)
   write(10,900)(RMSS (i)                           ,i=1,nobs)
   write(10,900)(mask_obs(i)                        ,i=1,nobs)
   close(10)
 900 format(10(1x,e12.5))

#if defined (MATLAB)
   ! Do the same for matlab -- only partly finished
   print *,'matlab'
   allocate(matio (nobs,1))
   mp=matopen('innovationstats.mat','w')

   matio(:,1)=obs(:)%lon
   pa1=mxCreateNumericMatrix(nobs,1,mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(matio,mxGetPr(pa1),nobs)
   status = matPutVariable(mp, 'lon', pa1)

   matio(:,1)=obs(:)%lat
   pa1=mxCreateNumericMatrix(nobs,1,mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(matio,mxGetPr(pa1),nobs)
   status = matPutVariable(mp, 'lat', pa1)

   matio(:,1)=meanD
   pa1=mxCreateNumericMatrix(nobs,1,mxClassIDFromClassName('double'),0)
   call mxCopyReal8ToPtr(matio,mxGetPr(pa1),nobs)
   status = matPutVariable(mp, 'mean_innov', pa1)
#endif

   ! Netcdf - safest bet
   ncfile='innovationstats.nc'
   if (NF90_CREATE(trim(ncfile),NF90_CLOBBER,ncid) /= NF90_NOERR) then
      print *,'An error occured when opening the netcdf file'
      stop '(obsstats)'
   end if
   ierr2=NF90_DEF_DIM(ncid,'nobs',nobs,obsdim)

   ierr2=NF90_DEF_VAR(ncid,'lon',NF90_Float,obsdim,var_id)
   ierr2=NF90_ENDDEF(ncid)
   ierr2=NF90_PUT_VAR(ncid,var_id,obs(:)%lon)

   ierr2=NF90_REDEF(ncid)
   ierr2=NF90_DEF_VAR(ncid,'lat',NF90_Float,obsdim,var_id)
   ierr2=NF90_ENDDEF(ncid)
   ierr2=NF90_PUT_VAR(ncid,var_id,obs(:)%lat)

   ierr2=NF90_REDEF(ncid)
   ierr2=NF90_DEF_VAR(ncid,'meaninnov',NF90_Float,obsdim,var_id)
   ierr2=NF90_ENDDEF(ncid)
   ierr2=NF90_PUT_VAR(ncid,var_id,meanD)

   ierr2=NF90_REDEF(ncid)
   ierr2=NF90_DEF_VAR(ncid,'RMSinnov',NF90_Float,obsdim,var_id)
   ierr2=NF90_ENDDEF(ncid)
   ierr2=NF90_PUT_VAR(ncid,var_id,RMSD)

   ierr2=NF90_REDEF(ncid)
   ierr2=NF90_DEF_VAR(ncid,'varobs',NF90_Float,obsdim,var_id)
   ierr2=NF90_ENDDEF(ncid)
   ierr2=NF90_PUT_VAR(ncid,var_id,RMSE)

   ierr2=NF90_REDEF(ncid)
   ierr2=NF90_DEF_VAR(ncid,'varmodel',NF90_Float,obsdim,var_id)
   ierr2=NF90_ENDDEF(ncid)
   ierr2=NF90_PUT_VAR(ncid,var_id,RMSS)

   ierr2=NF90_REDEF(ncid)
   ierr2=NF90_DEF_VAR(ncid,'maskobs',NF90_Float,obsdim,var_id)
   ierr2=NF90_ENDDEF(ncid)
   ierr2=NF90_PUT_VAR(ncid,var_id,mask_obs)

   ierr2=NF90_CLOSE(ncid)



! Global stats

   print *,'Global mean innovation      : ',sum(meanD(1 : nobs))/nobs
   print *,'Global mean innovation RMS  : ',sum(RMSD(1 : nobs))/nobs
   print *,'Global mean obs    variance : ',sum(RMSE(1 : nobs))/nobs
   print *,'Global mean model  variance : ',sum(RMSS(1 : nobs))/nobs


   ! Regional stats - only if the file below exists
   inquire(exist=ex,file='EnKF_regstats.in',iostat=ios)
   if (ios.ne.0) stop 'obsstat: error opening EnKF_regstats.in'
   if (ex) then
      allocate(pointinbox(nobs))
      print *,'Found EnKF_regstats.in - producing regional statistics'
      open(10,file='EnKF_regstats.in',status='old') 

      do while (ios==0) 

         read(10,'(a)',iostat=ios) regname
         if (ios/=0) cycle
         read(10,*,iostat=ios) numcorners
         print *, 'Region is ', regname, ' and has ', numcorners, ' corners'

         if (numcorners>maxcorners) then
            print *,'obsstats can only handle ',maxcorners,' corners'
            call exit(1)
         end if

         read(10,*,iostat=ios) crnlon(1:numcorners)
         read(10,*,iostat=ios) crnlat(1:numcorners)

         ! Find points in box
         pointinbox=0
         do i=1,nobs
            if ( inbox(crnlon,crnlat,numcorners,obs(i)%lon,obs(i)%lat)) &
            pointinbox(i)=1
         end do

         print *,'Total nobs     :',nobs
         print *,'Total in testbox:',sum(pointinbox)

         nrobsbox=sum(pointinbox)
         if (nrobsbox==0) cycle

         print *,trim(regname)//' mean innovation      : ',sum(meanD(1:nobs)*pointinbox(1:nobs))/nrobsbox
         print *,trim(regname)//' mean innovation RMS  : ',sum(RMSD*pointinbox )/nrobsbox
         print *,trim(regname)//' mean obs    variance : ',sum(RMSE*pointinbox )/nrobsbox
         print *,trim(regname)//' mean model  variance : ',sum(RMSS*pointinbox )/nrobsbox

      end do

      close (10)

   else
      print *,'EnKF_regstats.in not found - skipping regional statistics'
   end if


  print*,'obsstats: Finished'




end program 

