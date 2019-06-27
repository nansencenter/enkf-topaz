! File :         p_chk_obs.F90
!
!                Created: unknown
!

program p_obs_info
  use mod_measurement
  use nfw_mod
  implicit none


  real,parameter :: Tvar0=0.5
  real,parameter :: Svar0=0.02

  type (measurement), allocatable :: obs0(:)


  integer :: rsize,k,l
  integer :: k0
  logical :: exists = .false.
  type(measurement) :: record
  integer :: ios
  integer :: Iord
 
  type (measurement_info), allocatable :: obsinfo(:)
  type (measurement_info) :: demoinfo 
  integer                 :: ninfo

!  read nc file
  integer :: ncid,inn_id
  integer :: nobs 
  integer,allocatable :: sip(:),sjp(:)
  real,   allocatable :: svar(:),sdep(:), forevar(:)
  real,   allocatable :: oforvar(:),osvar(:)
  real,   allocatable :: innova(:), oinnovar(:)
  
  
  integer, external :: iargc
  character(512) :: option1,option2
  character*100  :: fileuf,filenc
  character*100  :: fileout

  real           :: VarScale

  if (iargc() /= 2) then
     print *, 'Usage: obs_info <all profile-uf>  <target profile.nc> '
     print *, 'Options:'
     print *, " obs_info first-uf last-nc"
     stop
  else
    call getarg(1, option1)
    call getarg(2, option2)
  end if
  fileuf=trim(option1)
  filenc=trim(option2)

  inquire(file =trim(fileuf), exist = exists)
  if (.not. exists) then
    print *, 'ERROR: file "????info.uf" does not exist'
    stop
  end if
  inquire(file =trim(filenc), exist = exists)
  if (.not. exists) then
    print *, 'ERROR: file "observation????.nc" does not exist'
    stop
  end if
 
  ! read the profile infomation
  ninfo=0;
  inquire(iolength = rsize) demoinfo 
  open(10, file = trim(fileuf), form = 'unformatted',&
         access = 'direct', recl = rsize, status = 'old')
    ! I guess there is no other way to work out the length other than read the
    ! file in fortran - PS
    k =1 
    do while (.true.)
       read(10, rec = k, iostat = ios) demoinfo 
       if (ios /= 0) then
          nobs = k - 1
          exit
       end if
       k = k + 1
    enddo
    print *,nobs,rsize
    if(k>1) then
      ninfo=nobs
      allocate(obsinfo(ninfo))
      do l=1,ninfo
        read(10,rec=l) obsinfo(l)
      end do
    endif
  close(10)  

  if(ninfo>0) then
    if(obsinfo(1)%id .eq. 'TEM') then
      VarScale=Tvar0
    elseif(obsinfo(1)%id .eq. 'SAL') then
      VarScale=Svar0
    endif
    ! read the thinning observed profile
    call nfw_open(filenc,nf_nowrite,ncid);
    if(nfw_dim_exists(ncid,'nsobs')) then
      call nfw_inq_dimid(filenc,ncid,'nsobs',inn_id)
      call nfw_inq_dimlen(filenc, ncid, inn_id, nobs)
      print *, nobs
      allocate(sip(nobs),sjp(nobs),svar(nobs),sdep(nobs))
      allocate(forevar(nobs))
      allocate(innova(nobs))
      call nfw_inq_varid(filenc, ncid, 'sipiv', inn_id)
      call nfw_get_var_int(filenc, ncid, inn_id, sip)
      call nfw_inq_varid(filenc, ncid, 'sjpiv', inn_id)
      call nfw_get_var_int(filenc, ncid, inn_id, sjp)
      call nfw_inq_varid(filenc, ncid, 'sdepth', inn_id)
      call nfw_get_var_double(filenc, ncid, inn_id, sdep)
      call nfw_inq_varid(filenc, ncid, 'svar', inn_id)
      call nfw_get_var_double(filenc, ncid, inn_id, svar)
      call nfw_inq_varid(filenc, ncid, 'forecast_variance', inn_id)
      call nfw_get_var_double(filenc, ncid, inn_id, forevar)
      call nfw_inq_varid(filenc, ncid, 'innovation', inn_id)
      call nfw_get_var_double(filenc, ncid, inn_id, innova)
    else
      print *, 'No NC file after thinning'
      stop
    endif

! compare the svar with 0.5 replaced by forecast variance 5 and 9 times
    allocate(oforvar(ninfo),osvar(ninfo),oinnovar(ninfo))
    oforvar=0; osvar=0; oinnovar=0;
    obsinfo(:)%signal=0
    do k=1,nobs
      where(obsinfo%ipiv==sip(k) .and. obsinfo%jpiv==sjp(k)) obsinfo%signal=1 
 !     where(obsinfo%ipiv==sip(k) .and. obsinfo%jpiv==sjp(k) &
 !           .and. svar(k)>=5*max(,forevar(k))) obsinfo%signal=2 
      where(obsinfo%ipiv==sip(k) .and. obsinfo%jpiv==sjp(k) & 
            .and. svar(k)>=3*max(VarScale,forevar(k))) 
        obsinfo%signal=3
        oforvar=forevar(k)
        osvar=svar(k)
        oinnovar=innova(k)

      end where
    end do
!    fileout="check_Prof_diag.txt"
!    open(10,file=trim(fileout))
!      do k=1,ninfo
!         if(obsinfo(k)==0) then
!
!
!         endif
!      end do 
!    close(10)
! 
!    stop

    fileout='diag.txt'
    open(10,file=trim(fileout)) 

      k0=0      
      do k=1,ninfo

      if(obsinfo(k)%signal/=1) then
         !print *,obsinfo(k)%signal, ' ',trim(obsinfo(k)%platnum),' ',trim(obsinfo(k)%inifile) 
         if(k0.eq.0) then
           write(10,'(a10,a7,2i5,2F9.3,F7.1,i3,3F9.3,2x,a)') trim(obsinfo(k)%platnum), &
               obsinfo(k)%id, obsinfo(k)%ipiv, obsinfo(k)%jpiv, &
               obsinfo(k)%lon,obsinfo(k)%lat, & 
               obsinfo(k)%maxdep, obsinfo(k)%signal,oforvar(k),osvar(k),oinnovar(k),  &
               trim(obsinfo(k)%inifile)
           demoinfo=obsinfo(k)
           k0=k0+1
         else
           if(demoinfo%platnum/=obsinfo(k)%platnum.and.demoinfo%ipiv/=obsinfo(k)%ipiv &
              .and.demoinfo%jpiv/=obsinfo(k)%jpiv) then
              write(10,'(a10,a7,2i5,2F9.3,F7.1,i3,3F9.3,2x,a)') trim(obsinfo(k)%platnum), & 
                 obsinfo(k)%id, obsinfo(k)%ipiv, obsinfo(k)%jpiv, &
                 obsinfo(k)%lon,obsinfo(k)%lat, & 
                 obsinfo(k)%maxdep, obsinfo(k)%signal,oforvar(k),osvar(k),oinnovar(k), &
                 trim(obsinfo(k)%inifile)
                 demoinfo=obsinfo(k)
           endif
         endif
      endif

      end do
    close(10)
    inquire(file=trim(fileout),exist=exists)
    if(exists) then
      call system("sort diag.txt | uniq > diag_sort.txt")
    endif

    deallocate(sip,sjp,svar,sdep)
    deallocate(forevar,innova)
    deallocate(obsinfo)
    deallocate(oforvar,osvar,oinnovar)
  endif

end program p_obs_info
