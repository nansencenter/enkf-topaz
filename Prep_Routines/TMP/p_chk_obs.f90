










! File :         p_chk_obs.F90
!
!                Created: unknown
!

program p_chk_obs
  use mod_measurement
  implicit none


  type (measurement), allocatable :: obs0(:)


  integer :: nobs, rsize,k
  logical :: exists = .false.
  type(measurement) :: record
  integer :: ios
  integer :: Iord
 
  integer, external :: iargc
  character(512) :: options

  if (iargc() /= 1) then
     print *, 'Usage: chk_obs <record number>'
     print *, '       chk_obs -h'
     print *, 'Options:'
     print *, '  -h -- describe parameter fie format'
     stop
  else
    call getarg(1, options)
    if (trim(options) == "-h") then
      print *, '  -h -- describe parameter fie format'
    end if
  end if


  
  k=1;


    inquire(file = 'observations.uf', exist = exists)
    if (.not. exists) then
      print *, 'ERROR: file "observations.uf" does not exist'
      stop
    end if
    inquire(iolength = rsize) record
    open(10, file = 'observations.uf', form = 'unformatted',&
         access = 'direct', recl = rsize, status = 'old')

    ! I guess there is no other way to work out the length other than read the
    ! file in fortran - PS
    !
    k = 1
    do while (.true.)
       read(10, rec = k, iostat = ios) record
       if (ios /= 0) then
          nobs = k - 1
          exit
       end if
       k = k + 1
    enddo
    print *,k

    call getarg(1, options)

    read(options,*) Iord  
    read(10,rec=Iord) record
    print *,Iord,record%d,record%id,record%date




end program p_chk_obs
