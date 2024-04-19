!KAL -- this module allows us to fine-tune the fields
!KAL -- we wish to include in tha analysis. The new
!KAL -- layout of the EnKF makes it possible to specify fields
!KAL -- to analyze at run-time rather than at compile-time
!KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!KAL --
!KAL -- Module variables:
!KAL --    numfields   - total number of fields to process
!KAL --    fieldnames  - the names of the fields we wish to analyze
!KAL --    fieldlevel  - the levels of the associated fields
!JPX --    numfields_hycom  - recording the total number of hycom file
!KAL --
!KAL -- Ex: If we only want to assimilate temperatures in layer
!KAL --     one and two, numfields, fieldnames and fieldlevel 
!KAL --     would look like:
!KAL --
!KAL --     numfields=2                                 
!KAL --     fieldnames (1)='temp', fieldnames (2)='temp'
!KAL --     fieldlevel (1)=     1, fieldlevel (2)=2
!KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!KAL -- The file "analysisfields.in" specifies the fields to 
!KAL -- inlude in the analysis. Format of one line is fieldname
!KAL -- first layer and last layer, example
!KAL --
!KAL -- fldname   1 22
!KAL -- 12345678901234567890123456789012345678901234567890
!KAL --
!KAL -- Fortran format for one line is '(a8,2i3)'
!KAL --
!KAL -- Example: to specify that we want temperature and salinity 
!KAL --          in layers 1..22 to be updated, as well as 
!KAL --          ice concentration (layer 0), specify:
!KAL --
!KAL -- saln      1 22
!KAL -- temp      1 22
!KAL -- hice      0  0
!KAL -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

module mod_analysisfields

character(len=*), parameter :: infile='analysisfields.in'
integer :: numfields
character(len=8), dimension(:), allocatable:: fieldnames
integer         , dimension(:), allocatable:: fieldlevel 
integer         , dimension(:), allocatable:: fieldindex

integer                     :: numfields_hycom
#if defined (HYCOM_CICE)
character(len=*), parameter :: infile_ice='analysisfields_ice.in'
#endif

contains

   integer function get_nrfields()
#if defined (QMPI)
   use qmpi
#else
   use qmpi_fake
#endif
   implicit none
   integer :: ios,first,last
   logical :: ex
   character(len=8) :: char8

   inquire(exist=ex,file=infile)
   if (.not. ex) then
      if (master) print *,'Could not find '//infile
      call stop_mpi()
   end if

   open(10,status='old',form='formatted',file=infile,action='read')
   ios=0
   get_nrfields=0
   do while (ios==0)
      read(10,100,iostat=ios) char8,first,last
      if (ios==0) then
         get_nrfields=get_nrfields+last-first+1
      end if
   end do
   close(10)
   100 format (a8,2i3)

   numfields_hycom=get_nrfields
#if defined (HYCOM_CICE)
   inquire(exist=ex,file=infile_ice)

!   if (.not. ex) then
!      if (master) print *,'Could not find '//infile_ice
!      call stop_mpi()
!   end if
   if (ex) then
     open(10,status='old',form='formatted',file=infile_ice,action='read')
     ios=0
     do while (ios==0)
       read(10,200,iostat=ios) char8,first,last
       if (ios==0) get_nrfields=get_nrfields+last-first+1
     end do
     close(10)
   else
     if (master) print *,'Could not find '//infile_ice
!      call stop_mpi()
   end if

   200 format (a8,2i3)
#endif

   end function

   subroutine get_analysisfields()
#if defined (QMPI)
   use qmpi
#else
   use qmpi_fake
#endif
   implicit none
   integer :: first,last,k,nfld,ios
   logical :: ex
   character(len=8) :: char8
   integer          :: k0

   numfields=get_nrfields()
   if (master) print *,'numfields is ',numfields
   if (numfields<=0 .or.numfields > 600) then !
      if (master) print *,'numfields is higher than max allowed setting or = 0'
      call stop_mpi()
   end if
   allocate(fieldnames(numfields))
   allocate(fieldlevel(numfields))
   allocate(fieldindex(numfields))


   inquire(exist=ex,file=infile)
   if (.not. ex) then
      if (master) print *,'Could not find '//infile
      call stop_mpi()
   end if

   open(10,status='old',form='formatted',file=infile,action='read')
   k0=0
   ios=0
   nfld=0
   do while (ios==0)
      read(10,100,iostat=ios) char8,first,last
      if (ios==0) then
         do k=first,last
            fieldnames (nfld+k-first+1)=char8
            fieldlevel (nfld+k-first+1)=k
            k0=k0+1
            fieldindex (nfld+k-first+1)=k0
         end do
         nfld=nfld+last-first+1
      end if
   end do
   close(10)
   100 format (a8,2i3)


#if defined (HYCOM_CICE)
   inquire(exist=ex,file=infile_ice)
   if (ex) then
     open(10,status='old',form='formatted',file=infile_ice,action='read')
     ios=0
     do while (ios==0)
        read(10,200,iostat=ios) char8,first,last
        if (ios==0) then
           do k=first,last
              fieldnames (nfld+k-first+1)=char8
              fieldlevel (nfld+k-first+1)=k
              k0=k0+1
              fieldindex (nfld+k-first+1)=k0
           end do
           nfld=nfld+last-first+1
        end if
     end do
     close(10)
   else
      if (master) print *,'Could not find '//infile_ice
    !  call stop_mpi()
   end if
   200 format (a8,2i3)
#endif

   if (nfld/=numfields) then
      if (master) print *,'An error occured when reading '//infile, ' ', nfld
      call stop_mpi()
   end if

   ! List fields used in analysis
   do k=1,numfields
      if (master) print *,fieldnames(k),fieldlevel(k),fieldindex(k)
   end do

   end subroutine
end module mod_analysisfields






