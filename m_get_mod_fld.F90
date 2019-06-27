module m_get_mod_fld
! KAL -- This routine reads one of the fields from the model, specified
! KAL -- by name, vertical level and time level 
! KAL -- This routine is really only effective for the new restart files.

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine get_mod_fld(fld,j,cfld,vlevel,tlevel,nx,ny)
#if defined (QMPI)
   use qmpi
#else
   use qmpi_fake
#endif
   implicit none
   integer,      intent(in)            :: nx,ny  ! Grid dimension
   integer,      intent(in)            :: j      ! Ensemble member to read
   real, dimension(nx,ny), intent(out) :: fld    ! output fld
   character(len=*), intent(in)        :: cfld   ! name of fld
   integer, intent(in)                 :: tlevel ! time level
   integer, intent(in)                 :: vlevel ! vertical level

   integer reclICE
   real*8, dimension(nx,ny) :: ficem,hicem,hsnwm,ticem,tsrfm

   logical ex

   character(len=*),parameter :: icefile='forecastICE.uf'

   ! KAL -- shortcut -- the analysis is for observation icec -- this little "if" 
   ! means the  analysis will only work for ice. Add a check though
   if ((trim(cfld)/='icec' .and. trim(cfld)/='hice')  .or. vlevel/=0 .or. tlevel/=1)then
      if (master) print *,'get_mod_fld only works for icec for now'
      call stop_mpi()
   end if

!###################################################################
!####################### READ  ICE   MODEL #########################
!###################################################################
#if defined (ICE)
#warning "COMPILING WITH ICE"
   inquire(exist=ex,file=icefile)
   if (.not.ex) then
      if (master) then
         print *,icefile//' does not exist!'
         print *,'(get_mod_fld)'
      end if
      call stop_mpi()
   end if
   inquire(iolength=reclICE)ficem,hicem,hsnwm,ticem,tsrfm  !,iceU,iceV
   open(10,file=icefile,form='unformatted',access='direct',recl=reclICE,action='read')
      read(10,rec=j)ficem,hicem,hsnwm,ticem,tsrfm !,iceU,iceV
      if (trim(cfld)=='icec') fld = ficem
      if (trim(cfld)=='hice') fld = hicem
   close(10)
#else
#warning "COMPILING WITHOUT ICE"
#endif


  return
end subroutine get_mod_fld



! KAL - This is for the new file type
subroutine get_mod_fld_new(memfile,fld,iens,cfld0,vlevel,tlevel,nx,ny)
   use mod_raw_io
#if defined (QMPI)
   use qmpi, only : qmpi_proc_num, master
#else
   use qmpi_fake
#endif
   implicit none
   integer,      intent(in)            :: nx,ny  ! Grid dimension
   integer,      intent(in)            :: iens   ! Ensemble member to read
   real, dimension(nx,ny), intent(out) :: fld    ! output fld
   character(len=*), intent(in)        :: memfile! base name of input files
   character(len=*), intent(in)        :: cfld0   ! name of fld
   integer, intent(in)                 :: tlevel ! time level
   integer, intent(in)                 :: vlevel ! vertical level

   real*8, dimension(nx,ny) :: readfldr8
   real*4, dimension(nx,ny) :: readfldr4
   real*4:: amin, amax,spval
   real :: bmin, bmax
   integer :: indx
!----------------------------------- 
   character*80 :: cfld
   cfld=trim(cfld0)
#if defined (SINGLE_RESTART)
   if (trim(cfld0) == 'icec') then
     cfld='ficem'
   elseif (trim(cfld0) == 'hice') then
     cfld='hicem'
   endif
#endif 
!----------------------------------- 
   ! Dette fordi is-variablane forelobig er paa gammalt format.
   if (trim(cfld) /= 'icec' .and. trim(cfld) /= 'hice') then

      ! KAL - 1) f kva index som skal lesast finn vi fraa .b fil (header)
      call rst_index_from_header(trim(memfile)//'.b', & ! filnavn utan extension
                                 cfld               , & ! felt som skal lesast fex saln,temp
                                 vlevel,              & ! vertikalnivaa
                                 tlevel,              & ! time level - kan vere 1 eller 2 - vi bruker 1 foreloepig
                                 indx,                & ! indexen som maa lesas fra data fila
                                 bmin,bmax,           & ! min/max - kan sjekkast mot det som er i datafila
                                 .true. )

      if (indx < 0) then
         if (master) then
            print *, 'ERROR: get_mod_fld_new(): ', trim(memfile), '.b: "',&
                 trim(cfld), '" not found'
         end if
         stop
      end if

      ! KAL -- les datafelt vi fann fraa header fila (indx)
      spval=0.
      call READRAW(readfldr4          ,& ! Midlertidig felt som skal lesast
                   amin, amax         ,& ! max/min fraa data (.a) fila 
                   nx,ny              ,& ! dimensjonar
                   .false.,spval      ,& ! dette brukast for  sette "no value" verdiar
                   trim(memfile)//'.a',& ! fil som skal lesast fraa
                   indx)                 ! index funne over

     ! Sjekk p at vi har lest rett - samanlign max/min fr filene
     if     (abs(amin-bmin).gt.abs(bmin)*1.e-4 .or. &
             abs(bmax-amax).gt.abs(bmax)*1.e-4     ) then
        print *,'Inconsistency between .a and .b files'
        print *,'.a : ',amin,amax
        print *,'.b : ',bmin,bmax
        print *,cfld,vlevel,tlevel
        print *,indx
        print *,'node ',qmpi_proc_num
        call exit(1)
     end if
     fld=readfldr4

   else ! fld = fice, hice
      ! Gammal rutine ja
      call get_mod_fld(readfldr8,iens,cfld,0,1,nx,ny)
      fld=readfldr8
   end if


end subroutine



end module m_get_mod_fld


