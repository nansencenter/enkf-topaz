module m_get_mod_fld_nc
! JPX -- This routine reads one of the fields from the model NC file,
! specified by name, vertical (categolory), time, and return the integer flag
! to indicate the reading is valid or not
! Currently this routine is effective for forecast.nc and
! ice_forecast.nc
   use nfw_mod
contains



subroutine get_mod_fld_nc(memfile,fld,cfld0,vlevel,tlevel,nx,ny)
#if defined (QMPI)
   use qmpi, only : qmpi_proc_num, master, stop_mpi
#else
   use qmpi_fake
#endif
   use nfw_mod
   implicit none
   integer,      intent(in)            :: nx,ny    ! Grid dimension
   real, dimension(nx,ny), intent(out) :: fld      ! output fld
   character(len=*), intent(in)        :: memfile  ! base name of input files
   character(len=*), intent(in)        :: cfld0    ! name of fld
   integer, intent(in)                 :: vlevel   ! vertical level  % 0 on surface
   integer, intent(inout)              :: tlevel   ! time level (also works as valid flag)

   real*8, dimension(nx,ny) :: readfldr8
   real*4, dimension(nx,ny) :: readfldr4

   integer              :: ncid
   integer              :: v_dimid, v_id
   integer              :: Ndim
   logical              :: Tex
   integer              :: vids(3)
   character*(NF_MAX_NAME) :: name
   integer              :: status,len1, len2, len3

   character*80 :: cfld

   cfld=trim(cfld0)
   call nfw_open(memfile,nf_nowrite,ncid)
   Tex=nfw_var_exists(ncid,cfld)
   if (Tex) then 
      call nfw_inq_varid(memfile, ncid, cfld, v_id)
      !print *,'Existing '//trim(memfile)//' the multi-cate.: ', trim(cfld), ' is ', v_id, vlevel
      if (vlevel==0) then   ! on surface
        call nfw_get_var_double(memfile,ncid,v_id,readfldr8)
        tlevel=1
        fld=readfldr8
      else                  ! for different category
        call nfw_get_var_double2D(memfile,ncid,v_id,vlevel,readfldr8,nx,ny)
        tlevel=1
        fld=readfldr8
      endif
   end if
   call nfw_close(memfile,ncid)
end subroutine

end module m_get_mod_fld_nc

