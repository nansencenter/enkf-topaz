module m_nf90_err
contains

   subroutine nf90_err(errcode,chars)
      use netcdf
      implicit none
      integer, intent(in) :: errcode
      character(len=*), optional :: chars
      character(len=80) :: hint


      hint =''
      if (present(chars)) hint=chars


      if (errcode/=NF90_NOERR) then
         write(6,'(a)') NF90_STRERROR(errcode)//'  '//trim(hint)
         stop '(handle_err)'
      end if

   end subroutine



end module m_nf90_err
