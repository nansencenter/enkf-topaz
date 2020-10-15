!
! File: nfw.f90
!
! Author: Pavel Sakov, CSIRO Marine Research
!
! Created: 17 March 2005
!
! Purpose: Contains wrappers to netcdf functions, mainly for easier
!          error handling.
! 
! Description:
!
!          Each subroutine in nfw.f90 is a simple wrapper of a similar
!          function in the NetCDF Fortran interface. The rules of use are
!          pretty simple: for a given NetCDF Fortran function, replace
!          prefix "nf_" by "nfw_" and add the NetCDF file name as the
!          first argument.
!
!          Here is the current list of subroutines in nfw_mod:
!
!          nfw_create(fname, mode, ncid)
!          nfw_open(fname, mode, ncid)
!          nfw_enddef(fname, ncid)
!          nfw_close(fname, ncid)
!          nfw_inq_unlimdim(fname, ncid, unlimdimid)
!          nfw_inq_dimid(fname, ncid, name, dimid)
!          nfw_inq_dimlen(fname, ncid, dimid, length)
!          nfw_def_dim(fname, ncid, name, length, dimid)
!          nfw_def_var(fname, ncid, name, type, ndims, dimids, varid)
!          nfw_inq_varid(fname, ncid, name, varid)
!          nfw_inq_varname(fname, ncid, varid, name)
!          nfw_inq_varndims(fname, ncid, varid, ndims)
!          nfw_inq_vardimid(fname, ncid, varid, dimids)
!          nfw_rename_var(fname, ncid, oldname, newname)
!          nfw_put_var_int(fname, ncid, varid, v)
!          nfw_put_var_double(fname, ncid, varid, v)
!          nfw_put_var_real(fname, ncid, varid, v)
!          nfw_get_var_int(fname, ncid, varid, v)
!          nfw_get_var_double(fname, ncid, varid, v)
!          nfw_put_vara_int(fname, ncid, varid, start, length, v)
!          nfw_put_vara_double(fname, ncid, varid, start, length, v)
!          nfw_get_vara_int(fname, ncid, varid, start, length, v)
!          nfw_get_vara_double(fname, ncid, varid, start, length, v)
!          nfw_get_att_int(fname, ncid, varid, attname, v)
!          nfw_get_att_real(fname, ncid, varid, attname, v)
!          nfw_get_att_double(fname, ncid, varid, attname, v)
!          nfw_put_att_text(fname, ncid, varid, attname, length, text)
!          nfw_put_att_int(fname, ncid, varid, attname, type, length, v)
!          nfw_put_att_real(fname, ncid, varid, attname, type, length, v)
!          nfw_put_att_double(fname, ncid, varid, attname, type, length, v)
!
!          Derived procedures:
!
!          nfw_get_var_double_firstrecord(fname, ncid, varid, v)
!          nfw_var_exists(ncid, name)
!          nfw_dim_exists(ncid, name)
! Modifications:
!
! 29/04/2008 PS: added nfw_rename_var(fname, ncid, oldname, newname)
! 21/10/2009 PS: added nfw_var_exists(ncid, name)
! 22/10/2009 PS: added nfw_put_att_double(fname, ncid, varid, attname, type, 
!
! 06/11/2009 PS: added nfw_dim_exists(ncid, name)
! 20/09/2019 PS: added nfw_var_att_exists(ncid, varid, attname, type,length)
!                nfw_put_att_real(fname, ncid, varid, attname, type, length, v)
!                nfw_get_att_real(fname, ncid, varid, attname, v)

module nfw_mod
  implicit none
  include 'netcdf.inc'

  character(*), private, parameter :: nfw_version = "0.03"
  integer, private, parameter :: logunit = 6
  character(*), private, parameter :: errprefix = "nfw: error:          "
  private quit1, quit2, quit3

contains

#if defined(F90_NOFLUSH)
  subroutine flush(dummy)
    integer, intent(in) :: dummy
  end subroutine flush
#endif

  ! Common exit point -- for the sake of debugging
  subroutine quit
    stop
  end subroutine quit

  subroutine quit1(fname, procname, status)
    character*(*), intent(in) :: fname
    character*(*), intent(in) :: procname
    integer, intent(in) :: status
 
    write(logunit, *)
    write(logunit, *) errprefix, '"', trim(fname), '": ', procname, '(): ',&
         nf_strerror(status)
    call flush(logunit)
    call quit
  end subroutine quit1

  subroutine quit2(fname, procname, name, status)
    character*(*), intent(in) :: fname
    character*(*), intent(in) :: procname
    character*(*), intent(in) :: name
    integer, intent(in) :: status

    write(logunit, *)
    write(logunit, *) errprefix, '"', trim(fname), '": ', procname, '(): "',&
         trim(name), '": ', nf_strerror(status)
    call flush(logunit)
    call quit
  end subroutine quit2

  subroutine quit3(fname, procname, name1, name2, status)
    character*(*), intent(in) :: fname
    character*(*), intent(in) :: procname
    character*(*), intent(in) :: name1
    character*(*), intent(in) :: name2
    integer, intent(in) :: status

    write(logunit, *)
    write(logunit, *) errprefix, '"', trim(fname), '": ', procname, '(): "',&
         trim(name1), '": "', trim(name2), '": ', nf_strerror(status)
    call flush(logunit)
    call quit
  end subroutine quit3

  subroutine nfw_create(fname, mode, ncid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: mode
    integer, intent(out) :: ncid

    integer :: status

    status = nf_create(trim(fname), mode, ncid)
    if (status /= 0) call quit1(fname, 'nf_create', status)
  end subroutine nfw_create

  subroutine nfw_open(fname, mode, ncid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: mode
    integer, intent(out) :: ncid

    integer :: status

    status = nf_open(trim(fname), mode, ncid)
    if (status /= 0) call quit1(fname, 'nf_open', status)
  end subroutine nfw_open

  subroutine nfw_enddef(fname, ncid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid

    integer :: status

    status = nf_enddef(ncid)
    if (status /= 0) call quit1(fname, 'nf_enddef', status)
  end subroutine nfw_enddef

  subroutine nfw_redef(fname, ncid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid

    integer :: status

    status = nf_redef(ncid)
    if (status /= 0) call quit1(fname, 'nf_redef', status)
  end subroutine nfw_redef

  subroutine nfw_close(fname, ncid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid

    integer :: status

    status = nf_close(ncid)
    if (status /= 0) call quit1(fname, 'nf_close', status)
  end subroutine nfw_close

  subroutine nfw_inq_unlimdim(fname, ncid, unlimdimid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(out) :: unlimdimid

    integer :: status
    
    status = nf_inq_unlimdim(ncid, unlimdimid)
    if (status /= 0) call quit1(fname, 'nf_inq_unlimdimid', status)
  end subroutine nfw_inq_unlimdim

  subroutine nfw_inq_dimid(fname, ncid, name, dimid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    character*(*), intent(in) :: name
    integer, intent(out) :: dimid

    integer :: status
    
    status = nf_inq_dimid(ncid, trim(name), dimid)
    if (status /= 0) call quit2(fname, 'nf_inq_dimid', name, status)
  end subroutine nfw_inq_dimid

  subroutine nfw_inq_dimlen(fname, ncid, dimid, length)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: dimid
    integer, intent(out) :: length

    integer :: status

    status = nf_inq_dimlen(ncid, dimid, length)
    if (status /= 0) call quit1(fname, 'nf_inq_dimlen', status)
  end subroutine nfw_inq_dimlen

  subroutine nfw_def_dim(fname, ncid, name, length, dimid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    character*(*), intent(in) :: name
    integer, intent(in) :: length
    integer, intent(out) :: dimid

    integer :: status

    status = nf_def_dim(ncid, name, length, dimid)
    if (status /= 0) call quit2(fname, 'nf_def_dim', name, status)
  end subroutine nfw_def_dim

  subroutine nfw_def_var(fname, ncid, name, type, ndims, dimids, varid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    character*(*), intent(in) :: name
    integer, intent(in) :: type
    integer, intent(in) :: ndims
    integer, intent(in) :: dimids(*)
    integer, intent(out) :: varid

    integer :: status

    status = nf_def_var(ncid, name, type, ndims, dimids, varid)
    if (status /= 0) call quit2(fname, 'nf_def_var', name, status)
  end subroutine nfw_def_var

  subroutine nfw_inq_varid(fname, ncid, name, varid)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    character*(*), intent(in) :: name
    integer, intent(out) :: varid

    integer :: status
   
    status = nf_inq_varid(ncid, trim(name), varid)
    if (status /= 0) call quit2(fname, 'nf_inq_varid', name, status)
  end subroutine nfw_inq_varid

  subroutine nfw_inq_varname(fname, ncid, varid, name)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(out) :: name

    integer :: status

    status = nf_inq_varname(ncid, varid, name)
    if (status /= 0) call quit1(fname, 'nf_inq_varname', status)
  end subroutine nfw_inq_varname

  subroutine nfw_inq_varndims(fname, ncid, varid, ndims)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(out) :: ndims

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_inq_varndims(ncid, varid, ndims)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_inq_varndims', name, status)
    end if
  end subroutine nfw_inq_varndims

  subroutine nfw_inq_vardimid(fname, ncid, varid, dimids)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(out) :: dimids(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_inq_vardimid(ncid, varid, dimids)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_inq_vardimid', name, status)
    end if
  end subroutine nfw_inq_vardimid


  subroutine nfw_inq_var_strlen(fname,ncid,varid,Sn1,Sn2)
    character*(*),  intent(in) :: fname
    integer,   intent(in)      :: ncid
    integer,   intent(in)      :: varid
    integer,      intent(out)  :: Sn1,Sn2 
 
    integer ::  ndims,dimid,len_var    
    integer :: dimids(NF_MAX_VAR_DIMS)
    integer :: status
    character*(NF_MAX_NAME) :: name

    call nfw_inq_varndims(fname, ncid, varid, ndims)
    if(ndims==1)then
       call nfw_inq_vardimid(fname,ncid,varid,dimids)
       dimid=dimids(1)
       status=nf_inq_dim(ncid,dimid,name,len_var)
       if (status /= 0) then
         write(*,*) 'Error in nfw_inq_var_strlen'
         call quit2(fname, 'nf_inq_dim', name, status)
       endif
       Sn1=len_var
       Sn2=1 
    elseif(ndims==2) then
       call nfw_inq_vardimid(fname,ncid,varid,dimids)
       dimid=dimids(1)
       status=nf_inq_dim(ncid,dimid,name,len_var)
       if (status /= 0) then
         write(*,*) 'Error in nfw_inq_var_strlen'
         call quit2(fname, 'nf_inq_var_dimlen-1', name, status)
       endif
       Sn1=len_var
       dimid=dimids(2)
       status=nf_inq_dim(ncid,dimid,name,len_var)
       if (status /= 0) then
         write(*,*) 'Error in nfw_inq_var_strlen'
         call quit2(fname, 'nf_inq_var_dimlen-2', name, status)
       endif
       Sn2=len_var
    endif
  end subroutine nfw_inq_var_strlen



  subroutine nfw_get_var_strings(fname,ncid,varid,varstring,Sn1,Sn2)
    character*(*),  intent(in) :: fname
    integer,   intent(in)      :: ncid
    integer,   intent(in)      :: varid
    integer,   intent(in)      :: Sn1,Sn2
    character, intent(out)     :: varstring (Sn1,Sn2)
 
    integer ::  ndims,dimid,len_var    
    integer :: dimids(NF_MAX_VAR_DIMS)
    integer :: status
    character*(NF_MAX_NAME) :: name

    call nfw_inq_varndims(fname, ncid, varid, ndims)

    if(ndims==1)then
       status = nf_get_var_text(ncid,varid,varstring)
       if (status /= 0) then
         call nfw_inq_varname(fname, ncid, varid, name)
         call quit2(fname, 'nf_get_var_strings', name, status)
       endif
    elseif(ndims==2) then
       status=nf_get_var_text(ncid,varid,varstring)
       !status=nf_get_var_text(ncid,varid,[1 1],[Sn1 Sn2],varstring)
       if (status /= 0) then
         call nfw_inq_varname(fname, ncid, varid, name)
         call quit2(fname, 'nf_get_var_strings', name, status)
       endif
    endif
  end subroutine

  subroutine nfw_rename_var(fname, ncid, oldname, newname)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    character*(*), intent(in) :: oldname
    character*(*), intent(in) :: newname

    integer :: varid
    integer :: status

    call nfw_inq_varid(fname, ncid, oldname, varid)
    status = nf_rename_var(ncid, varid, newname)
    if (status /= 0) then
       call quit2(fname, 'nf_rename_var', oldname, status)
    end if
  end subroutine nfw_rename_var

  subroutine nfw_put_var_int(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_put_var_int(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_put_var_double', name, status)
    end if
  end subroutine nfw_put_var_int

  subroutine nfw_put_var_double(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(8), intent(in) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_put_var_double(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_put_var_double', name, status)
    end if
  end subroutine nfw_put_var_double

  subroutine nfw_put_var_real(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(4), intent(in) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_put_var_real(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_put_var_real', name, status)
    end if
  end subroutine nfw_put_var_real

  subroutine nfw_get_var_int(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(out) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_get_var_int(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_var_int', name, status)
    end if
  end subroutine nfw_get_var_int

  subroutine nfw_get_var_int1(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer*1, intent(out) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_get_var_int1(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_var_int1', name, status)
    end if
  end subroutine nfw_get_var_int1

  subroutine nfw_get_var_int2(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer*2, intent(out) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_get_var_int2(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_var_int2', name, status)
    end if
  end subroutine nfw_get_var_int2

  subroutine nfw_get_var_real(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(4), intent(out) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_get_var_real(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_var_real', name, status)
    end if
  end subroutine nfw_get_var_real

  subroutine nfw_get_var_double(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(8), intent(out) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_get_var_double(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_var_double', name, status)
    end if
  end subroutine nfw_get_var_double

  subroutine nfw_get_var_text(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character, intent(out) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_get_var_text(ncid, varid, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_var_int', name, status)
    end if
  end subroutine nfw_get_var_text

  subroutine nfw_put_vara_int(fname, ncid, varid, start, length, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: start(*)
    integer, intent(in) :: length(*)
    integer, intent(in) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_put_vara_int(ncid, varid, start, length, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_put_vara_int', name, status)
    end if
  end subroutine nfw_put_vara_int

  subroutine nfw_put_vara_double(fname, ncid, varid, start, length, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: start(*)
    integer, intent(in) :: length(*)
    real(8), intent(in) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_put_vara_double(ncid, varid, start, length, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_put_vara_double', name, status)
    end if
  end subroutine nfw_put_vara_double

  subroutine nfw_get_vara_int(fname, ncid, varid, start, length, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: start(*)
    integer, intent(in) :: length(*)
    integer, intent(out) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_get_vara_int(ncid, varid, start, length, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_vara_int', name, status)
    end if
  end subroutine nfw_get_vara_int

  subroutine nfw_get_vara_double(fname, ncid, varid, start, length, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: start(*)
    integer, intent(in) :: length(*)
    real(8), intent(out) :: v(*)

    character*(NF_MAX_NAME) :: name
    integer :: status

    status = nf_get_vara_double(ncid, varid, start, length, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_vara_double', name, status)
    end if
  end subroutine nfw_get_vara_double

  subroutine nfw_get_att_int(fname, ncid, varid, attname, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(in) :: attname
    integer, intent(out) :: v(*)

    character*(NF_MAX_NAME) :: varname
    integer :: status

    status = nf_get_att_int(ncid, varid, attname, v)
    if (status /= 0) then
       if (varid /= nf_global) then
          call nfw_inq_varname(fname, ncid, varid, varname)
       else
          varname = 'NF_GLOBAL'
       end if
       call quit3(fname, 'nf_get_att_int', varname, attname, status)
    end if
  end subroutine nfw_get_att_int

  subroutine nfw_get_att_real(fname, ncid, varid, attname, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(in) :: attname
    real(4), intent(out) :: v(*)

    character*(NF_MAX_NAME) :: varname
    integer :: status

    status = nf_get_att_real(ncid, varid, attname, v)
    if (status /= 0) then
       if (varid /= nf_global) then
          call nfw_inq_varname(fname, ncid, varid, varname)
       else
          varname = 'NF_GLOBAL'
       end if
       call quit3(fname, 'nf_get_att_real', varname, attname, status)
    end if
  end subroutine nfw_get_att_real

  subroutine nfw_get_att_double(fname, ncid, varid, attname, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(in) :: attname
    real(8), intent(out) :: v(*)

    character*(NF_MAX_NAME) :: varname
    integer :: status

    status = nf_get_att_double(ncid, varid, attname, v)
    if (status /= 0) then
       if (varid /= nf_global) then
          call nfw_inq_varname(fname, ncid, varid, varname)
       else
          varname = 'NF_GLOBAL'
       end if
       call quit3(fname, 'nf_get_att_double', varname, attname, status)
    end if
  end subroutine nfw_get_att_double

  subroutine nfw_put_att_text(fname, ncid, varid, attname, length, text)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(in) :: attname
    integer, intent(in) :: length
    character*(*), intent(in) :: text

    integer :: status
    character*(NF_MAX_NAME) :: varname

    status = nf_put_att_text(ncid, varid, attname, length, trim(text))
    if (status /= 0) then
       if (varid /= nf_global) then
          call nfw_inq_varname(fname, ncid, varid, varname)
       else
          varname = 'NF_GLOBAL'
       end if
       call quit3(fname, 'nf_put_att_text', varname, attname, status)
    end if
  end subroutine nfw_put_att_text

  subroutine nfw_put_att_int(fname, ncid, varid, attname, type, length, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(in) :: attname
    integer, intent(in) :: type
    integer, intent(in) :: length
    integer, intent(in) :: v(*)

    integer :: status
    character*(NF_MAX_NAME) :: varname

    status = nf_put_att_int(ncid, varid, attname, type, length, v)
    if (status /= 0) then
       if (varid /= nf_global) then
          call nfw_inq_varname(fname, ncid, varid, varname)
       else
          varname = 'NF_GLOBAL'
       end if
       call quit3(fname, 'nf_put_att_int', varname, attname, status)
    end if
  end subroutine nfw_put_att_int

  subroutine nfw_put_att_real(fname, ncid, varid, attname, type, length, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(in) :: attname
    integer, intent(in) :: type
    integer, intent(in) :: length
    real(4), intent(in) :: v(*)

    integer :: status
    character*(NF_MAX_NAME) :: varname

    status = nf_put_att_real(ncid, varid, attname, type, length, v)
    if (status /= 0) then
       if (varid /= nf_global) then
          call nfw_inq_varname(fname, ncid, varid, varname)
       else
          varname = 'NF_GLOBAL'
       end if
       call quit3(fname, 'nf_put_att_real', varname, attname, status)
    end if
  end subroutine nfw_put_att_real

  subroutine nfw_put_att_double(fname, ncid, varid, attname, type, length, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    character*(*), intent(in) :: attname
    integer, intent(in) :: type
    integer, intent(in) :: length
    real(8), intent(in) :: v(*)

    integer :: status
    character*(NF_MAX_NAME) :: varname

    status = nf_put_att_double(ncid, varid, attname, type, length, v)
    if (status /= 0) then
       if (varid /= nf_global) then
          call nfw_inq_varname(fname, ncid, varid, varname)
       else
          varname = 'NF_GLOBAL'
       end if
       call quit3(fname, 'nf_put_att_double', varname, attname, status)
    end if
  end subroutine nfw_put_att_double

! Derived subroutines

  ! Reads the first record only
  subroutine nfw_get_var_double_firstrecord(fname, ncid, varid, v)
    character*(*), intent(in) :: fname
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    real(8), intent(out) :: v(*)

    integer :: ndims
    integer :: unlimdimid
    integer :: dimids(NF_MAX_VAR_DIMS)
    integer :: dimlen(NF_MAX_VAR_DIMS)
    integer :: dstart(NF_MAX_VAR_DIMS)
    integer :: i
    character*(NF_MAX_NAME) :: name
    integer :: status

    call nfw_inq_varndims(fname, ncid, varid, ndims)
    call nfw_inq_vardimid(fname, ncid, varid, dimids)
    call nfw_inq_unlimdim(fname, ncid, unlimdimid)
    
    do i = 1, ndims
       call nfw_inq_dimlen(fname, ncid, dimids(i), dimlen(i))
       dstart(i) = 1
    end do

    ! check size of v
    if (dimids(ndims) == unlimdimid) then
       dimlen(ndims) = 1 ! 1 record only
    end if

    status = nf_get_vara_double(ncid, varid, dstart, dimlen, v)
    if (status /= 0) then
       call nfw_inq_varname(fname, ncid, varid, name)
       call quit2(fname, 'nf_get_vara_double', name, status)
    end if
  end subroutine nfw_get_var_double_firstrecord

  logical function nfw_var_exists(ncid, name)
    integer, intent(in) :: ncid
    character*(*), intent(in) :: name

    integer :: varid
    integer :: status

    status = nf_inq_varid(ncid, trim(name), varid)
    nfw_var_exists = (status == 0)
  end function nfw_var_exists


  logical function nfw_var_att_exists(ncid,varid, name)
    integer, intent(in) :: ncid,varid
    character*(*), intent(in) :: name

    integer :: status
    integer :: cxtype,cnlen

    status = nf_inq_att(ncid, varid, trim(name),cxtype,cnlen)
    nfw_var_att_exists = (status == 0)
  end function nfw_var_att_exists


  logical function nfw_dim_exists(ncid, name)
    integer, intent(in) :: ncid
    character*(*), intent(in) :: name

    integer :: dimid
    integer :: status

    status = nf_inq_dimid(ncid, trim(name), dimid)
    nfw_dim_exists = (status == 0)
  end function nfw_dim_exists

end module nfw_mod
