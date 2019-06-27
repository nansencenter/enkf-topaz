module m_read_CLS_TSLA_grid
  ! Reads the CLS SST NetCDF dimensions

  integer, parameter, private :: STRLEN = 512

contains
  subroutine read_CLS_TSLA_grid(filename,gr)
    !use mod_dimensions
    use mod_grid
    use netcdf
    use nfw_mod
    implicit none

    character(len=80), intent(in) :: filename
    type(grid),        intent(out) :: gr
    character(len=80) :: fname
    logical :: ex
    !dimension ids
    integer :: data_ID,cycl_ID

    ! Array dimensions
    integer :: nb,cycl

    integer :: ncid,fcount
    character(STRLEN) :: Fpath
    character(STRLEN) :: ftemplate

    print *, 'read_CLS_TSLA_grid():'

    gr = default_grid
    gr%nx=0
    Fpath='./'
    ! Open file
    do fcount=1,7 !2 satellite Envissat,J2
       select case(fcount)
       case(1)
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_en*.nc'
       case(2)  
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_j1*.nc'
       case(3)
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_j2*.nc'
       case(4)  
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_e1*.nc'
       case(5)  
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_e2*.nc'
       case(6)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_tp*.nc'
       case(7)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_g2*.nc'
       end select
       call fname_fromtemplate(ftemplate, fname)
       inquire(file=trim(fname),exist=ex)
       if(ex) then
          call nfw_open(fname, nf_nowrite, ncid)
          print *, '  found "', trim(fname), '"...'

          ! Get dimension id in netcdf file ...
          call nfw_inq_dimid(fname, ncid, 'Data', data_ID)
          call nfw_inq_dimid(fname, ncid, 'Cycles', cycl_ID)
          ! Get dimension length from id
          call nfw_inq_dimlen(fname, ncid, data_ID, nb)
          call nfw_inq_dimlen(fname, ncid, cycl_ID, cycl)
          call nfw_close(fname, ncid)

          gr%nx=gr%nx+nb*cycl
          gr%ny=1
          gr%x0=0
          gr%y0=0
          gr%dx=0.1
          gr%dy=0.1
          gr%reg = .false.
          gr%order = 1
          gr%ux = 'm'
          gr%uy = 'm'
          gr%set = .true.
       endif
    enddo
  end subroutine read_CLS_TSLA_grid
  subroutine read_MYO_TSLA_grid(filename,gr)
    !use mod_dimensions
    use mod_grid
    use netcdf
    use nfw_mod
    implicit none

    character(len=80), intent(in) :: filename
    type(grid),        intent(out) :: gr
    character(len=80) :: fname
    logical :: ex
    !dimension ids
    integer :: time_ID

    ! Array dimensions
    integer :: nb

    integer :: ncid,fcount
    character(STRLEN) :: Fpath
    character(STRLEN) :: ftemplate

    print *, 'read_MYO_TSLA_grid():'

    gr = default_grid
    gr%nx=0
    Fpath='./'
    ! Open file
    do fcount=1,11   !2 satellite Envissat,J2
       select case(fcount)
       case(1)
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_en*.nc'
       case(2)  
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_j1*.nc'
       case(3)
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_j2*.nc'
       case(4)  
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_e1*.nc'
       case(5)  
          ftemplate=trim(Fpath)//'sla_'//trim(filename)//'_e2*.nc'
       case(6)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_tp*.nc'
       case(7)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_g2*.nc'
       case(8)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_c2*.nc'
       case(9)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_al*.nc'
       case(10)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_h2*.nc'
       case(11)
          ftemplate = trim(fpath)//'sla_'//trim(filename)//'_j3*.nc'
       end select
       call fname_fromtemplate(ftemplate, fname)
       inquire(file=trim(fname),exist=ex)
       if(ex) then
          call nfw_open(fname, nf_nowrite, ncid)
          print *, '  found "', trim(fname), '"...'

          ! Get dimension id in netcdf file ...
          call nfw_inq_dimid(fname, ncid, 'time', time_ID)
          ! Get dimension length from id
          call nfw_inq_dimlen(fname, ncid, time_ID, nb)
          call nfw_close(fname, ncid)

          gr%nx=gr%nx+nb
          gr%ny=1
          gr%x0=0
          gr%y0=0
          gr%dx=0.1
          gr%dy=0.1
          gr%reg = .false.
          gr%order = 1
          gr%ux = 'm'
          gr%uy = 'm'
          gr%set = .true.
       endif
    enddo
  end subroutine read_MYO_TSLA_grid
  
  

  subroutine fname_fromtemplate(ftemplate, fname)
    character(*), intent(in) :: ftemplate
    character(*), intent(inout) :: fname

    character(STRLEN) :: command ! (there may be a limit of 80 on some systems)
    integer :: ios

    command = 'ls '//trim(ftemplate)//' 2> /dev/null > tsla_files.txt'
    call system(trim(command));

    open(10, file = 'tsla_files.txt')
    read(10, fmt = '(a)', iostat = ios) fname
    close(10)
    if (ios /= 0) then
       fname = ""
    end if
  end subroutine fname_fromtemplate

end module m_read_CLS_TSLA_grid
