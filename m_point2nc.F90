! File:          m_point2nc.F90
!
! Created:       6 July 2010
!
! Last modified: 6/7/2010
!
! Author:        Pavel Sakov
!                NERSC
!
! Purpose:       Output of assimilation related information for selected points
!                to files in NetCDF format, 1 file per point.
!
! Description:   This module reads a list of points from a file "point2nc.txt"
!                in the working NetCDF directory. It then dumps the 
!                assimilation related information for these points in NetCDF
!                format to files named enkf_III,JJJ.nc, where III and JJJ - i
!                and j grid coordinates.
!
! Modifications: PS 4/8/2010 "point2nc.txt" now allows comments etc. E.g. put
!                            "#" in front of an entry to comment it out.

module m_point2nc
  use m_parameters
  implicit none

  integer, private :: FID = 31
  integer, parameter, private :: STRLEN = 512

  public p2nc_init
  public p2nc_testthiscell
  public p2nc_writeobs
  public p2nc_storeforecast
  public p2nc_writeforecast

  integer, private :: npoints
  integer, allocatable, dimension(:), private :: icoords, jcoords
  real(4), allocatable, dimension(:, :, :) :: forecast

contains

  ! Initialise the point output.
  !
  subroutine p2nc_init()
#if defined (QMPI)
    use qmpi
#else
    use qmpi_fake
#endif

    character(STRLEN) :: line
    integer :: iostatus
    integer :: i, j, n

    npoints = 0

    open(FID, file = trim(POINTFNAME), status = 'old', iostat = iostatus)
    if (iostatus /= 0) then
       if (master) then
          print *, 'WARNING: could not open "', trim(POINTFNAME), '" for reading'
          print *, '         no point output will be performed'
       end if
       return
    end if
    
    do while (.true.)
       read(FID, '(a)', iostat = iostatus) line
       if (iostatus == 0) then
          read(line, *, iostat = iostatus) i, j
          if (iostatus == 0) then
             npoints = npoints + 1
          end if
       else
          exit
       end if
    end do
    close(FID)

    if (master) then
       print '(a, i3, a)', ' p2nc: ', npoints, ' points specified'
    end if

    allocate(icoords(npoints), jcoords(npoints))

    open(FID, file = trim(POINTFNAME), status = 'old', iostat = iostatus)
    if (iostatus /= 0) then
       print *, 'ERROR: point2nc: I/O problem'
       stop
    end if
    
    n = 0
    do while (n < npoints)
       read(FID, '(a)', iostat = iostatus) line
       if (iostatus == 0) then
          read(line, *, iostat = iostatus) i, j
          if (iostatus == 0) then
             n = n + 1
             icoords(n) = i
             jcoords(n) = j
             if (master) then
                print '(a, i3, a, i4, a, i4)', '   point', n, ': i =', i, ', j =', j
             end if
          end if
       end if
    end do
    close(FID)
    if (master) then
       print *
    end if
  end subroutine p2nc_init

  
  ! Test if the output is requested for the point (i, j) 
  !
  function p2nc_testthiscell(i, j)
    logical :: p2nc_testthiscell
    integer, intent(in) :: i, j

    integer :: p

    p2nc_testthiscell = .false.
    do p = 1, npoints
       if (i == icoords(p) .and. j == jcoords(p)) then
          p2nc_testthiscell = .true.
          return
       end if
    end do
  end function p2nc_testthiscell


  ! Write the assimilation parameters (local observations and the X5 matrices)
  ! to the point output files.
  !
  subroutine p2nc_writeobs(i, j, nlobs, nens, X5, lon, lat, depth, rfactor,&
       ids, lobs, Hx, S, ss, lfactors)
    use mod_measurement
    use m_obs
    use nfw_mod

    integer, intent(in) :: i, j, nlobs, nens
    real, intent(in) :: X5(nens, nens)
    real, intent(in) :: lon(1), lat(1), depth(1)
    real, intent(in), optional :: rfactor(1)
    integer, intent(in), optional :: ids(nlobs)
    type(measurement), intent(in), optional :: lobs(nlobs)
    real, intent(in), optional :: Hx(nlobs)
    real, intent(in), optional :: S(nlobs, nens)
    real, intent(in), optional :: ss(nlobs), lfactors(nlobs)

    character(STRLEN) :: fname
    character(STRLEN) :: typename
    integer :: ncid
    integer :: dids(2)
    integer :: vid_ids, vid_lon, vid_lat, vid_val, vid_var, vid_hx, vid_s, vid_x5
    integer :: vid_a1, vid_a2, vid_a3, vid_a4, vid_otype, vid_ss, vid_lfactors
    integer :: otype(nlobs)
    integer :: o, ot
    
    write(fname, '(a, i0.3, ",", i0.3, ".nc")') 'enkf_', i, j
    call nfw_create(fname, nf_write, ncid)

    call nfw_def_dim(fname, ncid, 'p', nlobs, dids(2))
    call nfw_def_dim(fname, ncid, 'm', nens, dids(1))
    if (nlobs > 0) then
       call nfw_def_var(fname, ncid, 'obs_ids', nf_int, 1, dids(2), vid_ids)
       call nfw_def_var(fname, ncid, 'Hx', nf_double, 1, dids(2), vid_hx)
       call nfw_def_var(fname, ncid, 'lon', nf_double, 1, dids(2), vid_lon)
       call nfw_def_var(fname, ncid, 'lat', nf_double, 1, dids(2), vid_lat)
       call nfw_def_var(fname, ncid, 'obs_val', nf_double, 1, dids(2), vid_val)
       call nfw_def_var(fname, ncid, 'obs_var', nf_double, 1, dids(2), vid_var)
       call nfw_def_var(fname, ncid, 'a1', nf_double, 1, dids(2), vid_a1)
       call nfw_def_var(fname, ncid, 'a2', nf_double, 1, dids(2), vid_a2)
       call nfw_def_var(fname, ncid, 'a3', nf_double, 1, dids(2), vid_a3)
       call nfw_def_var(fname, ncid, 'a4', nf_double, 1, dids(2), vid_a4)
       call nfw_def_var(fname, ncid, 'obs_type', nf_int, 1, dids(2), vid_otype)
       call nfw_def_var(fname, ncid, 'S', nf_double, 2, dids, vid_s)
       call nfw_def_var(fname, ncid, 's', nf_double, 1, dids(2), vid_ss)
       call nfw_def_var(fname, ncid, 'lfactors', nf_double, 1, dids(2), vid_lfactors)
    end if
    dids(2) = dids(1)
    call nfw_def_var(fname, ncid, 'X5', nf_double, 2, dids, vid_x5)

    call nfw_put_att_double(fname, ncid, nf_global, 'lon', nf_double, 1, lon)
    call nfw_put_att_double(fname, ncid, nf_global, 'lat', nf_double, 1, lat)
    call nfw_put_att_double(fname, ncid, nf_global, 'depth', nf_double, 1, depth)
    call nfw_put_att_double(fname, ncid, nf_global, 'rfactor', nf_double, 1, rfactor)

    do ot = 1, nuobs
       write(typename, '(a, i1)') 'obstype', ot
       call nfw_put_att_text(fname, ncid, nf_global, typename, len_trim(unique_obs(ot)), trim(unique_obs(ot)))
    end do

    call nfw_enddef(fname, ncid)

    if (nlobs > 0) then
       call nfw_put_var_double(fname, ncid, vid_hx, Hx)
       call nfw_put_var_int(fname, ncid, vid_ids, ids)
       call nfw_put_var_double(fname, ncid, vid_lon, lobs % lon)
       call nfw_put_var_double(fname, ncid, vid_lat, lobs % lat)
       call nfw_put_var_double(fname, ncid, vid_val, lobs % d)
       call nfw_put_var_double(fname, ncid, vid_var, lobs % var)
       call nfw_put_var_double(fname, ncid, vid_a1, lobs % a1)
       call nfw_put_var_double(fname, ncid, vid_a2, lobs % a2)
       call nfw_put_var_double(fname, ncid, vid_a3, lobs % a3)
       call nfw_put_var_double(fname, ncid, vid_a4, lobs % a4)
       otype = 0
       do o = 1, nlobs
          do ot = 1, nuobs
             if (trim(lobs(o) % id) == trim(unique_obs(ot))) then
                otype(o) = ot
             end if
          end do
       end do

       call nfw_put_var_int(fname, ncid, vid_otype, otype)
       call nfw_put_var_double(fname, ncid, vid_s, transpose(S))
       call nfw_put_var_double(fname, ncid, vid_ss, ss)
       call nfw_put_var_double(fname, ncid, vid_lfactors, lfactors)
    end if

    call nfw_put_var_double(fname, ncid, vid_x5, transpose(X5))

    call nfw_close(fname, ncid)
  end subroutine p2nc_writeobs


  ! Store the values of the forecast field No. `fid' in each output point to
  ! the variable `forecast'.
  !
  subroutine p2nc_storeforecast(ni, nj, nrens, nfields, fid, field)
    integer, intent(in) :: ni, nj ! size of grid
    integer, intent(in) :: nrens
    integer, intent(in) :: nfields
    integer, intent(in) :: fid
    real(4), dimension(ni * nj, nrens), intent(in) :: field

    integer :: n

    if (npoints == 0) then
       return
    end if

    if (.not. allocated(forecast)) then
       allocate(forecast(nrens, npoints, nfields))
    end if

    do n = 1, npoints
       forecast(:, n, fid) = field((jcoords(n) - 1) * ni + icoords(n), :)
    end do
  end subroutine p2nc_storeforecast


  ! This procedure consolidates all forecast fields for each output point 
  ! together in the variable `forecast' on the master node, and then writes
  ! them to the appropriate NetCDF files.
  !
  subroutine p2nc_writeforecast
#if defined (QMPI)
    use qmpi
#else
    use qmpi_fake
#endif
    use distribute
    use nfw_mod
    use mod_analysisfields
    implicit none
    
    character(STRLEN) :: fname
    integer :: p, k, nf
    character(8) :: varname
    integer kstart
    integer ncid, dids(2), varid, nf2

#if defined(QMPI)
    if (.not. master) then
       call send(forecast(:, :, my_first_iteration : my_last_iteration), 0, 0)
       return ! leave writing to master
    else
       do p = 2, qmpi_num_proc ! here p is the MPI node ID
          call receive(forecast(:, :, first_iteration(p) : last_iteration(p)), p - 1, 0)
       end do
    end if
#endif

    ! only master keeps working here
    !
    do p = 1, npoints
       write(fname, '(a, i0.3, ",", i0.3, ".nc")') 'enkf_', icoords(p), jcoords(p)
       call nfw_open(fname, nf_write, ncid)
       call nfw_redef(fname, ncid)
       call nfw_inq_dimid(fname, ncid, 'm', dids(1))
       call nfw_enddef(fname, ncid)
    
       kstart = -1
       do k = 1, numfields
          if (kstart == -1) then
             kstart = k
             varname = fieldnames(k)
          end if

          ! check if there are more fields for this variable
          !
          if (k < numfields .and. fieldnames(k + 1) == varname) then
             cycle
          end if

          ! this is the last field for this variable - write the variable
          !
          nf = k - kstart + 1

          call nfw_redef(fname, ncid)

          if (nf == 1) then
             call nfw_def_var(fname, ncid, trim(varname), nf_float, 1, dids(1), varid)
          else
             if (.not. nfw_dim_exists(ncid, 'k')) then
                call nfw_def_dim(fname, ncid, 'k', nf, dids(2))
             else
                call nfw_inq_dimid(fname, ncid, 'k', dids(2))
                call nfw_inq_dimlen(fname, ncid, dids(2), nf2)
                if (nf /= nf2) then
                   print *, 'ERROR: p2nc_writeforecast(): varname = "', trim(varname),&
                        '", # levels = ', nf, '# levels in "', trim(fname), '" =', nf2
                   print *, 'ERROR: p2nc_writeforecast(): returning'
                end if
             end if
             call nfw_def_var(fname, ncid, trim(varname), nf_float, 2, dids, varid)
          end if

          call nfw_enddef(fname, ncid)

          call nfw_put_var_real(fname, ncid, varid, forecast(:, p, kstart : kstart + nf - 1))

          kstart = -1
       end do
       call nfw_close(fname, ncid)
    end do
  end subroutine p2nc_writeforecast

end module m_point2nc
