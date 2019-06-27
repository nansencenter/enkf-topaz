










! Bilinear coeffisients are calculated witin this program for all the 
! observation points.
! Only wet points are stored in the observation.uf file

module m_get_def_wet_point

  implicit none

  integer, parameter, private :: STRLEN = 512
  character(STRLEN), parameter, private :: MEANSSHFNAME = "meanssh.uf"
  
  private read_mean_ssh
  private land_nearby

contains 

  subroutine get_def_wet_point(obs, data, gr, depths, modlat, modlon, nrobs, nx, ny)
    ! Program converts to a general format readable by ENKF (observations.uf)
    use mod_measurement
    use mod_grid
    ! Functions to be used
    use m_oldtonew
    use m_bilincoeff
    use m_pivotp
    use m_confmap
    use m_spherdist

    integer, intent(in) :: nx, ny
    type (measurement), intent(in) :: data(:)
    type (measurement), intent(inout)   :: obs(:)
    type (grid), intent(in) :: gr      ! observations grid
    real, dimension(nx, ny), intent(in)  ::  depths, modlat, modlon
    integer, intent(out) :: nrobs
    integer, parameter :: maxobs = 1441 * 760 !2*400*600 ! maximum number of observations

    real, dimension(nx, ny) :: mean_ssh
    integer k, imin, imax, jmin, jmax
    integer ipiv, jpiv, nsupport, nsmin, nsmax
    real :: x0, y0
    real wetsill, griddiag, mingridsize, minobssize

    logical wet

    ! gr = default_grid
    nrobs = 0; 
    nsmin = maxobs; 
    nsmax = 0
    mingridsize = 1.E+10; 
    minobssize = 1.E+10    ! in meters

    call confmap_init(nx,ny) ! Initialize conformal mapping before calling

    !Calculate pivot points
    !Find wet points (all neigbours in water)
    !Find the points with defined data value
    !Put the data into the obs data structture
    !Compute bilinear coefficients
    print *, 'test: ',gridpoints(gr), ' !!'
    call read_mean_ssh(mean_ssh, nx, ny)

    do k = 1, gridpoints(gr)
       if (data(k) % id .eq. 'SLA' .or. data(k) % id .eq. 'sla' .or. &
            data(k) % id.eq. 'SSH' .or. data(k)%id .eq. 'ssh' .or.&
            data(k)%id.eq.'TSLA') then
          wetsill = 200.   ! Discarding data in shallow waters
       else
          wetsill=10.
       endif
       call oldtonew(data(k) % lat, data(k) % lon, y0, x0)
       call pivotp(x0, y0, ipiv, jpiv)
       ! Discard obs on model boundaries (TODO: cyclic domains) 
       ! Also valid if ns=0 
       imin = ipiv - data(k) % ns
       imax = ipiv + data(k) % ns + 1 
       jmin = jpiv - data(k) % ns
       jmax = jpiv + data(k) % ns + 1
       if ((imin .le. 0) .or. (jmin .le. 0) .or. (imax .ge. nx) .or. &
            (jmax .ge. ny)) cycle
       ! Is observation surrounded by wet grid points?
       if (any(depths(imin:imax, jmin:jmax) < wetsill .or. depths(imin:imax, jmin:jmax) == depths(imin:imax, jmin:jmax) + 1.0)) cycle
       wet = data(k) % status ! Discards inconsistent/Fill values
       if (data(k) % id .eq. 'SLA' .or. data(k) % id .eq. 'sla' .or.&
            data(k) % id .eq. 'TSLA') then
          wet = wet .and. .not. land_nearby(nx, ny, depths, modlon, modlat,&
               ipiv, jpiv, data(k) % lon, data(k) % lat,50000.0)
       endif
       if (data(k) % id .eq. 'HICE') then
          wet = wet .and. .not. land_nearby(nx, ny, depths, modlon, modlat,&
               ipiv, jpiv, data(k) % lon, data(k) % lat,30000.0)
       endif

       if(.not. undefined(data(k) % d, gr) .and. wet) then
          nrobs = nrobs + 1
          obs(nrobs) = data(k)
          obs(nrobs) % ipiv = ipiv
          obs(nrobs) % jpiv=  jpiv
          obs(nrobs) % status = .true. ! Wet obs
          if (data(k) % ns > 0) then   ! large support data: a1 is the obs support(m)
             griddiag = spherdist(modlon(ipiv, jpiv), modlat(ipiv, jpiv), &
                  modlon(ipiv + 1, jpiv + 1), modlat(ipiv + 1, jpiv + 1))
             !FC: 0.5 because m_Generate_element_Si runs from -ns to +ns 
             nsupport = floor(0.5 * data(k) % a1 / griddiag) ! both in meters
             obs(nrobs)%ns = nsupport  ! number of grid cells in the diagonal
             nsmin = min(nsmin, nsupport)
             nsmax = max(nsmax, nsupport)
             mingridsize = min(mingridsize, griddiag)
             minobssize = min(minobssize, data(k) % a1)
          else
             obs(nrobs) % ns = 0    ! point measurements have zero support
          endif
          call bilincoeff(modlon, modlat, nx, ny, obs(nrobs)%lon, &
               obs(nrobs) % lat, obs(nrobs) % ipiv, obs(nrobs) % jpiv, &
               obs(nrobs) % a1, obs(nrobs) % a2, obs(nrobs) % a3, &
               obs(nrobs) % a4)
       endif
    end do
    print*, 'Number of defined and wet observations: nrobs ', nrobs
    print*, 'Support (in nb of cells) between: ', nsmin, ' and ', nsmax
    print '(2(a,f8.3),a)', ' Minimum obs support: ', 0.001 * minobssize, &
         'km, min grid diagonal: ', 0.001 * mingridsize, ' km' 
  end subroutine get_def_wet_point


  subroutine read_mean_ssh(mean_ssh, nx, ny)
    integer, intent(in) :: nx, ny
    real, intent(out):: mean_ssh(nx, ny)
    logical :: exists

    inquire(file = trim(MEANSSHFNAME), exist = exists)
    if (.not. exists) then
       print *,'ERROR: read_mean_ssh(): file "', trim(MEANSSHFNAME), '" not found'
       stop
    end if
       
    open (10, file = trim(MEANSSHFNAME), status = 'unknown',form = 'unformatted')
    read (10) mean_ssh
    close (10)
  end subroutine read_mean_ssh


  logical function land_nearby(nx, ny, depths, modlon, modlat, ipiv, jpiv, obslon, obslat,Vdis)
    use m_spherdist
    implicit none
    integer, intent (in) :: nx, ny, ipiv, jpiv
    real, dimension(nx,ny), intent(in) :: depths, modlon, modlat
    real, intent (in) :: obslon,obslat 
    real, intent (in) :: Vdis    ! ~ meters minimum far from coastline

    integer :: ii, jj, ncells
    real :: griddist

    land_nearby = .false.
    ncells = ceiling(Vdis / spherdist(modlon(ipiv, jpiv), modlat(ipiv, jpiv),&
         modlon(ipiv, jpiv + 1), modlat(ipiv, jpiv + 1)))
    do jj = max(jpiv - ncells, 1), min(jpiv + ncells, ny)
       do ii = max(ipiv - ncells, 1), min(ipiv + ncells, nx)
          griddist = spherdist(modlon(ii, jj), modlat(ii, jj), obslon, obslat)
          if ((depths(ii,jj) < 1.or.depths(ii,jj)>20000) .and. griddist < Vdis) then
             print *, 'land_nearby: ',ipiv,jpiv
          !   print *,modlon(ii, jj), modlat(ii, jj), obslon, obslat
          !   print *,griddist,depths(ii,jj)
             land_nearby = .true.
             return
          end if
       enddo
    enddo
  end function land_nearby

end module m_get_def_wet_point
