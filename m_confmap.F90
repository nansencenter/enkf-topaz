module m_confmap
  implicit none

  logical :: confmap_initialised = .false.

  real :: pi_1
  real :: pi_2
  real :: deg
  real :: rad
  real :: theta_a
  real :: phi_a
  real :: theta_b
  real :: phi_b
  real :: di
  real :: dj
  complex :: imagone
  complex :: ac
  complex :: bc
  complex :: cmna
  complex :: cmnb
  real :: mu_s
  real :: psi_s
  real :: epsil
  logical :: mercator

  real :: lat_a, lon_a
  real :: lat_b, lon_b
  real :: wlim, elim
  real :: slim, nlim
  real :: mercfac
  integer :: ires, jres

contains

  ! This routine initializes constants used in the conformal mapping
  ! and must be called before the routines 'oldtonew' and 'newtoold'
  ! are called. The arguments of this routine are the locations of
  ! the two poles in the old coordiante system.
  !
  subroutine confmap_init(nx, ny)
    integer, intent(in) :: nx, ny

    real :: cx, cy, cz, theta_c, phi_c
    complex :: c, w
    logical :: ass, lold
    
    ! Read info file
    open(unit = 10, file = 'grid.info', form = 'formatted')
    read(10, *) lat_a, lon_a
    read(10, *) lat_b,lon_b
    read(10, *) wlim, elim, ires
    read(10, *) slim, nlim, jres
    read(10, *) ass
    read(10, *) ass
    read(10, *) ass
    read(10, *) mercator
    read(10, *) mercfac, lold
    close(10)
    if (ires /= nx .and. jres /= ny) then
       print *, 'initconfmap: WARNING -- the dimensions in grid.info are not'
       print *, 'initconfmap: WARNING -- consistent with nx and ny'
       print *, 'initconfmap: WARNING -- IGNORE IF RUNNING CURVIINT'
       stop '(initconfmap)'
    endif

    ! some constants
    !
    pi_1 = 3.14159265358979323846
    pi_2 = 0.5 * pi_1
    deg = 180.0 / pi_1
    rad = 1.0 / deg
    epsil = 1.0d-9

    di = (elim - wlim) / real(ires - 1)   ! delta lon'
    dj = (nlim - slim) / real(jres - 1)   ! delta lat' for spherical grid

    if (mercator) then
       dj = di
       if (lold) then
          print *, 'initconfmap: lold'
          slim = -mercfac * jres * dj
       else
          print *, 'initconfmap: not lold'
          slim = mercfac
       endif
    endif

    ! transform to spherical coordinates
    !
    theta_a = lon_a * rad
    phi_a = pi_2 - lat_a * rad
    theta_b = lon_b * rad
    phi_b = pi_2 - lat_b * rad

    ! find the angles of a vector pointing at a point located exactly
    ! between the poles
    !
    cx = cos(theta_a) * sin(phi_a) + cos(theta_b) * sin(phi_b)
    cy = sin(theta_a) * sin(phi_a) + sin(theta_b) * sin(phi_b)
    cz = cos(phi_a) + cos(phi_b)

    theta_c = atan2(cy, cx)
    phi_c = pi_2 - atan2(cz, sqrt(cx * cx + cy * cy))

    ! initialize constants used in the conformal mapping
    !
    imagone = (0.0, 1.0)
    ac = tan(0.5 * phi_a) * exp(imagone * theta_a)
    bc = tan(0.5 * phi_b) * exp(imagone * theta_b)
    c = tan(0.5 * phi_c) * exp(imagone * theta_c)
    cmna = c - ac
    cmnb = c - bc

    w = cmnb / cmna
    mu_s = atan2(aimag(w), real(w))
    psi_s = 2.0 * atan(abs(w))

    confmap_initialised = .true.
  end subroutine confmap_init

end module m_confmap
