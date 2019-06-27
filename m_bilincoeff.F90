module m_bilincoeff
  use m_oldtonew
  implicit none

contains

  ! This subroutine uses bilinear interpolation to interpolate the field
  ! computed by the model (MICOM) to the position defined by lon, lat
  ! The output is the interpolation coeffisients a[1-4]
  ! NB  NO locations on land.
  !
  subroutine bilincoeff(glon, glat, nx, ny, lon, lat, ipiv, jpiv, a1, a2, a3, a4)
    real, intent(in) :: glon(nx, ny), glat(nx, ny)
    integer, intent(in) :: nx ,ny
    real, intent(in) :: lon, lat
    integer, intent(in) :: ipiv, jpiv
    real, intent(out) :: a1, a2, a3, a4

    real :: t, u
    real :: lat1, lon1, lat2, lon2, latn, lonn


    call oldtonew(glat(ipiv, jpiv), glon(ipiv, jpiv), lat1, lon1)
    call oldtonew(glat(ipiv + 1, jpiv + 1), glon(ipiv + 1, jpiv + 1), lat2, lon2)
    call oldtonew(lat, lon, latn, lonn)

    t = (lonn - lon1) / (lon2 - lon1)
    u = (latn - lat1) / (lat2 - lat1)

    if (t < -0.1 .or. t > 1.1 .or. u < -0.1 .or. u > 1.1) then
       print *, 'ERROR: bilincoeff(): t, u = ', t, u, 'for lon, lat =', lon, lat
       stop
    end if

    a1 = (1.0 - t) * (1.0 - u)
    a2 = t * (1.0 - u)
    a3 = t * u
    a4 = (1.0 - t) * u
  end subroutine bilincoeff

  subroutine bilincoeff1(glon, glat, nx, ny, lon, lat, ipiv, jpiv, a1, a2, a3, a4)
    real, intent(in) :: glon(nx, ny), glat(nx, ny)
    integer, intent(in) :: nx ,ny
    real, intent(in) :: lon, lat
    integer, intent(in) :: ipiv, jpiv
    real, intent(out) :: a1, a2, a3, a4

    real :: xx(4), yy(4)
    real :: t, u

    xx(1) = glon(ipiv, jpiv)
    xx(2) = glon(ipiv + 1, jpiv)
    xx(3) = glon(ipiv + 1, jpiv + 1)
    xx(4) = glon(ipiv, jpiv + 1)
    yy(1) = glat(ipiv, jpiv)
    yy(2) = glat(ipiv + 1, jpiv)
    yy(3) = glat(ipiv + 1, jpiv + 1)
    yy(4) = glat(ipiv, jpiv + 1)
    call xy2fij(lon, lat, xx, yy, t, u)
    if (t < 0 .or. t > 1 .or. u < 0 .or. u > 1) then
       print *, 'ERROR: bilincoeff(): t, u = ', t, u, 'for lon, lat =', lon, lat
       !       stop
    end if

    a1 = (1.0 - t) * (1.0 - u)
    a2 = t * (1.0 - u)
    a3 = t * u
    a4 = (1.0 - t) * u
  end subroutine bilincoeff1

  subroutine xy2fij(x, y, xx, yy, fi, fj)
    real, intent(in) :: x, y
    real, intent(in) :: xx(4), yy(4)
    real, intent(out) :: fi, fj

    real :: a, b, c, d, e, f, g, h
    real :: aa, bb, cc
    real :: d1, d2

    a = xx(1) - xx(2) - xx(4) + xx(3)
    b = xx(2) - xx(1)
    c = xx(4) - xx(1)
    d = xx(1)
    e = yy(1) - yy(2) - yy(4) + yy(3)
    f = yy(2) - yy(1)
    g = yy(4) - yy(1)
    h = yy(1)

    aa = a * f - b * e;
    bb = e * x - a * y + a * h - d * e + c * f - b * g;
    cc = g * x - c * y + c * h - d * g;

    if (abs(aa) < 1d-5) then
       fi = -cc / bb * (1.0d0 + aa * cc / bb / bb);
    else
       fi = (-bb - sqrt(bb * bb - 4.0d0 * aa * cc)) / (2.0d0 * aa);
    end if
    d1 = a * fi + c
    d2 = e * fi + g
    if (abs(d2) > abs(d1)) then
       fj = (y - f * fi - h) / d2
    else
       fj = (x - b * fi - d) / d1
    end if
  end subroutine xy2fij

end module m_bilincoeff
