module m_get_mod_nrens

contains

  integer function get_mod_nrens(nx, ny)
#if defined (QMPI)
    use qmpi , only : stop_mpi, master
#else
    use qmpi_fake
#endif
    implicit none

    integer, intent(in) :: nx, ny

    integer j
    logical ex
    logical exa, exb
    character(len=3) :: cmem
    integer get_mod_ice, imem
    character(len=*), parameter :: icefile = 'forecastICE.uf'
    real*8, dimension(nx,ny) :: ficem, hicem, hsnwm, ticem, tsrfm

    !KAL predefined names memXXX.[ab] -- here we just check presence
    imem = 1
    exa = .true.
    exb = .true.
    do while (exa .and. exb)
       write(cmem, '(i3.3)') imem
       inquire(exist = exa, file = 'forecast' // cmem // '.a')
       inquire(exist = exb, file = 'forecast' // cmem // '.b')
       if (exa .and. exb) then
          imem = imem + 1
       end if
    end do
    get_mod_nrens = imem - 1

#if ! defined (SINGLE_RESTART)
#if defined (ICE)
    inquire(exist = ex, file = icefile)
    if (.not. ex) then
       if (master) then
          print *, 'ERROR: ' // icefile // ' does not exist!'
       end if
       call stop_mpi()
    end if

    inquire(iolength = j) ficem, hicem, hsnwm, ticem, tsrfm
    open(10, file = icefile, status = 'old', form = 'unformatted', access = 'direct', recl = j)
    do j = 1, 1000
       read(10, rec = j, err = 200) ficem, hicem, hsnwm, ticem, tsrfm
    enddo
200 close(10)
    get_mod_ice = j - 1
    if (master) then
       print '(a, i4, a)', ' EnKF: ',  get_mod_nrens, ' ocean states found'
       print '(a, i4, a)', ' EnKF: ',  get_mod_ice, ' ice  states found'
    end if
    get_mod_nrens = min(get_mod_nrens, get_mod_ice)
#endif
#endif
  end function get_mod_nrens

end module m_get_mod_nrens
