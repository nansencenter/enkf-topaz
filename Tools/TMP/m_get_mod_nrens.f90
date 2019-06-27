










module m_get_mod_nrens

contains

  integer function get_mod_nrens(nx, ny)
    use qmpi , only : stop_mpi, master
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

  end function get_mod_nrens

end module m_get_mod_nrens
