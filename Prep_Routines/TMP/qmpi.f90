












module qmpi_fake
  implicit none

  logical, parameter :: master = .true.
  integer, parameter :: qmpi_num_proc = 1
  integer, parameter :: qmpi_proc_num = 0

contains

  subroutine stop_mpi()
    stop 
  end subroutine stop_mpi

end module qmpi_fake

