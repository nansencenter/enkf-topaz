module distribute

#if defined(QMPI)
  use qmpi
#else
  use qmpi_fake
#endif

  !
  ! public stuff
  !
  integer, public :: my_number_of_iterations, my_first_iteration, my_last_iteration
  integer, dimension(:), allocatable, public :: number_of_iterations, first_iteration, last_iteration
  integer, dimension(:), allocatable, public :: randommap

contains

  subroutine distribute_iterations(nz)
    implicit none

    integer, intent(in) :: nz

    integer :: i, j
    real(8) :: num_procs_real, mean_iterations

    if (.not. allocated(number_of_iterations)) then
       allocate(number_of_iterations(qmpi_num_proc))
    end if
    if (.not. allocated(first_iteration)) then
       allocate(first_iteration(qmpi_num_proc))
    end if
    if (.not. allocated(last_iteration)) then
       allocate(last_iteration(qmpi_num_proc))
    end if

    if (master) then
       print *, 'Distribution of iterations:'
    end if

    num_procs_real = qmpi_num_proc
    mean_iterations = nz / num_procs_real

    j = -1
    if (int(mean_iterations) .eq. mean_iterations) then
       my_number_of_iterations = nz/qmpi_num_proc
       if (master) then
          number_of_iterations(:) = nz / qmpi_num_proc
          print *, 'All procs get ', number_of_iterations(1), 'iterations'
       endif
       j = qmpi_num_proc
    else
       do i = 1, qmpi_num_proc
          if (i * floor(mean_iterations) +&
               (qmpi_num_proc-i) * ceiling(mean_iterations) .eq. nz) then
             j = i
             exit
          endif
       end do

       if (qmpi_proc_num + 1 .le. j) then
          my_number_of_iterations = floor(mean_iterations)
       else
          my_number_of_iterations = ceiling(mean_iterations)
       endif

       if (master) then
          number_of_iterations(1:j) = floor(mean_iterations)
          number_of_iterations(j+1:qmpi_num_proc) = ceiling(mean_iterations)
          if ((j * floor(mean_iterations) +&
               (qmpi_num_proc - j) * ceiling(mean_iterations)) .ne. nz) then
             print *, 'ERROR in distribute_iteration()'
             stop
          endif
          if (nz .lt. qmpi_num_proc) then
             print *, 'Number of cells in z-direction than number of processors'
             stop
          endif
       endif
    endif

    if (master) then
       first_iteration(1) = 1; 
       last_iteration(1) = number_of_iterations(1)
       do i = 2, qmpi_num_proc
          first_iteration(i) = last_iteration(i - 1) + 1 
          last_iteration(i) = first_iteration(i) + number_of_iterations(i)-1
       end do
    endif

    if (qmpi_proc_num + 1 .le. j) then
       my_first_iteration = qmpi_proc_num*my_number_of_iterations + 1
    else
       my_first_iteration = j * (my_number_of_iterations - 1) +&
            (qmpi_proc_num - j) * my_number_of_iterations + 1
    endif
    my_last_iteration = my_first_iteration + my_number_of_iterations - 1

    print *, 'I am', qmpi_proc_num, ', my_first_ind =', my_first_iteration,&
         ', my_last_ind =', my_last_iteration
  end subroutine distribute_iterations

end module distribute

