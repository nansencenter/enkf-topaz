#if defined(QMPI)
module qmpi
!
! A module defining a minimalist interface to a subset of MPI.
! The first five primitives can in theory be used to parallelize
! any program. The module hides type specification, communicators,
! explicit error handling, the need to give explicit buffer size etc.
! Also provided are a few interfaces for often used broadcast and 
! reduction operations
!
! © Helge Avlesen <avle@ii.uib.no>, para//ab
!
! primitives: (optional arguments in brackets)
!
!   subroutine start_mpi()
!      starts the mpi subsystem. all processesors are assigned a number (myid).
!      the number of processors is numproc.
!   subroutine stop_mpi()
!      stops the mpi subsystem
!   subroutine barrier([label])
!      syncronization point for all processors. optionally prints a label on
!      the master processor (0).
!   subroutine send(data, target [,tag])
!      send object data to processor number target, tag is an optional integer
!      that defaults to 0. (if multiple messages are exchanged between a
!      pair of processors, a unique tag must be used for each exhange)
!   subroutine receive(data, source [,tag])
!      get object data from processor source, tag is optional and as for send
!      MPI will fail if the size of the object received is different from what
!      was sent.
!  
! The rest of the routines are included for convenience, they can be
! also be implemented using the above subroutines.
!
!   subroutine broadcast(data [,root])
!      broadcast data (any type) from processor root (default=0) to all
!      other processors.
!   subroutine mbroadcast(data [,data2,data3,data4,data5,data6] [,root])
!      broadcast up to 6 scalar variables of the same type, to all processors
!      from processor root (default=0)
!   subroutine reduce(type, data [,data2,data3,data4,data5,data6] [,root] )
!      reduce the scalar data, optionally also data2-data6, return result
!      on all processes. the operation can currently be of type 'sum', 'mul',
!      'min' or 'max' i.e. a sum or a product. data-data6 must be of the 
!      same type. if integer root is present, only return result on that 
!      processor (faster)
!
! Example: a program that sends a real from processor 0 to processor 1
!   use qmpi
!   real data
!   call start_mpi
!   data=myid
!   if(myid==0) call send(data, 1)
!   if(myid==1) then
!      call receive(data, 0)
!      print *,'hello, I am',myid,'got ',data,'from process 0'
!   end if
!   call stop_mpi
!   end
!
! More advanced usage example: to send a derived type from 0 to 1; 
! pack it in a string (could be packed into any array), send, receive, unpack.
! 
! type(any_type) var1
! character, allocatable :: buffer(:)
! ...
! N=size(transfer(var1,(/'x'/))))   !! compute size of type once
! allocate(buffer(N))
! if(myid==0)then
!     buffer = transfer(var1,buffer)
!     call send(buffer,1)
! end if
! if(myid==1)then
!     call receive(buffer,0)
!     var1 = transfer(buffer,var1)
! end if
! ...
!  
#warning "COMPILING WITH QMPI CODE"
  include 'mpif.h'
  integer, public :: qmpi_proc_num, qmpi_num_proc, ierr, errorcode, mpistatus(mpi_status_size)
  logical, public :: master=.false., slave=.false.

! some kinds. could use selected_real_kind(..) for this instead of hard coding
  integer, parameter :: dp=8, sp=4, long=8, short=2

  interface send
     module procedure            &
          qmpi_send_real4,       &
          qmpi_send_real4_1d,    &
          qmpi_send_real4_2d,    &
          qmpi_send_real4_3d,    &
          qmpi_send_real4_4d,    &
          qmpi_send_real8,       &
          qmpi_send_real8_1d,    &
          qmpi_send_real8_2d,    &
          qmpi_send_real8_3d,    &
          qmpi_send_real8_4d,    &
          qmpi_send_integer4,    &
          qmpi_send_integer4_1d, &
          qmpi_send_integer4_2d, &
          qmpi_send_integer4_3d, &
          qmpi_send_integer4_4d, &
          qmpi_send_integer8,    &
          qmpi_send_integer8_1d, &
          qmpi_send_integer8_2d, &
          qmpi_send_integer8_3d, &
          qmpi_send_integer8_4d, &
          qmpi_send_string,      &
          qmpi_send_character_1d,&
          qmpi_send_logical
  end interface

  interface receive
     module procedure &
          qmpi_recv_real4,       &
          qmpi_recv_real4_1d,    &
          qmpi_recv_real4_2d,    &
          qmpi_recv_real4_3d,    &
          qmpi_recv_real4_4d,    &
          qmpi_recv_real8,       &
          qmpi_recv_real8_1d,    &
          qmpi_recv_real8_2d,    &
          qmpi_recv_real8_3d,    &
          qmpi_recv_real8_4d,    &
          qmpi_recv_integer4,    & 
          qmpi_recv_integer4_1d, &
          qmpi_recv_integer4_2d, &
          qmpi_recv_integer4_3d, &
          qmpi_recv_integer4_4d, &
          qmpi_recv_integer8,    &
          qmpi_recv_integer8_1d, &
          qmpi_recv_integer8_2d, &
          qmpi_recv_integer8_3d, &
          qmpi_recv_integer8_4d, &
          qmpi_recv_string,      &
          qmpi_recv_character_1d,&
          qmpi_recv_logical
  end interface

  interface reduce
     module procedure &
          qmpi_integer_reduction, &
          qmpi_integer8_reduction,&
          qmpi_real_reduction,    &
          qmpi_real8_reduction
  end interface

  interface broadcast
     module procedure &
          qmpi_broadcast_logical,  &
          qmpi_broadcast_string,   &
          qmpi_broadcast_stringarr,&
          qmpi_broadcast_integer4, &
          qmpi_broadcast_integer4_array1d,  &
          qmpi_broadcast_integer4_array2d,  &
          qmpi_broadcast_integer8, &
          qmpi_broadcast_integer8_array1d, &
          qmpi_broadcast_integer8_array2d, &
          qmpi_broadcast_real4, &
          qmpi_broadcast_real4_array1d, &
          qmpi_broadcast_real4_array2d, &
          qmpi_broadcast_real4_array3d, &
          qmpi_broadcast_real4_array4d, &
          qmpi_broadcast_real8, &
          qmpi_broadcast_real8_array1d, &
          qmpi_broadcast_real8_array2d, &
          qmpi_broadcast_real8_array3d, &
          qmpi_broadcast_real8_array4d
  end interface

  interface mbroadcast
     module procedure &
          qmpi_broadcast_logicals, &
          qmpi_broadcast_real4s, &
          qmpi_broadcast_real8s, &
          qmpi_broadcast_integer4s, &
          qmpi_broadcast_integer8s
  end interface

contains

  subroutine start_mpi()
!
! initialize the core MPI subsystem
! this routine should be called as the first statement in the program.
! MPI does not specify what happen before MPI_init and after mpi_finalize
!
    implicit none

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world, qmpi_num_proc, ierr)
    call mpi_comm_rank(mpi_comm_world, qmpi_proc_num, ierr)

    master=.false.
    if(qmpi_proc_num==0) master=.true.
    if(qmpi_proc_num>0) slave=.true.
print*,'Inne i start_mpi: qmpi_proc_num =',qmpi_proc_num,' master =',master

    if(master) then
        write(*,'(a,i0,a)') 'MPI started with ',qmpi_num_proc,' processors'
    end if
  end subroutine start_mpi

  subroutine stop_mpi()
    implicit none
    call mpi_finalize(ierr)
    stop
  end subroutine stop_mpi

  subroutine barrier(label)
! makes all processes sync at this point, optionally print a label
    implicit none
    character(*), optional :: label
    call mpi_barrier(mpi_comm_world, ierr)
    if(master.and.present(label)) print *,'---barrier---',label,'---------'
  end subroutine barrier

  subroutine qmpi_send_logical(data, target, tag)
    implicit none
    logical data
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_send(data, counter, mpi_logical, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_logical count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_logical

  subroutine qmpi_send_string(data, target, tag)
    implicit none
    character(*) data
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=len(data)
    call mpi_send(data, counter, mpi_character, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_string count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_string

  subroutine qmpi_send_character_1d(data, target, tag)
    implicit none
    character data(:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_send(data, counter, mpi_character, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_character_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_character_1d
  
  subroutine qmpi_recv_character_1d(data, target, tag)
    implicit none
    character data(:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_recv(data, counter, mpi_character, target, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_character_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_character_1d
    
  subroutine qmpi_send_integer4(data, target, tag)
    implicit none
    integer(sp) data
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_send(data, counter, mpi_integer, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer4 count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer4

  subroutine qmpi_send_integer4_1d(data, target, tag)
    implicit none
    integer(sp) data(:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_send(data, counter, mpi_integer, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer4_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer4_1d

  subroutine qmpi_send_integer4_2d(data, target, tag)
    implicit none
    integer(sp) data(:,:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)
    call mpi_send(data, counter, mpi_integer, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer4_2d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer4_2d

  subroutine qmpi_send_integer4_3d(data, target, tag)
    implicit none
    integer(sp) data(:,:,:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)
    call mpi_send(data, counter, mpi_integer, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer4_3d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer4_3d

  subroutine qmpi_send_integer4_4d(data, target, tag)
    implicit none
    integer(sp) data(:,:,:,:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)
    call mpi_send(data, counter, mpi_integer, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer4_4d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer4_4d

  subroutine qmpi_send_integer8(data, target, tag)
    implicit none
    integer(long) data
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_send(data, counter, mpi_integer8, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer8 count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer8

  subroutine qmpi_send_integer8_1d(data, target, tag)
    implicit none
    integer(long) data(:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_send(data, counter, mpi_integer8, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer8_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer8_1d

  subroutine qmpi_send_integer8_2d(data, target, tag)
    implicit none
    integer(long) data(:,:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)
    call mpi_send(data, counter, mpi_integer8, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer8_2d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer8_2d

  subroutine qmpi_send_integer8_3d(data, target, tag)
    implicit none
    integer(8) data(:,:,:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)
    call mpi_send(data, counter, mpi_integer8, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer8_3d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer8_3d

  subroutine qmpi_send_integer8_4d(data, target, tag)
    implicit none
    integer(8) data(:,:,:,:)
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)
    call mpi_send(data, counter, mpi_integer8, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_integer8_4d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_integer8_4d

  subroutine qmpi_send_real4(data, target, tag)
    implicit none
    real(sp) data
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_send(data, counter, mpi_real, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_real4 count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real4

  subroutine qmpi_send_real8(data, target, tag)
    implicit none
    real(dp) data
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_send(data, counter, mpi_double_precision, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_real8 count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real8

  subroutine qmpi_send_real4_1d(data, target, tag)
    implicit none
    real(sp) data(:)
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_send(data, counter, mpi_real, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_real4_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real4_1d

  subroutine qmpi_send_real8_1d(data, target, tag)
    implicit none
    real(dp) data(:)
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_send(data, counter, mpi_double_precision, target, given_tag, mpi_comm_world, ierr)
    if(ierr.ne.0)then
        print *,'error send_real8_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real8_1d

  subroutine qmpi_send_real4_2d(data, target, tag)
    implicit none
    real(sp) data(:,:)
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)

    call mpi_send(data, counter, mpi_real, target, given_tag, mpi_comm_world, ierr)

    if(ierr.ne.0)then
        print *,'error send_real4_2d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real4_2d

  subroutine qmpi_send_real8_2d(data, target, tag)
    implicit none
    real(dp) data(:,:)
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)

    call mpi_send(data, counter, mpi_double_precision, target, given_tag, mpi_comm_world, ierr)

    if(ierr.ne.0)then
        print *,'error send_real8_2d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real8_2d

  subroutine qmpi_send_real4_3d(data, target, tag)
    implicit none
    real(sp) data(:,:,:)
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)

    call mpi_send(data, counter, mpi_real, target, given_tag, mpi_comm_world, ierr)

    if(ierr.ne.0)then
        print *,'error send_real4_3d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real4_3d

  subroutine qmpi_send_real8_3d(data, target, tag)
    implicit none
    real(dp) data(:,:,:)
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)

    call mpi_send(data, counter, mpi_double_precision, target, given_tag, mpi_comm_world, ierr)

    if(ierr.ne.0)then
        print *,'error send_real8_3d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real8_3d

  subroutine qmpi_send_real4_4d(data, target, tag)
    implicit none
    real(sp) data(:,:,:,:)
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)

    call mpi_send(data, counter, mpi_real, target, given_tag, mpi_comm_world, ierr)

    if(ierr.ne.0)then
        print *,'error send_real4_4d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real4_4d

  subroutine qmpi_send_real8_4d(data, target, tag)
    implicit none
    real(dp) data(:,:,:,:)
    integer target
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)

    call mpi_send(data, counter, mpi_double_precision, target, given_tag, mpi_comm_world, ierr)

    if(ierr.ne.0)then
        print *,'error send_real8_4d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_send_real8_4d

  subroutine qmpi_recv_integer4(data, source, tag)
    implicit none
    integer(sp) data
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_recv(data, counter, mpi_integer, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer4_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer4

  subroutine qmpi_recv_integer4_1d(data, source, tag)
    implicit none
    integer(sp) data(:)
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_recv(data, counter, mpi_integer, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer4_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer4_1d

  subroutine qmpi_recv_integer4_2d(data, source, tag)
    implicit none
    integer(sp) data(:,:)
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)
    call mpi_recv(data, counter, mpi_integer, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer4_2d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer4_2d

  subroutine qmpi_recv_integer4_3d(data, source, tag)
    implicit none
    integer(sp) data(:,:,:)
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)
    call mpi_recv(data, counter, mpi_integer, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer4_3d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer4_3d

  subroutine qmpi_recv_integer4_4d(data, source, tag)
    implicit none
    integer(sp) data(:,:,:,:)
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)
    call mpi_recv(data, counter, mpi_integer, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer4_4d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer4_4d

  subroutine qmpi_recv_integer8(data, source, tag)
    implicit none
    integer(long) data
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_recv(data, counter, mpi_integer8, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer8 count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer8

  subroutine qmpi_recv_integer8_1d(data, source, tag)
    implicit none
    integer(long) data(:)
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_recv(data, counter, mpi_integer8, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer8_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer8_1d

  subroutine qmpi_recv_integer8_2d(data, source, tag)
    implicit none
    integer(long) data(:,:)
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)
    call mpi_recv(data, counter, mpi_integer8, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer8_2d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer8_2d

  subroutine qmpi_recv_integer8_3d(data, source, tag)
    implicit none
    integer(8) data(:,:,:)
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)
    call mpi_recv(data, counter, mpi_integer8, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer8_3d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer8_3d

  subroutine qmpi_recv_integer8_4d(data, source, tag)
    implicit none
    integer(8) data(:,:,:,:)
    integer source, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)
    call mpi_recv(data, counter, mpi_integer8, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_integer8_4d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_integer8_4d

  subroutine qmpi_recv_real4(data, source, tag)
    implicit none
    real(sp) data
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_recv(data, counter, mpi_real, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_real4 count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real4

  subroutine qmpi_recv_real8(data, source, tag)
    implicit none
    real(dp) data
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_recv(data, counter, mpi_double_precision, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_real8 count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real8

  subroutine qmpi_recv_real4_1d(data, source, tag)
    implicit none
    real(sp) data(:)
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_recv(data, counter, mpi_real, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_real4_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real4_1d

  subroutine qmpi_recv_real8_1d(data, source, tag)
    implicit none
    real(dp) data(:)
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data)
    call mpi_recv(data, counter, mpi_double_precision, source, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_real8_1d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real8_1d

  subroutine qmpi_recv_real4_2d(data, source, tag)
    implicit none
    real(sp) data(:,:)
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)

    call mpi_recv(data, counter, mpi_real, source, given_tag, mpi_comm_world, mpistatus, ierr)

    if(ierr.ne.0)then
        print *,'error recv_real4_2d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real4_2d

  subroutine qmpi_recv_real8_2d(data, source, tag)
    implicit none
    real(dp) data(:,:)
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)

    call mpi_recv(data, counter, mpi_double_precision, source, given_tag, mpi_comm_world, mpistatus, ierr)

    if(ierr.ne.0)then
        print *,'error recv_real8_2d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real8_2d

  subroutine qmpi_recv_real4_3d(data, source, tag)
    implicit none
    real(sp) data(:,:,:)
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)

    call mpi_recv(data, counter, mpi_real, source, given_tag, mpi_comm_world, mpistatus, ierr)

    if(ierr.ne.0)then
        print *,'error recv_real4_3d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real4_3d

  subroutine qmpi_recv_real8_3d(data, source, tag)
    implicit none
    real(dp) data(:,:,:)
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)

    call mpi_recv(data, counter, mpi_double_precision, source, given_tag, mpi_comm_world, mpistatus, ierr)

    if(ierr.ne.0)then
        print *,'error recv_real8_3d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real8_3d

  subroutine qmpi_recv_real4_4d(data, source, tag)
    implicit none
    real(sp) data(:,:,:,:)
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)

    call mpi_recv(data, counter, mpi_real, source, given_tag, mpi_comm_world, mpistatus, ierr)

    if(ierr.ne.0)then
        print *,'error recv_real4_4d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real4_4d

  subroutine qmpi_recv_real8_4d(data, source, tag)
    implicit none
    real(dp) data(:,:,:,:)
    integer source
    integer, optional :: tag
    integer counter, given_tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)

    call mpi_recv(data, counter, mpi_double_precision, source, given_tag, mpi_comm_world, mpistatus, ierr)

    if(ierr.ne.0)then
        print *,'error recv_real8_4d count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_real8_4d

  subroutine qmpi_recv_logical(data, target, tag)
    implicit none
    logical data
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=1
    call mpi_recv(data, counter, mpi_logical, target, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_logical count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_logical

  subroutine qmpi_recv_string(data, target, tag)
    implicit none
    character(*) data
    integer target, counter, given_tag
    integer, optional :: tag

    given_tag=0
    if(present(tag)) given_tag=tag
    counter=len(data)
    call mpi_recv(data, counter, mpi_character, target, given_tag, mpi_comm_world, mpistatus, ierr)
    if(ierr.ne.0)then
        print *,'error recv_string count=',counter,'tag=',given_tag
        stop
    end if
  end subroutine qmpi_recv_string

  subroutine qmpi_broadcast_string(string,root)
!
! send string out to all processes. if not given
! process 0 will be used as the sender - root otherwise.
!
    implicit none
    character(len=*) string
    integer, optional :: root
    integer counter,boss

    counter=len(string)

    boss=0
    if(present(root)) then
        boss=root
    end if

    call mpi_bcast(string , counter, mpi_character, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_string

  subroutine qmpi_broadcast_stringarr(data,root)
    implicit none
    character(len=*) data(:)
    integer, optional :: root
    integer counter, boss

    counter=len(data(1))*size(data)

    boss=0
    if(present(root)) then
        boss=root
    end if

    call mpi_bcast(data, counter, mpi_character, boss, mpi_comm_world ,ierr)
  end subroutine qmpi_broadcast_stringarr

  subroutine qmpi_broadcast_real4(data,root)
    implicit none
    real(4) data
    integer, optional :: root
    integer counter,boss

    counter=1
    boss=0
    if(present(root)) boss=root
    call mpi_bcast(data , counter, mpi_real, boss, mpi_comm_world, ierr)
  end subroutine qmpi_broadcast_real4

  subroutine qmpi_broadcast_real8(data,root)
    implicit none
    real(8) data
    integer, optional :: root
    integer counter,boss

    counter=1
    boss=0
    if(present(root)) boss=root
    call mpi_bcast(data , counter, mpi_double_precision, boss, mpi_comm_world, ierr)
  end subroutine qmpi_broadcast_real8

  subroutine qmpi_broadcast_integer4(data,root)
    implicit none
    integer(4) data
    integer, optional :: root
    integer counter,boss

    counter=1
    boss=0
    if(present(root)) boss=root
    call mpi_bcast(data , counter, mpi_integer, boss, mpi_comm_world, ierr)
  end subroutine qmpi_broadcast_integer4

  subroutine qmpi_broadcast_integer8(data,root)
    implicit none
    integer(8) data
    integer, optional :: root
    integer counter,boss

    counter=1
    boss=0
    if(present(root)) boss=root
    call mpi_bcast(data , counter, mpi_integer8, boss, mpi_comm_world, ierr)
  end subroutine qmpi_broadcast_integer8

  subroutine qmpi_broadcast_logical(data, root)
    implicit none
    logical data
    integer, optional :: root
    integer counter,boss

    counter=1
    boss=0
    if(present(root)) boss=root
    call mpi_bcast(data , counter, mpi_logical, boss, mpi_comm_world, ierr)
  end subroutine qmpi_broadcast_logical


  subroutine qmpi_broadcast_integer4_array1d(data,root)
    implicit none
    integer(sp) data(:)
    integer, optional :: root
    integer counter,boss

    counter=size(data)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_integer, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_integer4_array1d
  
  subroutine qmpi_broadcast_integer8_array1d(data,root)
    implicit none
    integer(long) data(:)
    integer, optional :: root
    integer counter,boss

    counter=size(data)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_integer8, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_integer8_array1d

  subroutine qmpi_broadcast_integer4_array2d(data,root)
    implicit none
    integer(sp) data(:,:)
    integer, optional :: root
    integer counter,boss

    counter=size(data,1)*size(data,2)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_integer, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_integer4_array2d
    
  subroutine qmpi_broadcast_integer8_array2d(data,root)
    implicit none
    integer(long) data(:,:)
    integer, optional :: root
    integer counter,boss

    counter=size(data,1)*size(data,2)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_integer8, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_integer8_array2d

  subroutine qmpi_broadcast_real4_array1d(data,root)
    implicit none
    real(sp) data(:)
    integer, optional :: root
    integer counter, boss

    counter=size(data)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_real, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_real4_array1d

  subroutine qmpi_broadcast_real8_array1d(data,root)
    implicit none
    real(dp) data(:)
    integer, optional :: root
    integer counter, boss

    counter=size(data)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_double_precision, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_real8_array1d

  subroutine qmpi_broadcast_real4_array2d(data,root)
    implicit none
    real(sp) data(:,:)
    integer, optional :: root
    integer counter, boss

    counter=size(data,1)*size(data,2)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data, counter, mpi_real, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_real4_array2d

  subroutine qmpi_broadcast_real8_array2d(data,root)
    implicit none
    real(dp) data(:,:)
    integer, optional :: root
    integer counter, boss

    counter=size(data,1)*size(data,2)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data, counter, mpi_double_precision, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_real8_array2d

  subroutine qmpi_broadcast_real4_array3d(data,root)
    implicit none
    real(sp) data(:,:,:)
    integer, optional :: root
    integer counter, boss

    counter=size(data,1)*size(data,2)*size(data,3)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_real, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_real4_array3d

  subroutine qmpi_broadcast_real8_array3d(data,root)
    implicit none
    real(dp) data(:,:,:)
    integer, optional :: root
    integer counter, boss

    counter=size(data,1)*size(data,2)*size(data,3)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_double_precision, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_real8_array3d

  subroutine qmpi_broadcast_real4_array4d(data,root)
    implicit none
    real(sp) data(:,:,:,:)
    integer, optional :: root
    integer counter, boss

    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_real, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_real4_array4d

  subroutine qmpi_broadcast_real8_array4d(data,root)
    implicit none
    real(dp) data(:,:,:,:)
    integer, optional :: root
    integer counter, boss

    counter=size(data,1)*size(data,2)*size(data,3)*size(data,4)
    boss=0
    if(present(root)) then
        boss=root
    end if
    call mpi_bcast(data , counter, mpi_double_precision, boss, mpi_comm_world  ,ierr)
  end subroutine qmpi_broadcast_real8_array4d

  subroutine qmpi_broadcast_real4s(a,b,c,d,e,f,root)
!
! send a,b,c,d,e,f out to all processes. if not given
! process 0 will be used as the sender - root otherwise.
!
    implicit none
    real(sp) a
    real(sp), optional :: b,c,d,e,f
    integer, optional :: root
    integer counter,boss
    real(sp) rbuff(6)

    counter=0   ;  boss=0
    if(present(root)) then
        boss=root
    end if
!    if(present(a)) then
        counter=counter+1
        rbuff(counter)=a
!    end if
    if(present(b)) then
        counter=counter+1
        rbuff(counter)=b
    end if
    if(present(c)) then
        counter=counter+1
        rbuff(counter)=c
    end if
    if(present(d)) then
        counter=counter+1
        rbuff(counter)=d
    end if
    if(present(e)) then
        counter=counter+1
        rbuff(counter)=e
    end if
    if(present(f)) then
        counter=counter+1
        rbuff(counter)=f
    end if

    call mpi_bcast(rbuff , counter, mpi_real, boss, mpi_comm_world  ,ierr)

    counter=1
    a=rbuff(counter)
    if(present(b)) then
        counter=counter+1
        b=rbuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=rbuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=rbuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=rbuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=rbuff(counter)
    end if
  end subroutine qmpi_broadcast_real4s

  subroutine qmpi_broadcast_real8s(a,b,c,d,e,f,root)
!
! send a,b,c,d,e,f out to all processes. if not given
! process 0 will be used as the sender - root otherwise.
!
    implicit none
    real(dp) a
    real(dp), optional :: b,c,d,e,f
    integer, optional :: root
    integer counter,boss
    real(kind=8) rbuff(6)

    boss=0
    if(present(root)) then
        boss=root
    end if
    counter=1
    rbuff(counter)=a
    if(present(b)) then
        counter=counter+1
        rbuff(counter)=b
    end if
    if(present(c)) then
        counter=counter+1
        rbuff(counter)=c
    end if
    if(present(d)) then
        counter=counter+1
        rbuff(counter)=d
    end if
    if(present(e)) then
        counter=counter+1
        rbuff(counter)=e
    end if
    if(present(f)) then
        counter=counter+1
        rbuff(counter)=f
    end if

    call mpi_bcast(rbuff , counter, mpi_double_precision, boss, mpi_comm_world  ,ierr)

    counter=1
    a=rbuff(counter)
    if(present(b)) then
        counter=counter+1
        b=rbuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=rbuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=rbuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=rbuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=rbuff(counter)
    end if
  end subroutine qmpi_broadcast_real8s
  
  subroutine qmpi_broadcast_logicals(a,b,c,d,e,f,root)
!
! send a,b,c,d,e,f out to all processes. if not given
! process 0 will be used as the sender - root otherwise.
!
    implicit none
    logical a
    logical, optional :: b,c,d,e,f
    integer, optional :: root
    integer counter,boss
    logical lbuff(6)

    boss=0
    if(present(root)) then
        boss=root
    end if

    counter=1
    lbuff(counter)=a
    if(present(b)) then
        counter=counter+1
        lbuff(counter)=b
    end if
    if(present(c)) then
        counter=counter+1
        lbuff(counter)=c
    end if
    if(present(d)) then
        counter=counter+1
        lbuff(counter)=d
    end if
    if(present(e)) then
        counter=counter+1
        lbuff(counter)=e
    end if
    if(present(f)) then
        counter=counter+1
        lbuff(counter)=f
    end if

    call mpi_bcast(lbuff , counter, mpi_logical, boss, mpi_comm_world  ,ierr)

    counter=1
    a=lbuff(counter)

    if(present(b)) then
        counter=counter+1
        b=lbuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=lbuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=lbuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=lbuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=lbuff(counter)
    end if
  end subroutine qmpi_broadcast_logicals

  subroutine qmpi_broadcast_integer4s(a,b,c,d,e,f,root)
!
! send a,b,c,d,e,f out to all processes. if not given
! process 0 will be used as the sender - root otherwise.
!
    implicit none
    integer(sp) a
    integer(sp), optional :: b,c,d,e,f,root
    integer counter,boss
    integer ibuff(6)

    boss=0
    if(present(root)) then
        boss=root
    end if

    counter=1
!    if(present(a)) then
!        counter=counter+1
        ibuff(counter)=a
!    end if
    if(present(b)) then
        counter=counter+1
        ibuff(counter)=b
    end if
    if(present(c)) then
        counter=counter+1
        ibuff(counter)=c
    end if
    if(present(d)) then
        counter=counter+1
        ibuff(counter)=d
    end if
    if(present(e)) then
        counter=counter+1
        ibuff(counter)=e
    end if
    if(present(f)) then
        counter=counter+1
        ibuff(counter)=f
    end if

    call mpi_bcast(ibuff , counter, mpi_integer, boss, mpi_comm_world  ,ierr)

    counter=1
    a=ibuff(counter)

    if(present(b)) then
        counter=counter+1
        b=ibuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=ibuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=ibuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=ibuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=ibuff(counter)
    end if
  end subroutine qmpi_broadcast_integer4s

  subroutine qmpi_broadcast_integer8s(a,b,c,d,e,f,root)
!
! send a,b,c,d,e,f out to all processes. if not given
! process 0 will be used as the sender - root otherwise.
!
    implicit none
    integer(long) a
    integer(long), optional :: b,c,d,e,f,root
    integer counter,boss
    integer ibuff(6)

    boss=0
    if(present(root)) then
        boss=root
    end if

    counter=1
!    if(present(a)) then
!        counter=counter+1
        ibuff(counter)=a
!    end if
    if(present(b)) then
        counter=counter+1
        ibuff(counter)=b
    end if
    if(present(c)) then
        counter=counter+1
        ibuff(counter)=c
    end if
    if(present(d)) then
        counter=counter+1
        ibuff(counter)=d
    end if
    if(present(e)) then
        counter=counter+1
        ibuff(counter)=e
    end if
    if(present(f)) then
        counter=counter+1
        ibuff(counter)=f
    end if

    call mpi_bcast(ibuff , counter, mpi_integer8, boss, mpi_comm_world  ,ierr)

    counter=1
    a=ibuff(counter)

    if(present(b)) then
        counter=counter+1
        b=ibuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=ibuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=ibuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=ibuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=ibuff(counter)
    end if
  end subroutine qmpi_broadcast_integer8s

  subroutine qmpi_real_reduction(type,a,b,c,d,e,f,root)
!
! perform a reduction of 'type' on each of the given arguments a - f. 
! if type is:
!  'sum': for each argument, return the sum of the argument over all processors
!  'mul': the product
!  'min': the minimum value
!  'max': the maximum value
! root is an optional argument, if given only return the result on that processor (reduce)
!  the default is to return the result on all processors (allreduce)
!    
    implicit none
    character(3) type
    real(sp) a
    real(sp), optional, intent(inout) :: b,c,d,e,f
    integer, optional :: root
    integer counter,boss
    integer, parameter :: dp=8
    real(dp) rbuff(6),globrbuff(6)

    if( trim(type).ne.'sum' .and. trim(type).ne.'mul' .and. trim(type).ne.'min' .and. trim(type).ne.'max')then
        print *,'qmpi.f90 reduce error: reduction of type ',type,'not supported'
        stop
    end if
    
    boss=0
    if(present(root)) boss=root

    globrbuff(:)=0.0
    counter=0
!    if(present(a)) then
        counter=counter+1
        rbuff(counter)=real(a,dp)
!    end if
    if(present(b)) then
        counter=counter+1
        rbuff(counter)=real(b,dp)
    end if
    if(present(c)) then
        counter=counter+1
        rbuff(counter)=real(c,dp)
    end if
    if(present(d)) then
        counter=counter+1
        rbuff(counter)=real(d,dp)
    end if
    if(present(e)) then
        counter=counter+1
        rbuff(counter)=real(e,dp)
    end if
    if(present(f)) then
        counter=counter+1
        rbuff(counter)=real(f,dp)
    end if

    select case(type)
    case('sum')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_sum,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        end if
    case('mul')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_prod,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_prod,mpi_comm_world,ierr)
        end if
    case('min')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_min,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_min,mpi_comm_world,ierr)
        end if
    case('max')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_max,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
        end if
    end select

    counter=0
!    if(present(a)) then
        counter=counter+1
        a=globrbuff(counter)
!    end if
    if(present(b)) then
        counter=counter+1
        b=globrbuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=globrbuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=globrbuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=globrbuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=globrbuff(counter)
    end if
  end subroutine qmpi_real_reduction

  subroutine qmpi_real8_reduction(type,a,b,c,d,e,f,root)
!
! perform a reduction of 'type' on each of the given arguments a - f. 
! if type is:
!  'sum': for each argument, return the sum of the argument over all processors
!  'mul': the product
!  'min': the minimum value
!  'max': the maximum value
! root is an optional argument, if given only return the result on that processor (reduce)
!  the default is to return the result on all processors (allreduce)
!    
    implicit none
    integer, parameter :: dp=8
    character(3) type
    real(dp) a
    real(dp), optional, intent(inout) :: b,c,d,e,f
    integer, optional :: root
    integer counter,boss
    real(dp) rbuff(6),globrbuff(6)

    if( trim(type).ne.'sum' .and. trim(type).ne.'mul' .and. trim(type).ne.'min' .and. trim(type).ne.'max')then
        print *,'qmpi.f90 reduce error: reduction of type ',type,'not supported'
        stop
    end if

    boss=0
    if(present(root))boss=root

    globrbuff(:)=0.0
    counter=0
!    if(present(a)) then
        counter=counter+1
        rbuff(counter)=a
!    end if
    if(present(b)) then
        counter=counter+1
        rbuff(counter)=b
    end if
    if(present(c)) then
        counter=counter+1
        rbuff(counter)=c
    end if
    if(present(d)) then
        counter=counter+1
        rbuff(counter)=d
    end if
    if(present(e)) then
        counter=counter+1
        rbuff(counter)=e
    end if
    if(present(f)) then
        counter=counter+1
        rbuff(counter)=f
    end if
    
    select case(type)
    case('sum')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_sum,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        end if
    case('mul')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_prod,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_prod,mpi_comm_world,ierr)
        end if
    case('min')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_min,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_min,mpi_comm_world,ierr)
        end if
    case('max')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_max,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
        end if
    end select

    counter=0
!    if(present(a)) then
        counter=counter+1
        a=globrbuff(counter)
!    end if
    if(present(b)) then
        counter=counter+1
        b=globrbuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=globrbuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=globrbuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=globrbuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=globrbuff(counter)
    end if
  end subroutine qmpi_real8_reduction

  subroutine qmpi_integer_reduction(type,a,b,c,d,e,f,root)
!
! perform a reduction of 'type' on each of the given arguments a - f. 
! if type is:
!  'sum': for each argument, return the sum of the argument over all processors
!  'mul': the product
!  'min': the minimum value
!  'max': the maximum value
! root is an optional argument, if given only return the result on that processor (reduce)
!  the default is to return the result on all processors (allreduce)
!    
    implicit none
    character(3) type
    integer(sp) a
    integer(sp), optional, intent(inout) :: b,c,d,e,f
    integer, optional :: root
    integer counter,boss
    integer rbuff(6),globrbuff(6)

    if( trim(type).ne.'sum' .and. trim(type).ne.'mul' .and. trim(type).ne.'min' .and. trim(type).ne.'max')then
        print *,'qmpi.f90 reduce error: reduction of type ',type,'not supported'
        stop
    end if

    boss=0
    if(present(root))boss=root

    globrbuff(:)=0
    counter=0
    !if(present(a)) then
        counter=counter+1
        rbuff(counter)=a
    !end if
    if(present(b)) then
        counter=counter+1
        rbuff(counter)=b
    end if
    if(present(c)) then
        counter=counter+1
        rbuff(counter)=c
    end if
    if(present(d)) then
        counter=counter+1
        rbuff(counter)=d
    end if
    if(present(e)) then
        counter=counter+1
        rbuff(counter)=e
    end if
    if(present(f)) then
        counter=counter+1
        rbuff(counter)=f
    end if

    select case(type)
    case('sum')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter, MPI_INTEGER, mpi_sum,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter, MPI_INTEGER, mpi_sum,mpi_comm_world,ierr)
        end if
    case('mul')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter, MPI_INTEGER, mpi_prod,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter, MPI_INTEGER, mpi_prod,mpi_comm_world,ierr)
        end if
    case('min')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter, MPI_INTEGER, mpi_min,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter, MPI_INTEGER, mpi_min,mpi_comm_world,ierr)
        end if
    case('max')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter, MPI_INTEGER, mpi_max,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter, MPI_INTEGER, mpi_max,mpi_comm_world,ierr)
        end if
    end select

    counter=0
!    if(present(a)) then
        counter=counter+1
        a=globrbuff(counter)
!    end if
    if(present(b)) then
        counter=counter+1
        b=globrbuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=globrbuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=globrbuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=globrbuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=globrbuff(counter)
    end if
  end subroutine qmpi_integer_reduction

  subroutine qmpi_integer8_reduction(type,a,b,c,d,e,f,root)
!
! perform a reduction of 'type' on each of the given arguments a - f. 
! if type is:
!  'sum': for each argument, return the sum of the argument over all processors
!  'mul': the product
!  'min': the minimum value
!  'max': the maximum value
! root is an optional argument, if given only return the result on that processor (reduce)
!  the default is to return the result on all processors (allreduce)
!    
    implicit none
    character(3) type
    integer(long) a
    integer(long), optional, intent(inout) :: b,c,d,e,f
    integer, optional :: root
    integer counter,boss
    integer(long) rbuff(6),globrbuff(6)

    if(len(type).ne.3)then
        print *,'qmpi.f90 reduce error: type must be one of "mul","sum","min" or "max"'
        stop
    end if
    if( trim(type).ne.'sum' .and. trim(type).ne.'mul' .and. trim(type).ne.'min' .and. trim(type).ne.'max')then
        print *,'qmpi.f90 reduce error: reduction of type ',type,'not supported'
        stop
    end if

    boss=0
    if(present(root))boss=root

    globrbuff(:)=0_dp
    counter=0
!    if(present(a)) then
        counter=counter+1
        rbuff(counter)=a
!    end if
    if(present(b)) then
        counter=counter+1
        rbuff(counter)=b
    end if
    if(present(c)) then
        counter=counter+1
        rbuff(counter)=c
    end if
    if(present(d)) then
        counter=counter+1
        rbuff(counter)=d
    end if
    if(present(e)) then
        counter=counter+1
        rbuff(counter)=e
    end if
    if(present(f)) then
        counter=counter+1
        rbuff(counter)=f
    end if

    select case(type)
    case('sum')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter, MPI_INTEGER8, mpi_sum,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter, MPI_INTEGER8, mpi_sum,mpi_comm_world,ierr)
        end if
    case('mul')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter, MPI_INTEGER8, mpi_prod,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter, MPI_INTEGER8, mpi_prod,mpi_comm_world,ierr)
        end if
    case('min')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter, MPI_INTEGER8, mpi_min,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter, MPI_INTEGER8, mpi_min,mpi_comm_world,ierr)
        end if
    case('max')
        if(present(root))then
            call mpi_reduce(rbuff,globrbuff,counter, MPI_INTEGER8, mpi_max,boss,mpi_comm_world,ierr)
        else
            call mpi_allreduce(rbuff,globrbuff,counter, MPI_INTEGER8, mpi_max,mpi_comm_world,ierr)
        end if
    end select

    counter=1
    a=globrbuff(counter)
    
    if(present(b)) then
        counter=counter+1
        b=globrbuff(counter)
    end if
    if(present(c)) then
        counter=counter+1
        c=globrbuff(counter)
    end if
    if(present(d)) then
        counter=counter+1
        d=globrbuff(counter)
    end if
    if(present(e)) then
        counter=counter+1
        e=globrbuff(counter)
    end if
    if(present(f)) then
        counter=counter+1
        f=globrbuff(counter)
    end if
  end subroutine qmpi_integer8_reduction


! later?
! packing to reduce number of sends:
  
! call pack(u)
! call pack(eta(1,:))
! call pack(v)
! call send_pack(1)
! ...
! call receive_pack(0)
! call unpack(u)
! call unpack(eta(1,:)
!
  
end module qmpi

#else

#warning "COMPILING WITHOUT QMPI CODE"

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

#endif
