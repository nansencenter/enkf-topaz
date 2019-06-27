      subroutine swapendian(a)
      implicit none
      integer(kind=8), intent(inout) :: a  ! 4-bytes

      integer(kind=8) ii4,   io4     ! 4-bytes
      common/czioxe/  ii4,   io4     ! helps prevent unwanted optimization
      save  /czioxe/
      integer(kind=1) ii1(8),io1(8)  ! 1-byte
      equivalence    (ii4,ii1(1)), (io4,io1(1))  ! non-standard f90

        ii4 = a
        io1(1) = ii1(8)
        io1(2) = ii1(7)
        io1(3) = ii1(6)
        io1(4) = ii1(5)
        io1(5) = ii1(4)
        io1(6) = ii1(3)
        io1(7) = ii1(2)
        io1(8) = ii1(1)
        a = io4
      return
      end subroutine swapendian

      subroutine swapendian2(a,n)
      implicit none
      integer        , intent(in)    :: n    ! Size of input type to convert

      ! NB - input can be anything - can not be compiled with input argument checking
      integer(kind=1), intent(inout) :: a(n) 

      integer k
      integer(kind=1) ii4(16),   io4(16)     ! 16 bytes should beenough for everyone
      !common/czioxe/  ii4,   io4     ! helps prevent unwanted optimization
      !save  /czioxe/
      !integer(kind=1) ii1(16),io1(16)  ! 1-byte
      !equivalence    (ii4(1),ii1(1)), (io4(1),io1(1))  ! non-standard f90

        ii4(1:n) = a

        do k=1,n
           !io1(k) = ii1(n-k+1)
           io4(k) = ii4(n-k+1)
        end do

        a = io4(1:n)
      return
      end subroutine swapendian2

