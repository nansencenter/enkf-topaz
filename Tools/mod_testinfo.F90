module mod_testinfo


  type testinfo
     real             :: min, max 
     character(len=8) :: fldname
     integer          :: toosmall, toolarge
  end type testinfo

contains


  subroutine tests_init(tests, numtest)
    implicit none
    integer, intent(in) ::  numtest
    type(testinfo) , dimension(numtest), intent(out) :: tests
    integer             :: k
    tests(:)%toosmall=0
    tests(:)%toolarge=0
    tests(:)%fldname=''
    do k=1,numtest
       select case (k)
       case (1)
          tests(k)%fldname='temp'
          tests(k)%max    =40.
          tests(k)%min    =-4.8
       case (2)
          tests(k)%fldname='saln'
          tests(k)%max    =40.
          tests(k)%min    =0.   
       case (3)
          tests(k)%fldname='dp'
          tests(k)%max    =10000*9806.0
          tests(k)%min    =0.0
       case (4)
          tests(k)%fldname='u'
          tests(k)%max    =4
          tests(k)%min    =-4
       case (5)
          tests(k)%fldname='v'
          tests(k)%max    =4
          tests(k)%min    =-4
       case (6)
          tests(k)%fldname='ubavg'
          tests(k)%max    =2
          tests(k)%min    =-2
       case (7)
          tests(k)%fldname='vbavg'
          tests(k)%max    =2
          tests(k)%min    =-2
       case (8)
#if defined (SINGLE_RESTART)
          tests(k)%fldname='ficem'
#else
          tests(k)%fldname='icec'
#endif
          tests(k)%min    =0.0
          tests(k)%max    =1.0
       case (9)
#if defined (SINGLE_RESTART)
          tests(k)%fldname='hicem'
#else
          tests(k)%fldname='hice'
#endif
          tests(k)%min    =0.0
          tests(k)%max    =20.0
       case default
          print *,'Not set up test for k=',k
       end select
    end do
  end subroutine tests_init

  subroutine matchtest(cfld,tests,numtest,testindex)
    implicit none
    integer, intent(in) :: numtest
    character(len=*), intent(in) :: cfld
    type(testinfo), dimension(numtest), intent(in) :: tests
    integer, intent(out) :: testindex

    integer :: i

    testindex=-1
    do i=1,numtest
       if (trim(cfld)==trim(tests(i)%fldname)) testindex=i
    end do
  end subroutine matchtest




end module mod_testinfo

