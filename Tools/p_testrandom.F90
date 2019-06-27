program testrandom
! checks random numbers

! For parallelization
   use m_set_random_seed2
   implicit none
   integer seedsze, i, maxnb
   real, allocatable :: putseed(:)
   real rand, frac
   integer histo(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set a variable random seed
   call set_random_seed2

   maxnb=1000000
! Remove any randomness in the results for Parallab
!   call random_seed(size=seedsze)
!   allocate(putseed(seedsze))
!   putseed(:)=13
!   call random_seed(put=putseed)
!   deallocate(putseed)

histo=0
do i=1,maxnb
   call random_number(rand)
   do i=1,10 
     frac=1./float(i)
     if (rand>frac-0.1 .and. rand<frac) histo(i)=histo(i)+1 
   enddo

!   if (mod(i,100)==0) then 
!      print *, rand 
!   endif
enddo
print *,histo(:)/float(maxnb)

end program
