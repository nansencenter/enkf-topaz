      module mod_raw_io
      contains


! Modified from Alan Wallcraft's RAW routine by Knut Liseter @ NERSC
! So far only the "I" in "IO" is present
      SUBROUTINE READRAW(A,AMN,AMX,IDM,JDM,LSPVAL,SPVAL,CFILE1,K)
      IMPLICIT NONE
C
      REAL*4     SPVALH
      !PARAMETER (SPVALH=2**100)
      PARAMETER (SPVALH=1.0E30_4)
C
      REAL*4,        INTENT(OUT) :: A(IDM,JDM)
      REAL*4,        INTENT(OUT) :: AMN,AMX
      INTEGER,       INTENT(IN)  :: IDM,JDM
      LOGICAL,       INTENT(IN)  :: LSPVAL
      REAL*4,        INTENT(INOUT)  :: SPVAL
      INTEGER,       INTENT(IN)  :: K
      CHARACTER(len=*), INTENT(IN)  :: CFILE1
C
      REAL*4 :: PADA(4096)
C
C     MOST OF WORK IS DONE HERE.
C

      INTEGER      LEN_TRIM
      INTEGER      I,J,IOS,NRECL
      INTEGER NPAD
C
      IF(.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
C
!!! Calculate the number of elements padded!!!!!!!!!!!!!!!!!!!!!!!!
      NPAD=GET_NPAD(IDM,JDM)
C
      INQUIRE( IOLENGTH=NRECL) A,PADA(1:NPAD)
C     
C     
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='old',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS, ACTION='READ')
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      READ(11,REC=K,IOSTAT=IOS) A
      close(11)
C
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t read record ',K,
     &             ' from '//CFILE1(1:LEN_TRIM(CFILE1))
        CALL EXIT(4)
      ENDIF
C
      AMN =  SPVALH
      AMX = -SPVALH
      DO J= 1,JDM
      DO I=1,IDM
         IF     (A(I,J).LE.SPVALH) THEN
            AMN = MIN(real(AMN, 4), real(A(I,J), 4))
            AMX = MAX(real(AMX, 4), real(A(I,J), 4))
         ELSEIF (LSPVAL) THEN
            A(I,J) = SPVAL
         ENDIF
      END DO
      END DO
C                 
      RETURN
      END SUBROUTINE

! Modified from Alan Wallcraft's RAW routine by Knut Liseter @ NERSC
! This wll be the  "O" in "IO" is present
      SUBROUTINE WRITERAW(A,AMN,AMX,IDM,JDM,LSPVAL,SPVAL,CFILE1,K)
      IMPLICIT NONE
C
      REAL*4     SPVALH
      !PARAMETER (SPVALH=2**100)
      PARAMETER (SPVALH=1.0e30_4)
C
      REAL*4,        INTENT(INOUT) :: A(IDM,JDM)
      REAL*4,        INTENT(OUT)   :: AMN,AMX
      INTEGER,       INTENT(IN)    :: IDM,JDM
      LOGICAL,       INTENT(IN)    :: LSPVAL
      REAL*4,        INTENT(INOUT) :: SPVAL
      INTEGER,       INTENT(IN)    :: K
      CHARACTER(len=*), INTENT(IN) :: CFILE1
C
      REAL*4 :: PADA(4096)
C
C     MOST OF WORK IS DONE HERE.
C

      INTEGER      LEN_TRIM
      INTEGER      I,J,IOS,NRECL
      INTEGER NPAD
C
      IF(.NOT.LSPVAL) THEN
        SPVAL = SPVALH
      ENDIF
C
!!! Calculate the number of elements padded!!!!!!!!!!!!!!!!!!!!!!!!
      NPAD=GET_NPAD(IDM,JDM)
C
      PADA=0.
      INQUIRE( IOLENGTH=NRECL) A,PADA(1:NPAD)
C     
C     
      OPEN(UNIT=11, FILE=CFILE1, FORM='UNFORMATTED', STATUS='unknown',
     +         ACCESS='DIRECT', RECL=NRECL, IOSTAT=IOS)
      IF     (IOS.NE.0) THEN
        write(6,*) 'Error: can''t open ',CFILE1(1:LEN_TRIM(CFILE1))
        write(6,*) 'ios   = ',ios
        write(6,*) 'nrecl = ',nrecl
        CALL EXIT(3)
      ENDIF
C
      WRITE(11,REC=K,IOSTAT=IOS) A,PADA(1:NPAD)
      close(11)
C
      IF     (IOS.NE.0) THEN
        WRITE(6,*) 'can''t write record ',K,
     &             ' from '//CFILE1(1:LEN_TRIM(CFILE1))
        CALL EXIT(4)
      ENDIF
C
      AMN =  SPVALH
      AMX = -SPVALH
      DO J= 1,JDM
      DO I=1,IDM
         IF     (A(I,J).LE.SPVALH) THEN
            AMN = MIN(real(AMN, 4), real(A(I,J), 4))
            AMX = MAX(real(AMX, 4), real(A(I,J), 4))
         ELSEIF (LSPVAL) THEN
            A(I,J) = SPVAL
         ENDIF
      END DO
      END DO
C                 
      RETURN
      END SUBROUTINE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to get index of fields in data file (.a) from header file (.b)
      subroutine rst_index_from_header(fname,cfld,vlevel,tlevel,
     &                                 indx,bmin,bmax,skiphdr) 
      implicit none
      character(len=*), intent(in) :: fname     ! filename without extention
      character(len=*), intent(in) :: cfld      ! variable name
      integer         , intent(in) :: tlevel    ! time level
      integer         , intent(in) :: vlevel    ! vertical level
      integer         , intent(out):: indx      ! index in .a file
      real            , intent(out):: bmin,bmax ! min and max from b file
      logical         , intent(in) :: skiphdr

      integer :: itlevel, ivlevel
      character(len=8) :: icfld
      integer :: ios,i
      integer :: nskip_rst,nop
      logical :: match, ex

      nskip_rst=2
      nop = 999

      ! Open file
      inquire(exist=ex,file=trim(fname))
      if (.not. ex) then
         print *,'file '//trim(fname)//' is not present'
         call exit(1)
      end if
      open(nop,file=trim(fname),status='old',action='read')

      ! Skip first nskip lines
      if (skiphdr) then
         do i=1,nskip_rst
            read(nop,*)
         end do
      end if

      match=.false.
      indx=0
      ios=0
      do while (ios==0 .and. .not.match)
         read(nop,117,iostat=ios) icfld,ivlevel,itlevel,bmin,bmax
         match= icfld==cfld .and. ivlevel==vlevel .and. itlevel==tlevel
         indx=indx+1
         !print *,icfld,itlevel,ivlevel,bmin,bmax
      end do

      close(nop)

      if (.not.match) then
         !print *,'Could not find field '//cfld
         !print *,'Vertical level :',vlevel
         !print *,'Time     level :',tlevel
         indx=-1
         !call exit(1) ! Always return to caller
      endif

  117 format (a8,23x,i3,i3,2x,2e16.7)

      end subroutine

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to get field desc in  header file (.b) from index in data file (.a)
      subroutine rst_header_from_index(fname,cfld,vlevel,tlevel,
     &                                 indx,bmin,bmax,skiphdr) 
      implicit none
      character(len=*), intent(in)  :: fname     ! filename without extention
      character(len=8), intent(out) :: cfld      ! variable name
      integer         , intent(out) :: tlevel    ! time level
      integer         , intent(out) :: vlevel    ! vertical level
      integer         , intent(in)  :: indx      ! index in .a file
      real            , intent(out) :: bmin,bmax ! min and max from b file
      logical         , intent(in ) :: skiphdr   ! Skip header of .b file

      integer :: ios,i
      integer :: nskip_rst,nop
      logical :: ex


      nskip_rst=2
      nop = 999

      ! Open file
      inquire(exist=ex,file=trim(fname))
      if (.not. ex) then
         print *,'file '//trim(fname)//' not present'
         call exit(1)
      end if
      open(nop,file=trim(fname),status='old',action='read')


      ! Skip first nskip + index-1 lines
      !print *,'hei'
      if (skiphdr) then
         do i=1,nskip_rst
            read(nop,*)
         end do
      end if
      do i=1,indx-1
         read(nop,*)
      end do
      read(nop,117,iostat=ios) cfld,vlevel,tlevel,bmin,bmax
      close(nop)

      if (ios/=0) then
         !print *,'Could not get info from  index',indx
         !call exit(1)
         cfld=''
         tlevel=-1
         vlevel=-1
      endif

  117 format (a8,23x,i3,i3,2x,2e16.7)

      end subroutine
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to get index of fields in regional grid file (.a) from header file (.b)
      subroutine grid_index_from_header(fname,cfld,indx,bmin,bmax
     &                                  ,skiphdr) 
      implicit none
      character(len=*), intent(in) :: fname     ! filename without extention
      character(len=*), intent(in) :: cfld      ! variable name
      integer         , intent(out):: indx      ! index in .a file
      real            , intent(out):: bmin,bmax ! min and max from b file
      logical         , intent(in) :: skiphdr

      character(len=4) :: icfld
      character*80 :: cline
      integer :: ios,i
      integer :: nskip_grid,nop
      logical :: match, ex

      nskip_grid=3
      nop = 999

      ! Open file
      inquire(exist=ex,file=trim(fname))
      if (.not. ex) then
         print *,'file '//trim(fname)//' is not present'
         call exit(1)
      end if
      open(nop,file=trim(fname),status='old',action='read')


      ! Skip first nskip lines
      if (skiphdr) then
         do i=1,nskip_grid
            read(nop,*)
         end do
      end if

      match=.false.
      indx=0
      ios=0
      do while (ios==0 .and. .not.match)
         read(nop,'(a)') cline
         icfld=cline(1:4)
         i=index(cline,'=')
         read (cline(i+1:),*) bmin,bmax
         match= trim(icfld)==trim(cfld)
         indx=indx+1
      end do

      close(nop)

      if (.not.match) then
         indx=-1
      endif
      end subroutine grid_index_from_header


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to get index of fields in regional grid file (.a) from header file (.b)
      subroutine daily_index_from_header(fname,cfld,coord,indx,
     &   bmin,bmax,iflag)
      implicit none
      character(len=*), intent(in) :: fname     ! filename without extention
      character(len=*), intent(in) :: cfld      ! variable name
      integer         , intent(in) :: coord     ! vertical coordinate
      integer         , intent(out):: indx      ! index in .a file
      real            , intent(out):: bmin,bmax ! min and max from b file
      integer,          intent(in) :: iflag     ! used for debuging

      logical, parameter:: skiphdr=.false.
      character(len=5) :: char5
      character(len=8) :: char8
      character(len=80) :: tmpstr 
      integer :: ios
      integer :: nop
      logical :: match, ex
      real    :: dens,rday
      integer :: lcoord,nstep

      nop = 999

      ! Open file
      inquire(exist=ex,file=trim(fname))
      if (.not. ex) then
         print *,'file '//trim(fname)//' is not present'
         call exit(1)
      end if

      open(nop,file=trim(fname),status='old')
      ! Skip first nskip lines
      if (skiphdr) then
         do while (char5/='field' .and. ios==0)
            read(nop,'(a5)',iostat=ios) char5
         end do
      end if

      ! Read until we get the field we want
      indx=0
      ios=0
      char8=''
      lcoord=-1
      match=.false.
      do while(.not.match .and. ios==0)
         read(nop,118,iostat=ios) char8,tmpstr,lcoord,tmpstr,bmin,bmax
         match=(trim(cfld)==trim(char8) .and. lcoord==coord)
         if (match) then
            write(6,119) char8,tmpstr,lcoord,tmpstr,bmin,bmax,iflag
         endif
         indx=indx+1
      end do
      close(nop)

      if (.not.match) then
         indx=-1
      endif

  119 format (a8,' = ',a21,i3,a7,1p2e16.7,i5)
  118 format (a8,' = ',a21,i3,a7,1p2e16.7)
  117 format (a8,' = ',i11,f11.2,i3,f7.3,1p2e16.7)
      end subroutine daily_index_from_header


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







      INTEGER FUNCTION GET_NPAD(IDM,JDM)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDM,JDM
         GET_NPAD = 4096 - MOD(IDM*JDM,4096)
         GET_NPAD = mod(GET_NPAD,4096)
      END FUNCTION
      end module mod_raw_io
