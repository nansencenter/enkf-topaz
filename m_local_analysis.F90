! File:          m_local_analysis.F90
!
! Created:       L. Bertino, 2002
!
! Last modified: 13/04/2010
!
! Purpose:       Local analysis:
!                  -- calculation of X5
!                  -- update of the ensemble fields
!
! Description:   This module handles local analysis.
!
! Modifications:
!                20/9/2011 PS:
!                    - modified update_fields() to allow individual inflation
!                      for each of `nfields' fields - thanks to Ehouarn Simon
!                      for spotting this inconsistency
!                25/8/2010 PS:
!                    - "obs" and "nobs" are now global, stored in m_obs. 
!                      Accordingly, the local observations variables are now
!                      called "lobs" and "nlobs". Renamed "DD" to "D" and "d"
!                      to "dy". 
!                5/8/2010 PS:
!                    - moved applying inflation from calc_X5() to
!                      update_fields()
!                    - introduced "rfactor" argument to calc_X5() - increases
!                      obs. error variance for the update of anomalies.
!                29/7/2010 PS:
!                    - calc_X5(): updated the list of things that needed to be
!                      done for a point with no local obs
!                6/7/2010 PS:
!                    - moved ij2nc() to p2nc_writeobs() in m_point2nc.F90
!                19/6/2010 PS:
!                    - added X5 to the ij2nc() output
!                25/5/2010 PS:
!                    - modified to accommodate inflation
!                    - modified to calculate SRF (spread reduction factor)
!                13/4/2010 Alok Gupta: added open/close/barrier to ensure that
!                    X5tmp.uf exists before any node tries to access it.
!                8/4/2010 PS: replaced "X4" by "X5"; renamed "localanalysis()"
!                    to "update_fields()", and "pre_local_analysis()" by
!                    "calc_X5"
!                1/03/2010 PS:
!                  - Additional checks for file I/O, as the X4 read/write have
!                    been demonstrated to fail occasionally. A record is now
!                    written to X4tmp, then read back and compared until the
!                    two instances coincide (10 attempts max).
!                11/11/2009 PS:
!                  - Changed numerics. Now it is always assumed that R is 
!                    diagonal
!                  - Choice of two chemes: EnKF and DEnKF (for now)
!                  - X4 calculated either in ens or obs space, depending on
!                    relation between nobs (# of local observations) and nrens
!                  - dfs and nobs for each (i,j) are written to enkf_diag.nc
!                  - if TEST = .true. then local stuff for (I,J) around
!                    (TEST_I, TEST_J) is dumped to enkf_<I>,<J>.nc
!                6/3/2008 PS:
!                  - in pre_local_analysis():
!                    - introduced quick sort (O(n log n)) of pre-selected
!                      observations
!                    - reshuffled the interface
!                    - replaced output array of flags for local obs by an array
!                      of indices
!                  - in local_analysis():
!                      -- unified arrays subD and subS
!                      -- got rid of calls to getD()
!                      -- used matmul()
!                      -- introduced localisation function
!                      -- eliminated X2 and V
!                2007 K. A. Liseter and Ragnhild Blikberg:
!                      -- MPI parallelisation

module m_local_analysis
  implicit none

  !
  ! public stuff
  !
  real(4), allocatable, public :: X5(:,:,:)
  real(4), allocatable, public :: X5check(:,:,:)

  public calc_X5
  public update_fields

  integer, parameter, private :: STRLEN = 512
  integer, parameter, private :: MAXITER = 10

  integer, private :: nX5pad
  real(4), allocatable, private :: X5pad(:)

  private get_npad_la
  private locfun
  private get_local_obs
  private diag2nc
  private traceprod
 
  !
  ! available localisation functions
  !
  integer, parameter, private :: LOCFUN_NONE = 1
  integer, parameter, private :: LOCFUN_STEP = 2
  integer, parameter, private :: LOCFUN_GASPARI_COHN = 3

  !
  ! used localisation function
  !
  integer, private :: LOCFUN_USED = LOCFUN_GASPARI_COHN

  !
  ! available schemes
  !
  integer, parameter, private :: SCHEME_ENKF = 1
  integer, parameter, private :: SCHEME_ETKF = 2 ! not implemented
  integer, parameter, private :: SCHEME_DENKF = 3

  !
  ! used scheme
  !
  integer, private :: SCHEME_USED = SCHEME_DENKF

contains 

  ! This routine is called for each "field" (horizontal slab) after calcX5().
  ! It conducts the multiplication
  !   E^a(i, :) = E^f(i, :) * X5(i), i = 1,...,n,
  ! where n - state dimension.
  !
  ! In this package the localisation is conducted only horizontally, so that
  ! the local (nrens x nrens) ensemble transform matrix X5 is stored for each
  ! node of the horizontal model grid. In TOPAZ4 this requires  
  ! 880 x 800 x 100 x 100 x 4 = 28 GB of storage on disk for "tmpX5.uf". If the
  ! fileds were updated on one-by-one basis, this file would have to be read
  ! (in TOPAZ4) 146 times. Therefore, the fields are updated in bunches of
  ! `nfields' to reduce the load on disk.
  !
  subroutine update_fields(ni, nj, nrens, nfields, nobs_array, depths, fld, infls)
#if defined (QMPI)
    use qmpi
#else
    use qmpi_fake
#endif
    use mod_measurement
    implicit none

    integer, intent(in) :: ni, nj ! size of grid
    integer, intent(in) :: nrens ! size of ensemble
    integer, intent(in) :: nfields ! number of 2D fields to be updated
    integer, dimension(ni, nj), intent(in) :: nobs_array! number of local obs
    real, dimension(ni, nj), intent(in) :: depths 
    real(4), dimension(ni * nj, nrens * nfields), intent(inout) :: fld ! fields
    real, dimension(nfields), intent(in) :: infls ! inflation factors

    real(4), dimension(nrens, nrens) :: X5tmp
    real(4), dimension(nrens, nrens) :: IM ! inflation matrix

    integer :: m, i, j, f
    integer :: irecl, iostatus
    real(4) :: infl

    !KAL -- all nodes open for read access to temporary "X5" file 
    inquire(iolength = irecl) X5(1 : nrens, 1 : nrens, 1 : ni), X5pad
    open(17, file = 'tmpX5.uf', form = 'unformatted', access = 'direct',&
         status = 'old', recl = irecl)

    do j = 1, nj
       ! read X5 from disk
       read(17, rec = j, iostat = iostatus) X5
       if (iostatus /= 0) then 
          print *, 'ERROR: local_analysis(): I/O error at reading X5, iostatus = ', iostatus
          print *, 'ERROR: at j = ', j
          stop
       end if
 
       do i = 1, ni
          ! skip this cell if it is on land
          if (depths(i,j) <= 0.0) then
             cycle
          end if

          if (nobs_array(i, j) == 0 .and. all(infls == 1.0d0)) then
             cycle
          end if

          X5tmp = X5(:, :, i)
          do m = 1, nrens
             if (abs(1.0e0 - sum(X5tmp(:, m))) > 1.0e-5) then
                print *, 'ERROR: detected inconsistency in X5'
                print *, 'ERROR: at j = ', j, 'i = ', i
                print *, 'ERROR: sum(X5(:, ', m, ') = ', sum(X5tmp(:, m))
                stop
             end if
          enddo

          ! ensemble transformation, in real(4)
          !
          do f = 1, nfields
             infl = infls(f) ! conversion to real(4)
             if (infl /= 1.0) then
                IM = - (infl - 1.0) / real(nrens, 4)
                do m = 1, nrens
                   IM(m, m) = IM(m, m) + infl
                end do

                fld((j - 1) * ni + i, (f - 1) * nrens + 1 : f * nrens) =&
                     matmul(fld((j - 1) * ni + i, (f - 1) * nrens + 1 : f * nrens),&
                     matmul(X5tmp, IM))
             else
                fld((j - 1) * ni + i, (f - 1) * nrens + 1 : f * nrens) =&
                     matmul(fld((j - 1) * ni + i, (f - 1) * nrens + 1 : f * nrens), X5tmp)
             end if
          end do
       enddo
    enddo
    close(17)
  end subroutine update_fields


  ! This routine calculates X5 matrices involved in the EnKF analysis, 
  !   E^a(i, :) = E^f(i, :) * X5(i), i = 1,...,n,
  ! where n - state dimension.
  !
  ! X5(i) is calculated locally (for a given state element i) as 
  !   X5 = I + G s 1^T + T,
  ! where
  !   G = S^T (I + S S^T)^{-1} = (I + S^T S)^{-1} S^T
  !   T = I - 1/2 G S        (DEnKF)
  !   T = I + G(D - S)       (EnKF)
  !   T = (I + S^T S)^{-1/2} (ETKF)
  !   S = R^{-1/2} HA^f / sqrt(m - 1)
  !   s = R^{-1/2} (d - Hx^f) / sqrt(m - 1)
  !
  !   see Sakov et al. (2010): Asynchronous data assimilation with the EnKF,
  !   Tellus 62A, 24-29.
  !
  subroutine calc_X5(nrens, modlon, modlat, depths, mindx, meandx, dy, S,&
       radius, rfactor, nlobs_array, ni, nj)
#if defined (QMPI)
    use qmpi
#else
    use qmpi_fake
#endif
    use m_parameters
    use distribute
    use mod_measurement
    use m_obs
    use m_spherdist
    use m_random
    use m_point2nc
    implicit none

    ! Input/output arguments
    integer, intent(in) :: nrens
    real, dimension(ni, nj), intent(in) :: modlon, modlat
    real, dimension(ni, nj), intent(in) :: depths
    real, intent(in) :: mindx ! min grid size
    real, intent(in) :: meandx ! mean grid size
    real, dimension(nobs), intent(inout) :: dy ! innovations
    real, dimension(nobs, nrens), intent(inout) :: S ! HA
    real, intent(in) :: radius ! localisation radius in km
    real, intent(in) :: rfactor ! obs. variance multiplier for anomalies
    integer, dimension(ni, nj), intent(out) :: nlobs_array ! # of local obs
                                                           ! for each grid cell
    integer, intent(in) :: ni, nj ! horizontal grid size

    real, dimension(nrens, nrens) :: X5tmp
    integer, dimension(nobs) :: lobs ! indices of local observations

    real, allocatable, dimension(:,:) :: D ! observation perturbations
    real, allocatable, dimension(:) :: subdy
    real, allocatable, dimension(:) :: lfactors ! loc. coeffs stored for QC
    real, allocatable, dimension(:,:) :: subD, subS ! nobs x nrens
    real, allocatable, dimension(:,:) :: X1 ! nobs x nobs
    real, allocatable, dimension(:,:) :: G
    real, allocatable, dimension(:) :: x
    real :: sqrtm
    real :: tmp(1)

    integer :: iostatus
    integer, dimension(nj):: jmap, jmap_check
#if defined (QMPI)
    integer, allocatable, dimension(:, :) :: mpibuffer_int
    real(4), allocatable, dimension(:, :) :: mpibuffer_float1, mpibuffer_float2
#endif

    integer :: lapack_info

#if defined (QMPI)
    integer :: p
#endif
    integer :: nlobs ! # of local obs
    integer :: m, i, j, o, jj, iter
    logical :: testthiscell ! test analysis at a certain cell
    integer :: irecl
    integer :: nlobs_max ! maximal number of local obs
    real :: dist, lfactor
    type(measurement) :: obs0

    ! dfs calculation
    real :: dfs
    real(4) :: dfs_array(ni, nj)
    ! srf calculation
    real :: srf
    real(4) :: srf_array(ni, nj)

    ! "partial" dfs
    real :: pdfs(nuobs)
    real(4) :: pdfs_array(ni, nj, nuobs)
    ! "partial" srf
    real :: psrf(nuobs)
    real(4) :: psrf_array(ni, nj, nuobs)
    ! auxiliary variables for dfs and srf calculation, such as
    ! nobs for different obs types
    integer :: plobs(nobs, nuobs)
    integer :: pnlobs(nuobs)
    integer :: uo
    real    :: Typ     ! used by localization radius

    if (trim(METHODTAG) == "ENKF") then
       SCHEME_USED = SCHEME_ENKF
    elseif (trim(METHODTAG) == "DENKF") then
       SCHEME_USED = SCHEME_DENKF
    end if

    if (master) then
       if (SCHEME_USED == SCHEME_ENKF) then
          print *, 'using EnKF analysis scheme'
       elseif (SCHEME_USED == SCHEME_DENKF) then
          print *, 'using DEnKF analysis scheme'
       end if
    end if

    if (LOCRAD > 0.0d0) then
       if (trim(LOCFUNTAG) == "GASPARI-COHN"&
            .or. trim(LOCFUNTAG) == "GASPARI_COHN") then
          LOCFUN_USED = LOCFUN_GASPARI_COHN
       elseif (trim(LOCFUNTAG) == "STEP") then
          LOCFUN_USED = LOCFUN_STEP
       elseif (trim(LOCFUNTAG) == "NONE") then
          LOCFUN_USED = LOCFUN_NONE
       end if
    else
       LOCFUN_USED = LOCFUN_NONE
    end if

    if (master) then
       if (LOCFUN_USED ==  LOCFUN_GASPARI_COHN) then
          print *, 'using Gaspari-Cohn localisation'
       elseif (LOCFUN_USED ==  LOCFUN_STEP) then
          print *, 'using STEP localisation'
       elseif (LOCFUN_USED ==  LOCFUN_NONE) then
          print *, 'using NO localisation'
       end if
    end if

    sqrtm = sqrt(real(nrens) - 1.0d0)
    if (SCHEME_USED == SCHEME_ENKF) then
       allocate(D(nobs, nrens))
       do o = 1, nobs
          call randn(nrens, D(o, :))
          D(o, :) = D(o, :) / (rfactor * sqrtm)
       end do
    end if
    do o = 1, nobs
       S(o, :) = S(o, :) / (sqrt(obs(o) % var) * sqrtm)
       dy(o) = dy(o) / (sqrt(obs(o) % var) * sqrtm)
    end do

    ! Distribute loops across MPI nodes
    call distribute_iterations(nj)

    ! The binary file tmpX5.uf holds (ni x nj) local ensemble transform
    ! matrices X5, (nrens x nrens) each. They are used for creating the 
    ! analysed ensemble in local_analysis(). In TOPAZ3 tmpX5.uf takes about
    ! 30GB of the disk space.
    !
    nX5pad = get_npad_la(nrens * nrens, ni)
    allocate(X5pad(nX5pad))
    inquire(iolength = irecl) X5, X5pad

    if (master) then
       open(17, file = 'tmpX5.uf', form = 'unformatted', access = 'direct', status = 'unknown', recl = irecl)
       ! get the necessary space on disk, before starting simultaneous writing
       ! by different nodes
       write(17, rec = nj) X5
       close(17)
    end if
#if defined (QMPI)
    call barrier()
#endif
    open(17, file = 'tmpX5.uf', form = 'unformatted', access = 'direct',&
         status = 'old', recl = irecl)

    open(31, file = trim(JMAPFNAME), status = 'old', iostat = iostatus)
    if (iostatus /= 0) then
       if (master) then
          print *, 'WARNING: could not open jmap.txt for reading'
          print *, '         no re-mapping of grid rows performed'
       end if
       do j = 1, nj
          jmap(j) = j
       end do
    else
       read(31, *, iostat = iostatus) jmap
       if (iostatus /= 0) then
          print *, 'ERROR reading jmap.txt'
          stop
       end if
       close(31)
       jmap_check = 1
       jmap_check(jmap) = 0
       if (sum(jmap_check) /= 0) then
          print *, 'ERROR: non-zero control sum for jmap =', sum(jmap_check)
          stop
       end if
    end if

    ! main cycle (over horizontal grid cells)
    !
    dfs_array = 0.0
    pdfs_array = 0.0
    srf_array = 0.0
    psrf_array = 0.0
    nlobs_array = 0
    do jj = my_first_iteration, my_last_iteration
       j = jmap(jj)
       print *, 'calc_X5(): jj =', jj, 'j =', j

       do i = 1, ni
          ! data dumping flag
          testthiscell = p2nc_testthiscell(i, j)

          if (testthiscell) then
             print *, 'testthiscell: depth(,', i, ',', j, ') =', depths(i, j)
          end if

          if (depths(i, j) > 0.0d0) then
             nlobs = 0 ! no upper limit on the number of local observations
             call get_local_obs(i, j, radius * 1000.0, modlon, modlat,&
                  mindx, ni, nj, nlobs, lobs)
             nlobs_array(i, j) = nlobs
          else
             nlobs = 0
          end if

          if (testthiscell) then
             print *, 'testthiscell: nlobs(,', i, ',', j, ') =', nlobs
          end if

          if (nlobs == 0) then
             ! just in case
             X5(:, :, i) = 0.0
             X5tmp = 0.0d0
             do m = 1, nrens
                X5(m, m, i) = 1.0
                X5tmp(m, m) = 1.0d0
             enddo
             if (testthiscell) then
                tmp(1) = rfactor
                call p2nc_writeobs(i, j, nlobs, nrens, X5tmp, modlon(i, j),&
                     modlat(i, j), depths(i, j), tmp(1), lobs(1 : nlobs), &
                     obs(lobs(1 : nlobs)), x, subS, subdy, lfactors)
             end if
             dfs_array(i, j) = 0.0
             pdfs_array(i, j, :) = 0.0
             srf_array(i, j) = 0.0
             psrf_array(i, j, :) = 0.0
             cycle
          end if

          if (nlobs < 0) then ! an extra check on the C-Fortran interface
             print *, 'ERROR: nlobs =', nlobs, ' for i, j =', i, j
             call stop_mpi()
          end if

          ! Allocate local arrays
          if (SCHEME_USED == SCHEME_ENKF) then
             allocate(subD(nlobs, nrens))
          end if
          allocate(subdy(nlobs))
          allocate(lfactors(nlobs))
          allocate(subS(nlobs, nrens))
          ! ( BTW subS1 = subS / sqrt(rfactor) )
          allocate(G(nrens, nlobs))
          if (nlobs < nrens) then
             allocate(X1(nlobs, nlobs))
          else
             allocate(X1(nrens, nrens))
          end if

          if (SCHEME_USED == SCHEME_ENKF) then
             subD = D(lobs(1 : nlobs), :)
          end if
          subS = S(lobs(1 : nlobs), :)
          subdy = dy(lobs(1 : nlobs))

          ! taper ensemble observation anomalies and innovations
          !
          if (LOCFUN_USED /= LOCFUN_NONE) then
             do o = 1, nlobs
                obs0 = obs(lobs(o))
                dist = spherdist(modlon(i, j), modlat(i, j),&
                     obs0 % lon, obs0 % lat)
                Typ=Typobs(lobs(o))
                lfactor = locfun(dist / (radius*Typ) / 1000.0)
                subS(o, :) = subS(o, :) * lfactor
                subdy(o) = subdy(o) * lfactor
                lfactors(o) = lfactor
                
                if (SCHEME_USED == SCHEME_ENKF) then
                   subD(o, :) = subD(o, :) * lfactor
                end if
             end do
          else
             lfactors = 1
          end if

          ! first iteration - with rfactor = 1, for the update of the mean
          ! secons iteration - with the specified rfactorm for the update of
          ! the anomalies
          !
          do iter = 1,2
             if (iter == 2) then
                if (rfactor == 1.0d0) then
                   go to 10
                end if
                subS = subS / sqrt(rfactor)
             end if

             if (nlobs < nrens) then ! use observation space
                ! Construct matrix (S * S' + I) - to be inverted
                !
                X1 = matmul(subS, transpose(subS))
                do o = 1, nlobs
                   X1(o, o) = X1(o, o) + 1.0d0
                end do

                ! Inversion via Cholesky decomposition, done in two stages.
                !
                call dpotrf('U', nlobs, X1, nlobs, lapack_info)
                if (lapack_info /= 0) then
                   print *, '  ERROR: m_local_analysis(): LAPACK error in dpotrf: errno = '&
                        , lapack_info, 'i, j =', i, j
                   call stop_mpi
                endif
             
                call dpotri('U', nlobs, X1, nlobs, lapack_info)
                if (lapack_info /= 0) then
                   print *, '  ERROR: m_local_analysis(): LAPACK error in dpotri: errno = '&
                        , lapack_info, 'i, j =', i, j
                   call stop_mpi
                endif
             
                ! fill the lower triangular part of (symmetric) X1
                !
                do o = 2, nlobs
                   X1(o, 1 :  o - 1) = X1(1 : o - 1, o)
                end do

                G = matmul(transpose(subS), X1)
             else ! nlobs >= nrens:  use ensemble space
                X1 = matmul(transpose(subS), subS)
                do m = 1, nrens
                   X1(m, m) = X1(m, m) + 1.0d0
                end do

                ! Inversion
                !
                call dpotrf('U', nrens, X1, nrens, lapack_info)
                if (lapack_info /= 0) then
                   print *, '  ERROR: m_local_analysis(): LAPACK error in dpotrf: errno = '&
                        , lapack_info, 'i, j =', i, j
                   call stop_mpi
                endif
                call dpotri('U', nrens, X1, nrens, lapack_info)
                if (lapack_info /= 0) then
                   print *, '  ERROR: m_local_analysis(): LAPACK error in dpotri: errno = '&
                        , lapack_info, 'i, j =', i, j
                   call stop_mpi
                endif
             
                do m = 2, nrens
                   X1(m, 1 :  m - 1) = X1(1 : m - 1, m)
                end do

                G = matmul(X1, transpose(subS))
             end if

             if (iter == 1) then
                 do m = 1, nrens
                   X5tmp(m, :) = sum(G(m, :) * subdy(:))
                end do
             end if

             10 continue

             ! calculate DFS at iteration 1, SRF at iteration 2
             !
             if (iter == 1) then
                dfs = traceprod(G, subS, nrens, nlobs)
                dfs_array(i, j) = real(dfs, 4)
                pnlobs = 0
                do uo = 1, nuobs
                   do o = 1, nlobs
                      if (lobs(o) >= uobs_begin(uo) .and.&
                           lobs(o) <= uobs_end(uo)) then
                         pnlobs(uo) = pnlobs(uo) + 1
                         plobs(pnlobs(uo), uo) = o
                      end if
                   end do
                end do
                pdfs = 0.0d0
                psrf = 0.0d0
                do uo = 1, nuobs
                   if (pnlobs(uo) > 0) then
                      pdfs(uo) = traceprod(G(:, plobs(1 : pnlobs(uo), uo)),&
                           subS(plobs(1 : pnlobs(uo), uo), :), nrens, pnlobs(uo))
                   end if
                   pdfs_array(i, j, uo) = real(pdfs(uo), 4)
                end do
             else
                if (dfs /= 0.0d0) then
                   srf = sqrt(traceprod(subS, transpose(subS), nlobs, nrens)&
                        / traceprod(G, subS, nrens, nlobs)) - 1.0d0
                else
                   srf = 0.0d0
                end if
                srf_array(i, j) = real(srf, 4)
                do uo = 1, nuobs
                   if (pnlobs(uo) > 0) then
                      if (pdfs(uo) /= 0.0d0) then
                         psrf(uo) = sqrt(&
                              traceprod(subS(plobs(1 : pnlobs(uo), uo), :),&
                              transpose(subS(plobs(1 : pnlobs(uo), uo), :)),&
                              pnlobs(uo), nrens) /&
                              traceprod(G(:, plobs(1 : pnlobs(uo), uo)),&
                              subS(plobs(1 : pnlobs(uo), uo), :),&
                              nrens, pnlobs(uo))) - 1.0d0
                      else
                         psrf(uo) = 0.0d0
                      end if
                   end if
                   psrf_array(i, j, uo) = real(psrf(uo), 4)
                end do
             end if
          end do ! iter

          if  (SCHEME_USED == SCHEME_ENKF) then
             X5tmp = X5tmp + matmul(G, subD - subS)
          elseif (SCHEME_USED == SCHEME_DENKF) then
             X5tmp = X5tmp - 0.5d0 * matmul(G, subS)
          end if
          do m = 1, nrens
             X5tmp(m, m) = X5tmp(m, m) + 1.0d0
          enddo

          if (testthiscell) then
             ! ensemble mean
             allocate(x(nlobs))
             do o = 1, nlobs
                x(o) = obs(lobs(o)) % d - dy(lobs(o)) * sqrtm * sqrt(obs(lobs(o)) % var)
             end do
             tmp(1) = rfactor
             call p2nc_writeobs(i, j, nlobs, nrens, X5tmp, modlon(i, j),&
                  modlat(i, j), depths(i, j), tmp(1), lobs(1 : nlobs), &
                  obs(lobs(1 : nlobs)), x, subS, subdy, lfactors)
             deallocate(x)
          end if

          ! Put X5tmp into the final X5 matrix - to be written to a file
          !
          X5(:, :, i) = real(X5tmp, 4)

          deallocate(subS, subdy, lfactors, X1, G)
          if  (SCHEME_USED == SCHEME_ENKF) then
             deallocate(subD)
          end if
       end do ! i = 1, ni

       ! Write one "stripe" of the temporary matrix X5 to disk
       iter = 0
       do while (.true.)
          iter = iter + 1
          write(17, rec = j, iostat = iostatus) X5
          if (iostatus /= 0) then 
             print *, 'ERROR: calc_X5(): I/O error at writing X5, iostatus = ',&
                  iostatus
             print *, 'ERROR: at model line j =', j, ' counter jj = ', jj, 'iter =', iter
             if (iter < MAXITER) then
                cycle
             else
                print *, 'ERROR: max number of iterations reached, STOP'
                stop
             end if
          end if
          read(17, rec = j, iostat = iostatus) X5check
          if (iostatus /= 0) then 
             print *, 'ERROR: calc_X5(): I/O error at reading X5, iostatus = ',&
                  iostatus
             print *, 'ERROR: at j = ', j, ' jj = ', jj, 'iter =', iter
             if (iter < MAXITER) then
                cycle
             else
                print *, 'ERROR: max number of iterations reached, STOP'
                stop
             end if
          end if
          if (abs(maxval(X5 - X5check)) > 1.0e-6) then
             print *, 'ERROR: calc_X5(): inconsistency between written/read X5'
             print *, 'ERROR: j = ', j, ' jj = ', jj, 'iter =', iter,&
                  ' maxval(X5 - X5check) =', maxval(X5 - X5check)
             if (iter < MAXITER) then
                cycle
             else
                print *, 'ERROR: max number of iterations reached, STOP'
                stop
             end if
          end if
          exit ! OK
       end do
       print *, 'FINISHED j =', j, ' jj =', jj
    end do ! j = my_first_iteration, my_last_iteration

    close(17) ! X5 file

    if (SCHEME_USED == SCHEME_ENKF) then
       deallocate(D)
    end if

#if defined(QMPI)
    if (.not. master) then
       ! broadcast nlobs and dfs arrays to master
       call send(nlobs_array(:, jmap(my_first_iteration : my_last_iteration)), 0, 0)
       call send(dfs_array(:, jmap(my_first_iteration : my_last_iteration)), 0, 1)
       call send(srf_array(:, jmap(my_first_iteration : my_last_iteration)), 0, 1)
       allocate(mpibuffer_float1(ni, my_last_iteration - my_first_iteration + 1))
       allocate(mpibuffer_float2(ni, my_last_iteration - my_first_iteration + 1))
       do uo = 1, nuobs
          mpibuffer_float1 = pdfs_array(:, jmap(my_first_iteration : my_last_iteration), uo)
          call send(mpibuffer_float1, 0, uo + 1)
          mpibuffer_float2 = psrf_array(:, jmap(my_first_iteration : my_last_iteration), uo)
          call send(mpibuffer_float2, 0, uo + 1)
       end do
       deallocate(mpibuffer_float1)
       deallocate(mpibuffer_float2)
    else
       ! receive nlobs and dfs arrays
       do p = 2, qmpi_num_proc
          !
          ! PS: Ideally, it would be nice to be able to use a simple code like:
          !
          ! call receive(nlobs_array(&
          !              jmap(first_iteration(p) : last_iteration(p))), p - 1)
          !
          ! but this seems not to work, at least with the PGI compiler. 
          ! Perhaps, this is too much to expect from a call to a C function...
          ! The good news is that using a temporal array works fine.
          !
          allocate(mpibuffer_int(ni, last_iteration(p) - first_iteration(p) + 1))
          call receive(mpibuffer_int, p - 1, 0)
          nlobs_array(:, jmap(first_iteration(p) : last_iteration(p))) = mpibuffer_int
          deallocate(mpibuffer_int)
          allocate(mpibuffer_float1(ni, last_iteration(p) - first_iteration(p) + 1))
          call receive(mpibuffer_float1, p - 1, 1)
          dfs_array(:, jmap(first_iteration(p) : last_iteration(p))) = mpibuffer_float1
          allocate(mpibuffer_float2(ni, last_iteration(p) - first_iteration(p) + 1))
          call receive(mpibuffer_float2, p - 1, 1)
          srf_array(:, jmap(first_iteration(p) : last_iteration(p))) = mpibuffer_float2
          do uo = 1, nuobs
             call receive(mpibuffer_float1, p - 1, uo + 1)
             pdfs_array(:, jmap(first_iteration(p) : last_iteration(p)), uo) = mpibuffer_float1
             call receive(mpibuffer_float2, p - 1, uo + 1)
             psrf_array(:, jmap(first_iteration(p) : last_iteration(p)), uo) = mpibuffer_float2
          end do
          deallocate(mpibuffer_float1)
          deallocate(mpibuffer_float2)
       enddo
    endif
    ! broadcast nlobs array
    call broadcast(nlobs_array)
#endif

    if (master) then
       nlobs_max = maxval(nlobs_array)
       print *, 'maximal # of local obs =', nlobs_max,&
            ' reached for', count(nlobs_array == nlobs_max), 'grid cells'
       print *, 'average #(*) of local obs =', sum(nlobs_array(:, 1 : nj)) / real(count(nlobs_array(:, 1 : nj) > 0))
       print *, '  * over cells with non-zero number of local obs only'
       print *, 'localisation function of type', LOCFUN_USED, 'has been used'
       print *, 'analysis conducted in obs space in', count(nlobs_array(:, 1 : nj) > 0 .and. nlobs_array(:, 1 : nj) < nrens),&
            'cells'
       print *, 'analysis conducted in ens space in', count(nlobs_array(:, 1 : nj) >= nrens),&
            'cells'
       print *, 'maximal DFS =', maxval(dfs_array)
       print *, 'average(*) DFS =', sum(dfs_array) / real(count(dfs_array > 0))
       print *, '  * over cells with non-zero number of local obs only'
       print *, '# of cells with DFS > N / 2 =', count(dfs_array > real(nrens / 2, 4))

       call diag2nc(ni, nj, modlon, modlat, nlobs_array, dfs_array, pdfs_array,&
            srf_array, psrf_array)
    end if
  end subroutine calc_X5


  integer function get_npad_la(ni, nj)
    integer, intent(in) :: ni, nj

    get_npad_la = 4096 - mod(ni * nj, 4096)
    get_npad_la = mod(get_npad_la, 4096)
  end function get_npad_la


  real function locfun(x)
    real, intent(in) :: x

    real :: xx, xx2, xx3

    select case(LOCFUN_USED)

    case (LOCFUN_NONE)
       locfun = 1.0
    case (LOCFUN_STEP)
       if (x > 1.0) then
          locfun = 0.0
       else
          locfun = 1.0
       end if
    case (LOCFUN_GASPARI_COHN)
       if (x > 1.0) then
          locfun = 0.0
       else
          xx = x * 2.0
          xx2 = xx * xx
          xx3 = xx2 * xx
          if (xx < 1.0) then
             locfun = 1.0 + xx2 * (- xx3 / 4.0 + xx2 / 2.0)&
                  + xx3 * (5.0 / 8.) - xx2 * (5.0 / 3.0)
          else
             locfun = xx2 * (xx3 / 12.0 - xx2 / 2.0)&
                  + xx3 * (5.0 / 8.0) + xx2 * (5.0 / 3.0)&
                  - xx * 5.0 + 4.0 - (2.0 / 3.0) / xx
          end if
          locfun = max(locfun, 0.0)
       end if
    case default
       print *, 'ERROR: m_local_analysis.F90: locfun(): LOCFUN_USED =', LOCFUN_USED, 'is unknown'
       stop
    end select
  end function locfun


  ! - Sort observations by their distance to the given grid point (i, j).  
  ! - Identify observations within a given radius `rmax'.
  ! - Select `nlobs' nearest observations; update `nlobs' if there are not
  ! enough observations within the radius.
  !
  ! Note that because all observations are parsed for each 2D grid point, this
  ! subroutine may become a bottleneck if the total number of observations
  ! grows substantially from the current point... If this happens, we may
  ! consider putting all observations in a K-D tree like in Szyonykh et. al
  ! (2008), A local ensemble transform Kalman filter data assimilation system
  ! for the NCEP global model (2008). Tellus 60A, 113-130.
  !
  subroutine get_local_obs(i, j, rmax0, modlon, modlat, mindx,&
       ni, nj, nlobs, lobs)
    use mod_measurement
    use m_obs
    use m_spherdist

    implicit none
    integer, intent(in) :: i, j
    real, intent(in) :: rmax0 ! maximal allowed distance
    real, intent(in) :: modlon(ni, nj)
    real, intent(in) :: modlat(ni, nj)
    real, intent(in) :: mindx
    integer, intent(in) :: ni, nj
    integer, intent(inout) :: nlobs ! input : max allowed # of local obs
                                   ! output: actual # of local obs for this
                                   !         point
    integer, intent(out) :: lobs(nobs) ! indices of local observations

    integer :: ngood
    integer :: sorted(nobs)
    real :: dist(nobs)
    integer :: o
    real :: rmax,rmax2

    lobs = 0
    ngood = 0
    do o = 1, nobs
       rmax=rmax0*Typobs(o)
       rmax2 = (rmax / mindx) ** 2
       if ((obs(o) % ipiv - i) ** 2 + (obs(o) % jpiv - j) ** 2 > rmax2) then
          cycle
       end if

       dist(o) = spherdist(obs(o) % lon, obs(o) % lat, modlon(i, j), modlat(i, j))
       if (dist(o) <= rmax) then
          ngood = ngood + 1
          lobs(ngood) = o
       end if
    end do

    if (nlobs <= 0 .or. nlobs >= ngood) then
       !
       ! use all observations within localisation support radius
       !
       nlobs = ngood
    else
       !
       ! use `nlobs' closest observations
       !
       call order(dble(nobs), dist, dble(ngood), lobs, sorted)
       lobs(1 : nlobs) = sorted(1 : nlobs)
    end if
  end subroutine get_local_obs


  ! This subroutine writes (1) the number of local observations, (2)
  ! the number of degrees of freedom of signal (DFS), and (3) spread reduction
  ! factor (SRF) to file "enkf_diag.nc"
  !
  subroutine diag2nc(ni, nj, lon, lat, nlobs_array, dfs_array, pdfs_array, &
       srf_array, psrf_array)
    use mod_measurement
    use m_obs
    use nfw_mod
    implicit none

    integer, intent(in) :: ni
    integer, intent(in) :: nj
    real, intent(in) :: lon(ni, nj)
    real, intent(in) :: lat(ni, nj)
    integer, intent(in) :: nlobs_array(ni, nj)
    real(4), intent(in) :: dfs_array(ni, nj)
    real(4), intent(in) :: pdfs_array(ni, nj, nuobs)
    real(4), intent(in) :: srf_array(ni, nj)
    real(4), intent(in) :: psrf_array(ni, nj, nuobs)

    character(STRLEN) :: fname
    character(STRLEN) :: varname
    integer :: ncid
    integer :: dimids(2)
    integer :: lon_id, lat_id, nlobs_id, dfs_id, pdfs_id(nuobs), srf_id,&
         psrf_id(nuobs)
    integer :: uo

    fname = 'enkf_diag.nc'
    call nfw_create(fname, nf_clobber, ncid)
    
    call nfw_def_dim(fname, ncid, 'i', ni, dimids(1))
    call nfw_def_dim(fname, ncid, 'j', nj, dimids(2))
    call nfw_def_var(fname, ncid, 'lon', nf_float, 2, dimids, lon_id)
    call nfw_def_var(fname, ncid, 'lat', nf_float, 2, dimids, lat_id)
    call nfw_def_var(fname, ncid, 'nobs', nf_int, 2, dimids, nlobs_id)
    call nfw_def_var(fname, ncid, 'dfs', nf_float, 2, dimids, dfs_id)
    do uo = 1, nuobs
       write(varname, '(a, a)') 'dfs_', trim(unique_obs(uo))
       call nfw_def_var(fname, ncid, trim(varname), nf_float, 2, dimids, pdfs_id(uo))
    end do
    call nfw_def_var(fname, ncid, 'srf', nf_float, 2, dimids, srf_id)
    do uo = 1, nuobs
       write(varname, '(a, a)') 'srf_', trim(unique_obs(uo))
       call nfw_def_var(fname, ncid, trim(varname), nf_float, 2, dimids, psrf_id(uo))
    end do
    call nfw_enddef(fname, ncid)

    call nfw_put_var_double(fname, ncid, lon_id, lon)
    call nfw_put_var_double(fname, ncid, lat_id, lat)
    call nfw_put_var_int(fname, ncid, nlobs_id, nlobs_array)
    call nfw_put_var_real(fname, ncid, dfs_id, dfs_array)
    call nfw_put_var_real(fname, ncid, srf_id, srf_array)
    do uo = 1, nuobs
       call nfw_put_var_real(fname, ncid, pdfs_id(uo), pdfs_array(:, :, uo))
       call nfw_put_var_real(fname, ncid, psrf_id(uo), psrf_array(:, :, uo))
    end do
    call nfw_close(fname, ncid)
  end subroutine diag2nc


  ! Calculates the trace of a product of two matrices. (Does not calculate
  ! the off-diagonal elements in the process.)
  !
  real function traceprod(A, B, n, m)
    real, intent(in) :: A(n, m), B(m, n)
    integer, intent(in) :: n, m

    integer :: i

    traceprod = 0.0d0
    do i = 1, n
       traceprod = traceprod + sum(A(i, :) * B(:, i))
    end do
  end function traceprod

end module m_local_analysis
