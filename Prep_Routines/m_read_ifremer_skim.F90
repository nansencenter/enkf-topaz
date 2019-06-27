
module m_read_ifremer_skim
  implicit none

  integer, parameter, private :: STRLEN =2180
  integer,parameter, private :: Msamp0=5000

  public read_ifremer_skim
  
  private data_inquire
  private data_readfile
  private grid_readxyz
  private land_nearby 


contains

  subroutine read_ifremer_skim(fnames,obstype,variance, nx,ny,data, datainfo)
    use mod_measurement
    use m_oldtonew
    use m_confmap
    use m_bilincoeff
    use m_pivotp
    use nfw_mod
    character(*), intent(in) :: fnames
    character(*), intent(in) :: obstype
    real, intent(in) :: variance
    integer, intent(in) :: nx, ny
    type(measurement), allocatable, intent(out) :: data(:)
    type(measurement_info), allocatable, intent(out),optional :: datainfo(:)

    character(STRLEN) :: fname
    integer :: nfile, nbeam, nsamp

    real(8), allocatable :: juld(:,:),lon(:,:), lat(:, :)
    real(8), allocatable :: radia_an(:,:),beam_an(:,:)
    integer, allocatable :: ipiv(:,:), jpiv(:,:), mask(:,:)

    real(8), dimension(nx, ny) :: modlat, modlon
    real(8), dimension(nx, ny) :: depths

    integer :: ndata 
    integer :: f, l, p, np

    real :: latnew, lonnew
    real :: wetsill
    integer :: imin,imax,jmin,jmax
    logical :: ex1
    real :: rvar


    print *, 'BEGIN read_ifremer_skim(): ',trim(fnames)
    call data_inquire(fnames, nfile, nbeam, nsamp)
    print *, '  overall: nbeam =', nbeam, ', nsamp =', nsamp,' nfile=', nfile
    if ( nsamp > Msamp0 ) then
       nsamp = Msamp0
       print *, '  Limited overall: nbeam =', nbeam, ', nsamp =', nsamp
    endif


    allocate(juld(nsamp,nbeam))
    allocate(lon(nsamp,nbeam))
    allocate(lat(nsamp,nbeam))
    allocate(radia_an(nsamp,nbeam))
    allocate(beam_an(nsamp,nbeam))
    radia_an=9999.
    beam_an=9999.
    juld=9999.
    lat=9999.
    lon=9999.
    p=1
    do f = 1, nfile
      call data_readfile(f, trim(obstype), np, juld(1:nsamp,p : nbeam),&
           lon(1:nsamp,p : nbeam), lat(1:nsamp,p : nbeam),&
           radia_an(1 : nsamp,p:nbeam), beam_an(1:nsamp,p : nbeam))
      p = p + np
    end do


    allocate(ipiv(nsamp,nbeam))
    allocate(jpiv(nsamp,nbeam))
    allocate(mask(nsamp,nbeam))
    ipiv=0;  jpiv=0;  mask=1;

    where (lat > 90) mask=0
    print *, '    ', count(mask == 1), ' initial locations.'
    call confmap_init(nx, ny)
    do p = 1, nbeam
      do l =1, nsamp
        if (mask(l,p)==1) then
          call oldtonew(real(lat(l,p)),real(lon(l,p)),latnew,lonnew)
          call pivotp(lonnew, latnew, ipiv(l,p), jpiv(l,p))
        endif
      end do
    end do
    where (ipiv < 2 .or. jpiv < 2 .or. ipiv > nx - 1 .or. jpiv > ny - 1) mask = 0
    print *, '  after calculaling pivot points in model domain:'
    print *, '    ', count(mask == 1), ' valid locations'


    ! Check for the observation being wet
    ! deeper than 200 m and at least 50 km far from coastline
    wetsill=200
    call grid_readxyz(nx, ny, modlat, modlon, depths)
    do p = 1, nbeam
      do l = 1, nsamp
        if (mask(l,p) == 1) then
          imin= max(ipiv(l,p)-1,1)
          imax= min(ipiv(l,p)+1,nx)
          jmin= max(jpiv(l,p)-1,1)
          jmax= min(jpiv(l,p)+1,ny)
          if (any(depths(imin:imax,jmin:jmax)< wetsill.or. &
depths(imin:imax,jmin:jmax) == depths(imin:imax, jmin:jmax) &
       + 1.0)) then
            mask(l,p)=0
          else
            ex1=land_nearby(nx, ny, depths, modlon, modlat,&
               ipiv(l,p), jpiv(l,p), lon(l,p), lat(l,p),50000.0)
            if (ex1==.true.) mask(l,p)=0
          endif
        endif
      end do
    end do
    print *, '  after removing locations on land:'
    print *, '    ', count(mask == 1), ' valid locations'


    ndata = 0
    allocate(data(count(mask==1)))
    do p = 1, nbeam
      do l =1, nsamp
        if (mask(l,p) == 1) then
          ndata = ndata + 1
         ! call random_seed()
          call random_number(rvar)
          data(ndata) % d =2*rvar
          if (beam_an(l,p)<=6) then
            data(ndata) % var = 0.04 
          elseif (beam_an(l,p)>=12) then
            data(ndata) % var = 0.17**2 
          endif
          data(ndata) % id = trim(obstype)
          data(ndata) % lon = real(lon(l,p))
          data(ndata) % lat = real(lat(l,p))
          data(ndata) % ipiv = ipiv(l,p)
          data(ndata) % jpiv = jpiv(l,p)
          data(ndata) % ns   = 0
          data(ndata) % date = 0
          data(ndata) % depth = 1 

          call bilincoeff(real(modlon), real(modlat), nx, ny, &
               real(lon(l,p)), real(lat(l,p)), ipiv(l,p), jpiv(l,p), &
               data(ndata) % a1, data(ndata) % a2, data(ndata) % a3,&
               data(ndata) % a4)
!
          data(ndata) % status = .true. ! (active)
          data(ndata) % i_orig_grid = l 
          data(ndata) % j_orig_grid = p 

          if (mod(ndata,10000)==0) print *, ndata
        endif
      end do
    end do


   deallocate(juld)
   deallocate(lon)
   deallocate(lat)
   deallocate(radia_an)
   deallocate(beam_an)
   deallocate(ipiv,jpiv,mask)

  end subroutine

  logical function land_nearby(nx, ny, depths, modlon, modlat, ipiv, jpiv, obslon, obslat,Vdis)
    use m_spherdist
    implicit none
    integer, intent (in) :: nx, ny, ipiv, jpiv
    real*8, dimension(nx,ny), intent(in) :: depths, modlon, modlat
    real, intent (in) :: obslon,obslat 
    real, intent (in) :: Vdis    ! ~ meters minimum far from coastline

    integer :: ii, jj, ncells
    real :: griddist

    land_nearby = .false.
    ncells = ceiling(Vdis / spherdist(modlon(ipiv, jpiv), modlat(ipiv, jpiv),&
         modlon(ipiv, jpiv + 1), modlat(ipiv, jpiv + 1)))
    do jj = max(jpiv - ncells, 1), min(jpiv + ncells, ny)
       do ii = max(ipiv - ncells, 1), min(ipiv + ncells, nx)
          griddist = spherdist(modlon(ii, jj), modlat(ii, jj), obslon, obslat)
          if ((depths(ii,jj) < 1.or.depths(ii,jj)>20000) .and. griddist < Vdis) then
            ! print *, 'land_nearby: ',ipiv,jpiv
          !   print *,modlon(ii, jj), modlat(ii, jj), obslon, obslat
          !   print *,griddist,depths(ii,jj)
             land_nearby = .true.
             return
          end if
       enddo
    enddo
  end function land_nearby



  subroutine data_readfile(fid, obstype, mbeam, juld_all, lon_all, & 
    lat_all, d1_all, d2_all)
    use nfw_mod
    implicit none
    integer, intent(in) :: fid
    character(*), intent(in) :: obstype
    integer, intent(inout) :: mbeam
    real(8), intent(inout), dimension(:,:) :: juld_all
    real(8), intent(inout), dimension(:,:) :: lon_all, lat_all
    real(8), intent(inout), dimension(:,:) :: d1_all,d2_all

    character*80 :: fname
    integer :: f
    integer :: ncid
    integer :: id
    integer :: msamp,m_samp

    real, dimension(:,:),allocatable :: tmp_all


    open(10, file = 'infiles.txt')
    do f = 1, fid
       read(10, fmt = '(a)') fname
    end do
    close(10)
    !print *, '  reading "', trim(fname), '"'

    call nfw_open(trim(fname), nf_nowrite, ncid)

    ! nbeam
    !
    call nfw_inq_dimid(fname, ncid, 'nbeam', id)
    call nfw_inq_dimlen(fname, ncid, id, mbeam)
    !print *, '    nbeam = ', mbeam

    ! nsamp
    !
    call nfw_inq_dimid(fname, ncid, 'sample', id)
    call nfw_inq_dimlen(fname, ncid, id, msamp)

    allocate(tmp_all(msamp,mbeam))
    m_samp=min(Msamp0,msamp);

    !juld 
    !
    call nfw_inq_varid(fname, ncid, 'time', id)
    call nfw_get_var_double(fname, ncid, id, tmp_all(1:msamp,1:mbeam))
    do f = 1,mbeam 
      juld_all(1:m_samp,f)=tmp_all(1:m_samp,f)
    end do
    !lon 
    !
    call nfw_inq_varid(fname, ncid, 'lon', id)
    call nfw_get_var_double(fname, ncid, id, tmp_all(1:msamp,1:mbeam))
    do f = 1, mbeam
      lon_all(1:m_samp,f)=tmp_all(1:m_samp,f)
    end do
    !lat 
    !
    call nfw_inq_varid(fname, ncid, 'lat', id)
    call nfw_get_var_double(fname, ncid, id, tmp_all(1:msamp,1:mbeam))
    do f = 1, mbeam
      lat_all(1:m_samp,f)=tmp_all(1:m_samp,f)
    end do
    !radial_angle 
    !
    call nfw_inq_varid(fname, ncid, 'radial_angle', id)
    call nfw_get_var_double(fname, ncid, id, tmp_all(1:msamp,1:mbeam))
    do f = 1, mbeam
      d1_all(1:m_samp,f)=tmp_all(1:m_samp,f)
    end do
    !beam_angle 
    !
    call nfw_inq_varid(fname, ncid, 'beam_angle', id)
    call nfw_get_var_double(fname, ncid, id, tmp_all(1:msamp,1:mbeam))
    do f = 1, mbeam
      d2_all(1:m_samp,f)=tmp_all(1:m_samp,f)
    end do

!    print *,'ok: tmp_all(1,1)', &
!juld_all(1,1),lon_all(1,1),lat_all(1,1),d1_all(1,1),d2_all(1,1)
    

  end subroutine





  subroutine data_inquire(fnames,nfile,nbeam,nsamp)
    use nfw_mod
    character(*), intent(in) :: fnames
    integer, intent(inout) :: nfile, nbeam, nsamp

    character(STRLEN) :: command ! (there may be a limit of 80 on some systems)
    character(STRLEN) :: fname
    integer :: ios
    integer :: ncid
    integer :: id

    integer :: nbeam_this, nsamp_this
    integer :: ns1,ns2,samp_this

    nfile = 0
    nbeam = 0
    nsamp = 0

    call system("ls "//trim(fnames)//" > infiles.txt");

    nfile = 0
    open(10, file = 'infiles.txt')
    do while (.true.)
       read(10, fmt = '(a)', iostat = ios) fname
       if (ios /= 0) then
          exit
       end if

       nfile = nfile + 1
       print *, '  file #', nfile, ' = "', trim(fname), '"'

       call nfw_open(trim(fname), nf_nowrite, ncid)

       ! nbeam
       !
       call nfw_inq_dimid(fname, ncid, 'nbeam', id)
       call nfw_inq_dimlen(fname, ncid, id, nbeam_this)
       print *, '    nbeam = ', nbeam_this

       ! nsamp
       !
       call nfw_inq_dimid(fname, ncid, 'sample', id)
       call nfw_inq_dimlen(fname, ncid, id, nsamp_this)
       
       nbeam = nbeam + nbeam_this
       if (nsamp_this > nsamp) then
          nsamp = nsamp_this
          print *, 'Increasing .. ',trim(fname)
       end if

       call nfw_close(fname, ncid)
    end do
    close(10)

  end subroutine


  subroutine grid_readxyz(nx, ny, lat, lon, depth)
    integer, intent(in) :: nx, ny
    real(8), dimension(nx, ny), intent(inout) :: lat, lon, depth

    logical :: exists
    character(len = 128) :: fname
    
    fname = 'newpos.uf'
    inquire(file = fname, exist = exists)
    if (.not. exists) then
       print *, 'grid_readxyz(): ERROR: "', trim(fname), '" does not exist'
       stop
    end if
    open(10, file = fname, form = 'unformatted', status = 'old')
    print *, '  grid_readxyz(): reading "', trim(fname), '"...'
    read(10) lat, lon
    close(10)

    write(fname, '(a, i3.3, a, i3.3, a)') 'depths', nx, 'x', ny, '.uf'
    inquire(file = fname, exist = exists)
    if (.not. exists) then
       print*, 'grid_readxyz(): ERROR: "', trim(fname), '" does not exist'
       stop
    end if
    open (unit = 10, file = fname, status = 'old', form = 'unformatted')
    print *, '  grid_readxyz(): reading "', trim(fname), '"...'
    read(10) depth
    close(10)
  end subroutine grid_readxyz




end module m_read_ifremer_skim
