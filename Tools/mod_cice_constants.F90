module mod_cice_constants
   integer,  parameter :: Nilayer  = 7   ! vertical layers for ice
   integer,  parameter :: Nslayer  = 1   ! vertical layers for snow

   real,     parameter :: rhos     = 330.0   ! density of snow (kg/m^3)
   real,     parameter :: rhoi     = 917.0   ! density of ice (kg/m^3)
   real,     parameter :: rhow     = 1026.0  ! density of seawater (kg/m^3)
   real,     parameter :: cp_air   = 1005.0  ! specific heat of air (J/kg/K)
   ! (Briegleb JGR 97 11475-11485  July 1992)
   real,     parameter :: cp_ice   = 2106.0   ! specific heat of fresh ice (J/kg/K)
   real,     parameter :: cp_ocn   = 4218.0   ! specific heat of ocean (J/kg/K)
   real,     parameter :: Lsub     = 2.835e6  ! latent heat, sublimation freshwater (J/kg)
   real,     parameter :: Lvap     = 2.501e6  ! latent heat, vaporization freshwater (J/kg)
                                              
   real,     parameter :: Lfresh   = 0.334e6  ! latent heat of melting fresh ice (j/kg)
                                              ! =Lsub-Lvap
   real,     parameter :: pi       = 3.1415926
   real,     parameter :: depressT = 0.054
   ! ice_constants.F90:
   real,     parameter :: Tocnfrz  = -1.8    ! freezing temp of seawater (C),
   real,     parameter :: Tffresh  = 273.15  ! freezing temp of fresh ice (K)
   real,     parameter :: Tsmelt   = 0.0     ! melting temperature, snow top surface

   real,     parameter :: saltmax  = 3.2     ! max salt in ice
   real,     parameter :: min_salin=0.1
   real,     parameter :: hs_min   = 1.e-4     ! min snow thickness for computing zTsn (m)

   real,     parameter :: msal     = 0.573 
   real,     parameter :: nsal     = 0.407 

   type iced_varname
      character(len=:), allocatable :: var
   end type iced_varname

  ! three types of variables in iced should be updated 
   integer,  parameter :: Nv3=3, Nv22=18, Nv32=25
   type(iced_varname)  :: Var3d1(Nv3),Var3d2(Nv32),Var2d2(Nv22)

contains


! define the ice-thickness category as reference to ice_itd.F90
subroutine hi_cate(hilim,ncat,kcatb,kitd)
   implicit none
   integer, intent(in) :: ncat, kcatb,kitd
   real,   intent(out) :: hilim(ncat+1)
   real    :: cc1, cc2, cc3
   real    :: x1,h_min   ! minmum ice thickness for thermo
   integer :: k
   h_min=0.01;
   if (kcatb==0) then
      if (kitd==1) then  ! linear remapping itd category limits
        cc1=3.0/real(ncat) 
        cc2=15.0*cc1 
        cc3=3.0
        hilim(1)=0.0    ! minimum ice thickness, m
      else
       ! delta function itd category limits
        h_min = 0.1    ! minimum ice thickness allowed (m) for therm
        cc1=max(1.1/real(ncat),h_min*1.0)
        cc2=25*cc1
        cc3=2.25
        hilim(1)=h_min
      end if
      do k = 2, ncat+1
        x1=real(k-2)/real(ncat) 
        hilim(k)=hilim(k-1)+cc1+cc2*(1.0+tanh(cc3*(x1-1.0)))
      end do
   elseif (kcatb==1) then
      cc1=3.0/real(ncat)
      cc2=0.5/real(ncat)
      hilim(1)=0.0
      do k = 2, ncat+1
         x1 = real(k-1)
         hilim(k) = x1 * (cc1+(k-2)*cc2)
      end do
  elseif (kcatb==2)  then
      if (ncat==5) then
         ! WMO5 thinnest 3 categories combined
         hilim(1)=0.0; hilim(2)=0.3;
         hilim(3)=0.7; hilim(4)=1.2;
         hilim(5)=2.0; hilim(6)=999.;
      elseif (ncat==6) then
         ! WMO5 thinnest 2 categories combined
         hilim(1)=0.0; hilim(2)=0.15;
         hilim(3)=0.3; hilim(4)=0.7;
         hilim(5)=1.2; hilim(6)=2.0;
         hilim(7)=999.;
      elseif (ncat==7) then
         ! all thickness categories 
         hilim(1)=0.0; hilim(2)=0.1;
         hilim(3)=0.15; hilim(4)=0.3;
         hilim(5)=0.7; hilim(6)=1.2;
         hilim(7)=2.0; hilim(8)=999.;
      else
         print *, 'kcatbound=2 (WMO) must have ncat=5, 6 or 7'
         stop
      endif 
   elseif (kcatb==-1)  then   ! single category
     hilim(1)=0.0;  hilim(2)=100.0
   else
     print *, 'Error: no support by current version!'
     print *, 'kcatbound=1 (-1, 0 or 2) '
     stop
   end if
end subroutine hi_cate

! reference to ice_init.F90
! initializing aicen at new ice grid
! hbar=3.0 using 5 categories (Mean SIT 2.52 m)also used by cice 5.1
! hbar=2.0 using 4 categories (mean SIT 1.73 m)
! hbar=1.6 using 3 categories (mean SIT 1.30 m)
subroutine Initial_aicen(Ncat,aicen,hilim,hbar)
   implicit none
   integer, intent(in) :: Ncat
   real,    intent(in) :: hbar,hilim(Ncat+1)
   real,   intent(out) :: aicen(Ncat)
  
   integer  :: i
   real     :: hinit,suma0

   suma0=0
   do i=1,Ncat
      if (i<Ncat) then
         hinit=0.5*(hilim(i)+hilim(i+1))
      else
         hinit=hilim(i)+1.0
      endif
      aicen(i)=max(0.,2.*hbar*hinit-hinit**2)
   end do
   suma0=sum(aicen(:))
   do i=1,Ncat
      aicen(i)=aicen(i)/suma0
   end do
end subroutine Initial_aicen


! fixed salinity profile and melting tempertuare in ice layers
! only working under ktherm=1 for BL99 thermo.
! Nlay=Nilayer+1
subroutine fix_zsinprofile(zSin,Tmlt,Nlay)
   implicit none
   integer,   parameter :: ktherm=1   ! for BL99
   integer,  intent(in) :: Nlay
   real,    intent(out) :: zSin(Nlay), Tmlt(Nlay)
   real    :: z
   integer :: k,Nlay0
   Nlay0=Nlay-1
   do k=1,Nlay-1
      z=(k-.5)/Nlay0
      zSin(k)=0.5*saltmax*(1-cos(pi*z**(nsal/(z+msal))))
      Tmlt(k)=-depressT*zSin(k)
      !print *, 'k=',k,z,zSin(k)
   end do
   zSin(Nlay)=saltmax
   Tmlt(Nlay)=-depressT*saltmax
   !print *, 'k=',k,' 1 ',zSin(k)
end subroutine fix_zsinprofile



end module
