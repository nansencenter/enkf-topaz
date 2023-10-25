module mod_mush_ktherm
   use mod_cice_constants 

contains
   ! calculation the frezing temperature which relates with the model setting
   real function Tfrz(sss,ktherm) 
      use mod_cice_constants
      implicit none
      integer,    intent(in) :: ktherm
      real,       intent(in) :: sss
      if(ktherm==1)  then        ! BL99
        Tfrz=sss*(-0.054)
      elseif(ktherm==2) then     ! therm_mushy 
        Tfrz=sss/(0.01849*sss-18.48)
      end if
   end function 


! prescribe vertical profile of salinity and melting temperature
! and then to calcualte for mush layers in cice 
subroutine ice_vertical_layer(zqin, zSin,Nilayer,sss,Tfrz,Tsf,ktherm)
   implicit none
 
   ! refer to ice_therm_itd.F90 used by mushy layer
   real,                     parameter :: dSin0_frz=3.0 ! bulk salinity reduction of newly formed frazil
   real,                     parameter :: phi_ini=0.75  ! initial frazil liquid fraction

   integer,                 intent(in) :: ktherm
   integer,                 intent(in) :: Nilayer
   real,                    intent(in) :: sss
   real,                    intent(in) :: Tfrz    ! freezing temperature(C)
   real,                    intent(in) :: Tsf    ! ice surface temperature, Tsfcn
   real,dimension(Nilayer),intent(out) :: zqin,zSin
   
   real                    :: slope,Ti,Tmlt,Sicez
   real                    :: zn
   integer                 :: k

   real                    :: Si0,Ti0
   real                    :: Qi0
   slope=Tfrz-Tsf
   if (ktherm==1) then  ! only used for BL99 thermodynamics
     do k=1,Nilayer
        zn=(real(k)-0.5)/real(Nilayer)
        Sicez=(saltmax/2.)*(1-cos(pi*zn**(nsal/msal+zn)))
        Sicez=max(Sicez,min_salin)
        Tmlt=-Sicez*depressT
        Ti=min(Tsf+slope*zn,Tsf+(Tmlt-Tsf)*zn)
        !l_brine         ! if true, treat brine pocket effects
        zqin(k) = -rhoi*(cp_ice*(Tmlt-Ti) + &
                   (Lsub-Lvap)*(1.0-Tmlt/Ti)-cp_ocn*Tmlt)
        zSin(k)=Sicez
     end do

   elseif (ktherm==2) then  ! used for mushy referring to ice_therm_itd.F90
     if(sss>2.0*dSin0_frz) then
        Si0=sss-dSin0_frz
     else
        Si0=sss**2/(4.0*dSin0_frz)
     end if
     ! liquidus_temperature_mush
     call liquidus_temperature_mush(Si0/phi_ini,Ti0)
     Ti0=min(Ti0,-0.1)
     zSin=Si0
     call enthalpy_mush(Ti0,Si0,Qi0)
     zqin=Qi0

   else
     print '(a30,i5)', 'Wrong input for ktherm=', ktherm
     stop
   endif 

end subroutine 


subroutine liquidus_temperature_mush(Sbr, zTin)
    ! liquidus relation: equilibrium temperature as function of brine salinity
    ! based on empirical data from Assur (1958)

    real,  parameter :: Sb_liq =  123.667028
    real,  parameter :: &
       az1_liq = -18.48         , &
       bz1_liq = 0.0            , &
       az2_liq = -10.3085       , &
       bz2_liq = 62.4           
    real,  intent(in) :: Sbr    ! ice brine salinity (ppt)
    real, intent(out) :: zTin   ! ice layer temperature (C)

    real :: t_high ! mask for high temperature liquidus region
    real :: M1_liq, M2_liq, N1_liq, N2_liq, O1_liq, O2_liq

    M1_liq = az1_liq            
    N1_liq = -az1_liq/1000.0         
    O1_liq = -bz1_liq/az1_liq 
    M2_liq = az2_liq            
    N2_liq = -az2_liq/1000.0          
    O2_liq = -bz2_liq / az2_liq
    t_high = merge(1.0, 0.0, (Sbr <= Sb_liq))

    zTin = ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + &
          ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (1.0 - t_high)

end subroutine liquidus_temperature_mush




!%  doubl check ice_therm_mushy.F90:
subroutine enthalpy_mush(zTin, zSin, zqin)
    ! enthalpy of mush from mush temperature and bulk salinity
    real,  intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real, intent(out) :: &
         zqin    ! ice layer enthalpy (J m-3) 

    real :: &
         phi     ! ice liquid fraction 

    phi = liquid_fraction(zTin, zSin)
    
    zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
           rhoi * cp_ice * zTin - (c1 - phi) * rhoi * Lfresh

end subroutine enthalpy_mush




!=======================================================================
! Mushy Layer Formulation - Assur (1958) liquidus
!=======================================================================

  function liquidus_brine_salinity_mush(zTin) result(Sbr)

    ! liquidus relation: equilibrium brine salinity as function of temperature
    ! based on empirical data from Assur (1958)

    real, intent(in) :: &
         zTin         ! ice layer temperature (C)

    real :: &
         Sbr          ! ice brine salinity (ppt)

    real :: &
         t_high   , & ! mask for high temperature liquidus region
         lsubzero     ! mask for sub-zero temperatures

    t_high   = merge(c1, c0, (zTin > Tb_liq))
    lsubzero = merge(c1, c0, (zTin <= c0))

    Sbr = ((zTin + J1_liq) / (K1_liq * zTin + L1_liq)) * t_high + &
          ((zTin + J2_liq) / (K2_liq * zTin + L2_liq)) * (c1 - t_high)

    Sbr = Sbr * lsubzero

  end function liquidus_brine_salinity_mush


function liquid_fraction(zTin, zSin) result(phi)

    ! liquid fraction of mush from mush temperature and bulk salinity

    real, intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real :: &
         phi , & ! liquid fraction
         Sbr     ! brine salinity (ppt)

    Sbr = max(liquidus_brine_salinity_mush(zTin),puny)
    phi = zSin / max(Sbr, zSin)

  end function liquid_fraction





end module
