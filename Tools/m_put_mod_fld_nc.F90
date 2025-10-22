module m_put_mod_fld_nc
! JX-This routine has two main functions:
! firstly postprocess the concerned variables in iced file for CICE
! secondly ouput the final fields into the CICE restart file
   use nfw_mod
   use mod_cice_constants
   use m_get_mod_fld_nc
!  emissivity= 0.95_dbl_kind    ,&! emissivity of snow and ice

   ! 'minus1p8'         Tf = -1.8 C (default)
   ! ice_forcing.F90: ! 'linear_salt'   Tf = -depressT * sss
   ! ice_forcing.F90: ! 'mushy'         Tf conforms with mushy layer thermo (ktherm=2)

contains

!   add parameter :: Vadjust = 0        ! skip to adjust the ice/snow volume
!  Used for BL99 in CICE
subroutine fix_cice(fice,hice,sss,nx,ny,ncat,restart,icerestart,yearday,Ahice,Vadjust)
   use mod_cice_constants
   use mod_mush_ktherm   
   implicit none
   integer,                    parameter :: ktherm   = 1   ! Note: match with the cice model
   integer,                   intent(in) :: nx,ny,ncat
   real, dimension(nx,ny), intent(inout) :: fice,hice     ! mandatory fields in analysisfields_ice.in
   real, dimension(nx,ny),    intent(in) :: sss  
   character(*),              intent(in) :: restart, icerestart  ! old restarts
   integer,                   intent(in) :: yearday, Vadjust
   real,                      intent(in) :: Ahice

   !real,     parameter :: spalv0=0.01
   real,     parameter :: thresh=1.0e-4         ! threshold for enthalpy update
   real,     parameter :: Athresh=5.0e-4        ! threshold total contrentration 
   real,     dimension(nx,ny,ncat) :: aicen,vicen,vsnon  

   ! internal variables in subtouine
   integer                      :: Levlayer 
   integer                      :: i,j,k
   real                         :: Tair0,trcrn_i,trcrn_s
   real                         :: Tmz,Ti
   real                         :: tmp0, tmp1, htmp
   logical                      :: Tupdate
   
   real, dimension(nx,ny)       :: ficem,hicem,hsnwm
   real, dimension (nx,ny)      :: fld2d,fld2d1,fld2d2 
   integer, dimension (nx,ny)   :: Imask 
   integer                      :: tlevel,ll
   character(len=20)            :: cfld
   
   real, dimension (nx,ny,ncat) :: fld3d
   character(len=80)            :: bfix_restart, new_icerestart

   real, allocatable, dimension(:,:,:) :: aicen_f,vicen_f,vsnon_f 
   real, allocatable, dimension(:,:,:) :: qicen_f,qsnon_f 
   real, allocatable,   dimension(:,:) :: fice_f,hice_f 
   real, allocatable,   dimension(:,:) :: fld2d_f,fld2d_Tf
   real, allocatable,   dimension(:,:) :: Tf 
   real, allocatable,     dimension(:) :: bkcat,bkaice      ! thickness threshholds
   real, allocatable,     dimension(:) :: ModzSin,ModTmlt

   ! for debugging
   integer,parameter  ::  ist=510,jst=580
   integer            ::  ibor,jbor
   real, allocatable,   dimension(:,:) :: Vtemp ! weight adjusted by HICE

   ! active 3D variables
   Var3d1(1)%var='aicen';    Var3d1(2)%var='vicen'; Var3d1(3)%var='vsnon'  

   ! other 3D variables delimited by aicen/vicen/vsnon
   Var3d2(1)%var='iage';  Var3d2(2)%var='FY';      Var3d2(3)%var='alvl';
   Var3d2(4)%var='vlvl';   Var3d2(5)%var='apnd';    Var3d2(6)%var='hpnd';
   Var3d2(7)%var='ipnd';    Var3d2(8)%var='dhs';     Var3d2(9)%var='ffrac';
   Var3d2(10)%var='Tsfcn';
   Var3d2(11)%var='sice001'; Var3d2(12)%var='sice002';Var3d2(13)%var='sice003';
   Var3d2(14)%var='sice004';  Var3d2(15)%var='sice005';Var3d2(16)%var='sice006';
   Var3d2(17)%var='sice007';  Var3d2(18)%var='qice001'; Var3d2(19)%var='qice002';
   Var3d2(20)%var='qice003'; Var3d2(21)%var='qice004'; Var3d2(22)%var='qice005';
   Var3d2(23)%var='qice006';Var3d2(24)%var='qice007'; Var3d2(25)%var='qsno001';

   ! other 2D variables delimited by ficem/hicem/hsnwm
   Var2d2(1)%var='uvel';    Var2d2(2)%var='vvel';      Var2d2(3)%var='strocnxT';
   Var2d2(4)%var='strocnyT'; Var2d2(5)%var='stressp_1'; Var2d2(6)%var='stressp_2';
   Var2d2(7)%var='stressp_3'; Var2d2(8)%var='stressp_4'; Var2d2(9)%var='stressm_1';
   Var2d2(10)%var='stressm_2'; Var2d2(11)%var='stressm_3'; Var2d2(12)%var='stressm_4';
   Var2d2(13)%var='stress12_1'; Var2d2(14)%var='stress12_2';Var2d2(15)%var='stress12_3';
   Var2d2(16)%var='stress12_4';  Var2d2(17)%var='iceumask'; Var2d2(18)%var='frz_onset'


   bfix_restart='bfix_'//trim(restart)
   call system('cp '//trim(restart)//' '//trim(bfix_restart));

   new_icerestart='fix_'//trim(icerestart)
   call system('cp '//trim(icerestart)//' '//trim(new_icerestart));

   Levlayer=Nilayer+1

   Imask=0   
   where(sss>0.and.sss<99)
      Imask=1         ! sea water point
   end where

   ! bakup the analyzed ice statment for debugging
   call replace_var_ncfile(fice,nx,ny,1,'ficem',bfix_restart,2,Imask)
   call replace_var_ncfile(hice,nx,ny,1,'hicem',bfix_restart,2,Imask)

   select case(yearday)
      case (60:150)
          Tair0=263.15
      case (151:330)
          Tair0=273.0
      case default
          Tair0=253.0
   end select
   trcrn_s=min(0.,Tair0-273.15)
   trcrn_i=0.

   ! read the background of aicen/vicen
   allocate(aicen_f(nx,ny,ncat),vicen_f(nx,ny,ncat),vsnon_f(nx,ny,ncat))
   allocate(qicen_f(nx,ny,ncat),qsnon_f(nx,ny,ncat))
   allocate(fice_f(nx,ny),hice_f(nx,ny))
   allocate(fld2d_f(nx,ny),fld2d_Tf(nx,ny))
   allocate(ModzSin(Levlayer),ModTmlt(Levlayer),bkcat(Levlayer)) 
   allocate(Tf(nx,ny),bkaice(ncat)) 
   allocate(Vtemp(nx,ny)) 

   Vtemp=1.0
   ! refer to ice_in
   call hi_cate(bkcat,ncat,0,1)
   call Initial_aicen(ncat,bkaice,bkcat,2.0)  ! first four categories for new ice as default

   !Ocean freezing temperature 
   Tf=0.
   do j=1,ny
     do i=1,nx
       if (Imask(i,j)==1) then  ! filter the land point
          if (ktherm==1) then
             Tf(i,j)=-depressT*sss(i,j)   ! It can be tuned if not stable
          elseif (ktherm==2) then
             Tf(i,j)=sss(i,j)/(-18.48+(18.48*sss(i,j)/1000.))
          else
             Tf(i,j)=-1.8
          endif
       endif
     end do
   end do

   !subroutine fix_zsinprofile(zSin,Tmlt,Nlay)
   call fix_zsinprofile(ModzSin,ModTmlt,Levlayer)
   !print '(a15,8f8.4)', 'Kcate bound SIT:',bkcat(1:ncat+1)
   !print '(a15,8f8.4)', '    Ice melt T.:',ModTmlt(:)

   ! for some case to keep the same ice statement
   if (Vadjust==0) then
      print *, 'Turn back with the previous ice statement!'
      ! file writing test:
      tlevel=1
      do k=1,ncat
         call get_mod_fld_nc(trim(icerestart), fld2d,'qice001', k, tlevel, nx, ny)
         aicen_f(:,:,k)=fld2d;
         call get_mod_fld_nc(trim(icerestart), fld2d,'qice002', k, tlevel, nx, ny)
         vicen_f(:,:,k)=fld2d;
         call get_mod_fld_nc(trim(icerestart), fld2d,'qice003', k, tlevel, nx, ny)
         vsnon_f(:,:,k)=fld2d;
      end do
      call replace_var_ncfile(aicen_f,nx,ny,ncat,'qice001',new_icerestart,3)
      print  '(a8)','qice002:' 
      call replace_var_ncfile(vicen_f,nx,ny,ncat,'qice002',new_icerestart,3)
      print  '(a8)','qice003:' 
      call replace_var_ncfile(vsnon_f,nx,ny,ncat,'qice003',new_icerestart,3)
      return
   endif

   tlevel=1
   do k=1,ncat
      call get_mod_fld_nc(trim(icerestart), fld2d,'aicen', k, tlevel, nx, ny)
      aicen_f(:,:,k)=fld2d;
      call get_mod_fld_nc(trim(icerestart), fld2d,'vicen', k, tlevel, nx, ny)
      vicen_f(:,:,k)=fld2d;
      call get_mod_fld_nc(trim(icerestart), fld2d,'vsnon', k, tlevel, nx, ny)
      vsnon_f(:,:,k)=fld2d;
      call get_mod_fld_nc(trim(icerestart), fld2d,'qsno001', k, tlevel, nx, ny)
      qsnon_f(:,:,k)=fld2d;
   end do

   aicen=aicen_f; vicen=vicen_f; vsnon=vsnon_f; qicen_f=qsnon_f;

   ! Step 0:  QC for the analyzed active 3-D ice variables: aicen, vicen,vsnon
   ficem=0
   ! adjusting aicen:
   do j=1,ny
     do i=1,nx
       if (Imask(i,j)==1.and.fice(i,j)>Athresh) then  ! filter the land point
          do k=1, ncat
             aicen(i,j,k)=max(0.,aicen(i,j,k))
             aicen(i,j,k)=min(1.,aicen(i,j,k))
             aicen_f(i,j,k)=max(0.,aicen_f(i,j,k))
             aicen_f(i,j,k)=min(1.,aicen_f(i,j,k))
          end do
          ficem(i,j)=sum(aicen(i,j,:))
          fice_f(i,j)=sum(aicen_f(i,j,:))
          ! require the minimal threshold for the distribution 
          do k=1, ncat
             ! require the minimal threshold for the distribution 
             if (aicen_f(i,j,k)<thresh.or.fice_f(i,j)<Athresh) aicen_f(i,j,k)=0.
          end do
          fice_f(i,j)=sum(aicen_f(i,j,:))
  
          if (fice_f(i,j)>Athresh) then
             ! --- rescaled by analyzed total ice concentration
             do k=1, ncat
                aicen(i,j,k)=aicen_f(i,j,k)*fice(i,j)/fice_f(i,j)
             end do
             ficem(i,j)=sum(aicen(i,j,:))
             if (ficem(i,j)>1) then
                do k=1, ncat
                   aicen(i,j,k)=aicen(i,j,k)/ficem(i,j)
                end do
             endif
          else  ! Generate the new ice
             if (Vadjust==1) then
                aicen(i,j,:)=0.; vicen(i,j,:)=0.; vsnon(i,j,:)=0.
                ! ! some idea: dealing with the grid with the new ice after DA
                ! using the ice distribution at the nearest neighbor grid 
                ! or using the new ice equation to introduce more new ice in
                ! future
             elseif (Vadjust==2) then
                ! ! some idea: dealing with the grid with the new ice after DA
                ! using the ice distribution at the nearest neighbor grid 
                ! No real test for that:
                call ij_neighborgrid(i,j,Imask,nx,ny,fice_f,Athresh,ibor,jbor)
                if (ibor>0.and.jbor>0) then
                   aicen(i,j,:)=aicen_f(ibor,jbor,:)*fice(i,j)/fice_f(ibor,jbor)
                   vicen(i,j,:)=vicen(ibor,jbor,:)*fice(i,j)/fice_f(ibor,jbor)
                   vsnon(i,j,:)=vsnon(ibor,jbor,:)*fice(i,j)/fice_f(ibor,jbor)
                   qicen_f(i,j,:)=qsnon_f(ibor,jbor,:)
                else
                   aicen(i,j,:)=0.; vicen(i,j,:)=0.; vsnon(i,j,:)=0.
                endif
       
             endif
          endif
          ficem(i,j)=sum(aicen(i,j,:))
       else
          aicen(i,j,:)=0.; vicen(i,j,:)=0.; vsnon(i,j,:)=0.
          aicen_f(i,j,:)=0.; vicen_f(i,j,:)=0.;  ficem(i,j)=0.; fice_f(i,j)=0.
       endif
     end do
   end do

   ! dealing with vicen / vsnon:
   do j=1,ny
      do i=1,nx
         if (Imask(i,j)==1.and.ficem(i,j)>0) then  ! filter the land point
            if (fice_f(i,j)>Athresh) then
               do k=1,ncat
                  if (aicen_f(i,j,k)>0.and.aicen(i,j,k)>0) then
                     vicen(i,j,k)=vicen(i,j,k)*aicen(i,j,k)/aicen_f(i,j,k)
           !         vsnon(i,j,k)=vsnon(i,j,k)*aicen(i,j,k)/aicen_f(i,j,k)
                  endif
               enddo

               ! When ice pack area assimilating HICE
               if (ficem(i,j)>.75.and.Ahice>0) then
                   Vtemp(i,j)=ficem(i,j)*(hice(i,j)*Ahice+sum(vicen(i,j,:))*(1-Ahice))/sum(vicen(i,j,:))
                   do k=1,ncat
                      if (aicen(i,j,k)>0) then
                         vicen(i,j,k)=vicen(i,j,k)*Vtemp(i,j)
                      else
                         vicen(i,j,k)=0.
                      endif
                   enddo
               endif 

               ! last check:
               do k=1,ncat
                  if (aicen(i,j,k)>0) then
                     vicen(i,j,k)=max(vicen(i,j,k),aicen(i,j,k)*bkcat(k))
                     vicen(i,j,k)=min(vicen(i,j,k),aicen(i,j,k)*bkcat(k+1))
                  else
                     vicen(i,j,k)=0; vsnon(i,j,k)=0; aicen(i,j,k)=0
                  endif
               enddo
            else
               if (Vadjust==1) then
                  vicen(i,j,:)=0.; vsnon(i,j,:)=0.
               elseif (Vadjust==2) then   ! considering the new ice generation by DA
                  do k=1,ncat
                     if (aicen(i,j,k)>0) then
                        vicen(i,j,k)=max(vicen(i,j,k),aicen(i,j,k)*bkcat(k))
                        vicen(i,j,k)=min(vicen(i,j,k),aicen(i,j,k)*bkcat(k+1))
                     endif
                     ! new fresh snow
            !         vsnon(i,j,k)=min(0.2*vicen(i,j,k),hs_min*aicen(i,j,k))
                  enddo
               endif
            endif
         else
            vicen(i,j,:)=0; vsnon(i,j,:)=0; aicen(i,j,:)=0; ficem(i,j)=0
         endif 
      end do
   end do

   !========================================================================
   ! finial control and replacing processing: 
   ! update the variables in iced file
   call replace_var_ncfile(aicen,nx,ny,ncat,'aicen',new_icerestart,3)
   call replace_var_ncfile(vicen,nx,ny,ncat,'vicen',new_icerestart,3)
   call replace_var_ncfile(vsnon,nx,ny,ncat,'vsnon',new_icerestart,3)

   ! for other 2D variables delimited by ficem
   do ll=1,Nv22
      cfld=trim(Var2d2(ll)%var); tlevel=1
      call get_mod_fld_nc(trim(icerestart), fld2d,cfld, 0, tlevel, nx, ny)
      where(ficem<=0.)
         fld2d=0.0
      end where           
      call replace_var_ncfile(fld2d,nx,ny,1,cfld,new_icerestart,2)
      print '(a8,i3)', trim(cfld),tlevel
   end do

   ! for other 3D(ncat) variables (different category) 
   do ll=1,Nv32
   !do ll=1,17
      cfld=trim(Var3d2(ll)%var); tlevel=1; Tupdate=.true.
      if (ll>17.and.ll<25) then
         Tupdate=.false.
         continue 
      endif
      do k=1,ncat
         call get_mod_fld_nc(trim(icerestart), fld2d,cfld, k, tlevel, nx, ny)
         if (ll>10) then
            call get_mod_fld_nc(trim(new_icerestart), fld2d_Tf,'Tsfcn', k, tlevel, nx, ny)
         endif
         fld2d1=aicen(:,:,k);  fld2d_f=aicen_f(:,:,k)
         if (ll/=10) then
            where(fld2d1<=0.or.ficem<=0)
               fld2d=0.0
            end where
         else
            where(fld2d1<=0.or.ficem<=0)
               fld2d=-1.8
            end where
         endif
         do j=1,ny
            do i=1,nx
               if (aicen(i,j,k)>0) then
                  select case (ll)
                     case(11:17)
                        fld2d(i,j)=ModzSin(ll-10)
!                        case(18:24)
!                           if (aicen_f(i,j,k)<=thresh) then
!                              Tmz=ModTmlt(ll-17)   ! Tmlt
!                              ! replace the previous usage from Madlen: Tmz-fld2d_Tf(i,j)
!                              ! at 6th May 2022
!                              tmp1=Tf(i,j)+1.8
!                              !tmp1=Tmz-fld2d_Tf(i,j)
!                              Ti=-1.8+tmp1*(real(ll-17)-0.5)/real(Nilayer)
!
!                              fld2d(i,j)=-rhoi*(cp_ice*(Tmz-Ti)+Lfresh*(1.0-Tmz/Ti) &
!                                      -cp_ocn*Tmz)/real(Nilayer)
!                           endif
              !       case(18:24)
              !          if (ficem(i,j)>.75.and.Ahice>0) then
              !             fld2d(i,j)=fld2d(i,j)/Vtemp(i,j)
              !          endif
              !       case(25)    !snow
               !         if (fice_f(i,j)<=Athresh.and.Vadjust==2) then
               !            fld2d(i,j)=qicen_f(i,j,k)
               !         endif
              !             if (aicen_f(i,j,k)<thresh) then
              !                Ti=min(0.,fld2d_Tf(i,j))
              !                fld2d(i,j)=-rhos*(Lfresh-cp_ice*Ti)
!             !                 fld2d(i,j)=-rhos*Lfresh
              !             endif
                  end select
               else
                  if (ll/=10) then
                     fld2d(i,j)=0
                  else
                     fld2d(i,j)=-1.8
                  endif

               endif
            end do
         end do

         if (ll>17.and.ll<24) then
            print *, 'final: '
            print '(a8,i3,a10,e15.5,a8,i3)','qice',ll,'(510,580)',fld2d(ist,jst),' at k=',k
         endif
         fld3d(:,:,k)=fld2d
      end do
      if (Tupdate==.true.) then
         print  '(a8,i3)', trim(cfld),tlevel
         call replace_var_ncfile(fld3d,nx,ny,ncat,cfld,new_icerestart,3)
      endif
   end do

   deallocate(aicen_f,vicen_f,vsnon_f,fice_f,hice_f)
   deallocate(fld2d_f,fld2d_Tf,Tf)
   deallocate(ModzSin,ModTmlt,bkcat,bkaice)
   deallocate(qicen_f,qsnon_f)
   deallocate(Vtemp)
end subroutine fix_cice



!----------------------------------------------------------------------
subroutine ij_neighborgrid(ipiv,jpiv,Vmask,idm,jdm,fldr,Vthresh,ii,jj)
   implicit none
   integer,     parameter :: Grid_radius=5
   integer,    intent(in) :: ipiv,jpiv
   integer,    intent(in) :: idm,jdm
   integer,    intent(in) :: Vmask(idm,jdm)  ! Valid water mask if equals 1
   real,       intent(in) :: fldr(idm,jdm)
   real,       intent(in) :: Vthresh 
   integer,   intent(out) :: ii,jj
   !------------------------------
   integer :: i,j,k
   integer :: i1,i2,j1,j2
   logical :: Lgo
   ii=0;  jj=0
   k=1;  Lgo=.true.
   do while (k<=Grid_radius)
      j1=max(1,jpiv-k); j2=min(jdm,jpiv+k)
      ! Western
      i=max(1,ipiv-k); j=j1
      do while (j<=j2.and.Lgo)
         if (Vmask(i,j)==1.and.fldr(i,j)>Vthresh) then
            ii=i; jj=j
            j=j2;  Lgo=.false.
         endif
         j=j+1
      end do
      ! Eastern
      i=min(idm,ipiv+k); j=j1
      do while (j<=j2.and.Lgo)
         if (Vmask(i,j)==1.and.fldr(i,j)>Vthresh) then
            ii=i; jj=j
            j=j2;  Lgo=.false.
         endif
         j=j+1
      end do

      i1=max(1,ipiv-k);  i2=min(idm,ipiv+k); 
      ! Southern
      j=max(1,jpiv-k);  i=i1
      do while (i<=i2.and.Lgo)
         if (Vmask(i,j)==1.and.fldr(i,j)>Vthresh) then
            ii=i; jj=j
            i=i2;  Lgo=.false.
         endif
         i=i+1
      end do
      ! Northern
      j=min(jdm,jpiv+k);  i=i1
      do while (i<=i2.and.Lgo)
         if (Vmask(i,j)==1.and.fldr(i,j)>Vthresh) then
            ii=i; jj=j
            i=i2;  Lgo=.false.
         endif
         i=i+1
      end do
      if (ii>0.and.jj>0) then
         k=Grid_radius
      endif
      k=k+1
   end do
end subroutine



! itype equals 2 for 2-D variable; 3 for 3-D variable
subroutine replace_var_ncfile(fldr,nx,ny,nz,cfld,ncfile,itype,Lmask)
   implicit none
   integer,                    parameter :: itst=40, jtst=40
   integer,                   intent(in) :: nx,ny,nz
   integer,                   intent(in) :: itype
   real, dimension(nx,ny,nz), intent(in) :: fldr
   character(len=*),          intent(in) :: cfld
   character(len=80),         intent(in) :: ncfile
   integer, optional,         intent(in) :: Lmask(nx,ny)
   
   integer   :: ncid, vtmp_id
   integer   :: ns2(2), nc2(2)
   integer   :: ns3(3), nc3(3)

   integer   :: i,j,k
   real, dimension(nx,ny)    :: fld2
   real, dimension(nx,ny,nz) :: fld0,fld3
   
   ns2=1; nc2(1)=nx; nc2(2)=ny; 
   ns3=1; 
   nc3(1)=nx; nc3(2)=ny;  nc3(3)=nz; 

   call nfw_open(trim(ncfile),nf_write, ncid)
   call nfw_inq_varid(trim(ncfile),ncid,trim(cfld),vtmp_id)

   if(present(Lmask)) then
      fld0=fldr
      if (itype==2) then
         call nfw_get_vara_double(trim(ncfile), ncid, vtmp_id, ns2, nc2, fld2)
         !print *, 'get mask ', trim(cfld),fld2(itst,jtst),Lmask(itst,jtst)
      elseif (itype==3) then   ! 3D variable
         call nfw_get_vara_double(trim(ncfile), ncid, vtmp_id, ns3, nc3, fld3)
      end if
      do j=1,ny
         do i=1,nx
            if (Lmask(i,j)/=1) then  ! land point 
               if (itype==2) then
                  fld0(i,j,:)=fld2(i,j)
               elseif (itype==3) then
                  fld0(i,j,:)=fld3(i,j,:)
               end if
            end if 
         end do
      end do
   endif

   if (itype==2) then   ! 2D variable
     if (present(Lmask)) then
!        print *, 'replace mask ', trim(cfld),fld0(itst,jtst,1)
        call nfw_put_vara_double(trim(ncfile), ncid, vtmp_id, ns2, nc2, fld0)
     else
        call nfw_put_vara_double(trim(ncfile), ncid, vtmp_id, ns2, nc2, fldr)
     endif
   elseif (itype==3) then   ! 3D variable
     if (present(Lmask)) then
        call nfw_put_vara_double(trim(ncfile), ncid, vtmp_id, ns3, nc3, fld0)
     else
        call nfw_put_vara_double(trim(ncfile), ncid, vtmp_id, ns3, nc3, fldr)
     endif
   else
     print *, 'Wrong input replace_var_ncfile()'
     stop
   end if
   call nfw_close(trim(ncfile), ncid)

end subroutine replace_var_ncfile


end module m_put_mod_fld_nc

