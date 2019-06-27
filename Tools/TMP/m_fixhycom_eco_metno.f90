










module m_fixhycom_eco_metno
!Ehouarn: March 2011
!fixanalysis: remapping of tracers after physical analysis!
!use of remapping subroutines embedded in hycom to interpolate!
!biogeochemical tracer on the analysis grid (dp)!
!Remapping is realized after correction of negative anlaysis dp!

contains
      subroutine hybgen_weno_coefs(s,dp,lc,ci,kk,ks,thin)
      implicit none

      integer kk,ks
      logical lc(kk)
      real    s(kk,ks),dp(kk),ci(kk,ks,2),thin
!
!-----------------------------------------------------------------------
!  1) coefficents for remaping from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the 
!             interfacial values 
!
!     REFERENCE?
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       dp    - initial layer thicknesses (>=thin)
!       lc    - use PCM for selected layers
!       kk    - number of layers
!       ks    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       ci    - coefficents for hybgen_weno_remap
!                ci.1 is value at interface above
!                ci.2 is value at interface below
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
!-----------------------------------------------------------------------
!
      real, parameter :: dsmll=1.0e-8
!
      integer j,i
      real    q,q01,q02,q001,q002
      real    qdpjm(kk),qdpjmjp(kk),dpjm2jp(kk)
      real    zw(kk+1,3)

      !compute grid metrics
      do j=2,kk-1
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
        qdpjmjp(j) = 1.0/(dp(j-1) +     dp(j) + dp(j+1))
        dpjm2jp(j) =      dp(j-1) + 2.0*dp(j) + dp(j+1)
      enddo !j
      j=kk
        qdpjm(  j) = 1.0/(dp(j-1) +     dp(j))
!
      do i= 1,ks
        do j=2,kk
          zw(j,3) = qdpjm(j)*(s(j,i)-s(j-1,i))
        enddo !j
          j = 1  !PCM first layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
        do j=2,kk-1
          if     (lc(j) .or. dp(j).le.thin) then  !use PCM
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0
          else
            q001 = dp(j)*zw(j+1,3)
            q002 = dp(j)*zw(j,  3)
            if (q001*q002 < 0.0) then
              q001 = 0.0
              q002 = 0.0
            endif
            q01 = dpjm2jp(j)*zw(j+1,3)
            q02 = dpjm2jp(j)*zw(j,  3)
            if     (abs(q001) > abs(q02)) then
              q001 = q02
            endif
            if     (abs(q002) > abs(q01)) then
              q002 = q01
            endif
            q    = (q001-q002)*qdpjmjp(j)
            q001 = q001-q*dp(j+1)
            q002 = q002+q*dp(j-1)

            ci(j,i,2) = s(j,i)+q001
            ci(j,i,1) = s(j,i)-q002
            zw(  j,1) = (2.0*q001-q002)**2
            zw(  j,2) = (2.0*q002-q001)**2
          endif  !PCM:WEND
        enddo !j
          j = kk  !PCM last layer
            ci(j,i,1) = s(j,i)
            ci(j,i,2) = s(j,i)
            zw(j,  1) = 0.0
            zw(j,  2) = 0.0

        do j=2,kk
          q002 = max(zw(j-1,2),dsmll)
          q001 = max(zw(j,  1),dsmll)
          zw(j,3) = (q001*ci(j-1,i,2)+q002*ci(j,i,1))/(q001+q002)
        enddo !j
          zw(   1,3) = 2.0*s( 1,i)-zw( 2,3)  !not used?
          zw(kk+1,3) = 2.0*s(kk,i)-zw(kk,3)  !not used?

        do j=2,kk-1
          if     (.not.(lc(j) .or. dp(j).le.thin)) then  !don't use PCM
            q01  = zw(j+1,3)-s(j,i)
            q02  = s(j,i)-zw(j,3)
            q001 = 2.0*q01
            q002 = 2.0*q02
            if     (q01*q02 < 0.0) then
              q01 = 0.0
              q02 = 0.0
            elseif (abs(q01) > abs(q002)) then
              q01 = q002
            elseif (abs(q02) > abs(q001)) then
              q02 = q001
            endif
            ci(j,i,1) = s(j,i)-q02
            ci(j,i,2) = s(j,i)+q01
          endif  !PCM:WEND
        enddo !j
      enddo !i
      return
      end subroutine hybgen_weno_coefs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hybgen_weno_remap(si,pi,dpi,ci,&
                                  so,po,ki,ko,ks,thin)
      implicit none
!
      integer ki,ko,ks
      real    si(ki,ks),pi(ki+1),dpi(ki),ci(ki,ks,2),&
             so(ko,ks),po(ko+1),thin
!
!-----------------------------------------------------------------------
!  1) remap from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the 
!             interfacial values 
!             the output is the average of the interpolation
!             profile across each output cell.
!
!     REFERENCE?
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       pi    - initial layer interface depths (non-negative)
!                  pi(   1) is the surface
!                  pi(ki+1) is the bathymetry
!                  pi(k+1) >= pi(k)
!       dpi   - initial layer thicknesses (dpi(k)=pi(k+1)-pi(k))
!       ci    - coefficents from hybgen_weno_coefs
!                ci.1 is value at interface above
!                ci.2 is value at interface below
!       ki    - number of  input layers
!       ko    - number of output layers
!       ks    - number of fields
!       po    - target interface depths (non-negative)
!                  po(   1) is the surface
!                  po(ko+1) is the bathymetry (== pi(ki+1))
!                  po(k+1) >= po(k)
!       thin  - layer thickness (>0) that can be ignored
!
!  3) output arguments:
!       so    - scalar fields in po-layer space
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
      integer i,k,l,lb,lt
      real    dpb,dpt,qb0,qb1,qb2,qt0,qt1,qt2,xb,xt,zb,zt,zx,o
      real*8  sz
!
      zx=pi(ki+1) !maximum depth
      zb=max(po(1),pi(1))
      lb=1
      do while (pi(lb+1).lt.zb .and. lb.lt.ki)
        lb=lb+1
      enddo
      do k= 1,ko  !output layers
        zt = zb
        zb = min(po(k+1),zx)
!       write(lp,*) 'k,zt,zb = ',k,zt,zb
        lt=lb !top will always correspond to bottom of previous
              !find input layer containing bottom output interface
        do while (pi(lb+1).lt.zb .and. lb.lt.ki)
          lb=lb+1
        enddo
        if     (zb-zt.le.thin .or. zt.ge.zx) then
          if     (k.ne.1) then
!
! ---       thin or bottomed layer, values taken from layer above
!
            do i= 1,ks
              so(k,i) = so(k-1,i)
            enddo !i
          else !thin surface layer
            do i= 1,ks
              so(k,i) = si(k,i)
            enddo !i
          endif
        else

!         form layer averages.
!
!         if     (pi(lb).gt.zt) then
!           write(lp,*) 'bad lb = ',lb
!           stop
!         endif
          xt=(zt-pi(lt))/max(dpi(lt),thin)
          xb=(zb-pi(lb))/max(dpi(lb),thin)
          if     (lt.ne.lb) then  !multiple layers
            dpt = pi(lt+1)-zt
            dpb = zb-pi(lb)
            qt1 = xt*(xt-1.0)
            qt2 = qt1+xt
            qt0 = 1.0-qt1-qt2
            qb1 = (xb-1.0)**2
            qb2 = qb1-1.0+xb
            qb0 = 1.0-qb1-qb2
            do i= 1,ks
              o = si((lt+lb)/2,i)  !offset to reduce round-off
              sz = dpt*(qt0*(si(lt,i)  -o) + &
                       qt1*(ci(lt,i,1)-o) + &
                       qt2*(ci(lt,i,2)-o)  )
              do l=lt+1,lb-1
                sz = sz+dpi(l)*(si(l,i) - o)
              enddo !l
              sz  = sz + dpb*(qb0*(si(lb,i)  -o) + &
                             qb1*(ci(lb,i,1)-o) + &
                             qb2*(ci(lb,i,2)-o)  )
              so(k,i) = o + sz/(zb-zt)  !zb-zt>=thin
            enddo !i
          else !single layer
            qt1 = xb**2 + xt**2 + xb*xt + 1.0 - 2.0*(xb+xt)
            qt2 = qt1 - 1.0 + (xb+xt)
            qt0 = 1.0 - qt1 - qt2
            do i= 1,ks
              sz=qt0*(si(lt,i)  -o) + &
                qt1*(ci(lt,i,1)-o) + &
                qt2*(ci(lt,i,2)-o) 
              so(k,i) = o + sz
            enddo !i
          endif !layers
        endif !thin:std layer
      enddo !k
      return
      end subroutine hybgen_weno_remap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      integer function tracr_get_incr(char2) 
      !function returns the number of the tracer.
      !char => integer
      implicit none        
      character(len=2) :: char2
      
      tracr_get_incr=-1
      select case (char2)
         case ('01')
           tracr_get_incr=1
	 case ('02')
           tracr_get_incr=2
	 case ('03')
           tracr_get_incr=3
	 case ('04')
           tracr_get_incr=4
	 case ('05')
           tracr_get_incr=5
	 case ('06')
           tracr_get_incr=6
	 case ('07')
           tracr_get_incr=7
	 case ('08')
           tracr_get_incr=8
	 case ('09')
           tracr_get_incr=9 
	 case ('10')
           tracr_get_incr=10
	 case ('11')
           tracr_get_incr=11
	 case default
           print *,'tracer unknown',char2
      end select
      return
      
      end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      integer function compute_kisop(temp,sal,nz) 
      !function defines which layers are isopycnal 
      implicit none             
      integer::nz
      real,dimension(1:nz)::temp,sal
      real::eps,tmp
      integer::k
      
      eps=0.1
      
      do k=2,nz-1
        if (sig(temp(k),sal(k))-sig_ref(k).lt.0.)then
	  tmp=(sig(temp(k),sal(k))-sig_ref(k))/(sig_ref(k)-sig_ref(k-1))
	  if (tmp.gt.eps)then
	    compute_kisop=k
	    return
	  else
	    tmp=(sig(temp(k),sal(k))-sig_ref(k))/(sig_ref(k+1)-sig_ref(k))
	    if(tmp.gt.eps)then
	     compute_kisop=k
	     return
	    endif 
	  endif
	endif
      enddo
      compute_kisop=nz
      return
      
      end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1    
      
      real function sig(t,s)
      !function returns the value of sigma_0
      !according to T and S
      implicit none
      real::t,s
      real :: c1,c2,c3,c4,c5,c6,c7
        
      c1=-1.36471E-01  
      c2= 4.68181E-02  
      c3= 8.07004E-01  
      c4=-7.45353E-03  
      c5=-2.94418E-03  
      c6= 3.43570E-05  
      c7= 3.48658E-05  
      
      sig=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))         
      return	 
      
      end function
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      real function sig_ref(k)
      !function return the value of the target density for a given layer
      implicit none
      integer::k
      
      select case (k) 
        case (1)
	  sig_ref=0.1
	case (2)
	  sig_ref=0.2 
	case (3)
	  sig_ref=0.3 	  
	case (4)
	  sig_ref=0.4 	   
 	case (5)
	  sig_ref=0.5 
        case (6)
	  sig_ref=24.05
	case (7)
	  sig_ref=24.96 
	case (8)
	  sig_ref=25.68 	  
	case (9)
	  sig_ref=26.05 	   
 	case (10)
	  sig_ref=26.30 	  
        case (11)
	  sig_ref=26.60
	case (12)
	  sig_ref=26.83 
	case (13)
	  sig_ref=27.03 	  
	case (14)
	  sig_ref=27.20 	   
 	case (15)
	  sig_ref=27.33 
        case (16)
	  sig_ref=27.46
	case (17)
	  sig_ref=27.55 
	case (18)
	  sig_ref=27.66 	  
	case (19)
	  sig_ref=27.74 	   
 	case (20)
	  sig_ref=27.82 	  
	case (21)
	  sig_ref=27.90
	case (22)
	  sig_ref=27.97 
	case (23)
	  sig_ref=28.01 	  
	case (24)
	  sig_ref=28.04 	   
 	case (25)
	  sig_ref=28.07 
        case (26)
	  sig_ref=28.09
	case (27)
	  sig_ref=28.11 
	case (28)
	  sig_ref=28.13 
       
      end select    
	       
      return	 
      
      end function
 
end module
