module modal_aero_wateruptake

!  Portable science module for modal aerosol water uptake.
!  Contains Kohler theory wet radius calculation and polynomial solvers.
!
!  RCE 07.04.13:  Adapted from MIRAGE2 code

use ccpp_kinds, only: kind_phys

implicit none
private
save

public :: &
   modal_aero_wateruptake_init,  &
   modal_aero_wateruptake_sub,   &
   modal_aero_wateruptake_diag

real(kind_phys), parameter :: third = 1._kind_phys/3._kind_phys
real(kind_phys) :: pi43

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_wateruptake_init(pi, errmsg, errflg)
   real(kind_phys),         intent(in)  :: pi
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   errmsg = ''
   errflg = 0

   pi43 = pi*4.0_kind_phys/3.0_kind_phys

end subroutine modal_aero_wateruptake_init

subroutine modal_aero_wateruptake_sub( &
   aero_props, aero_state, &
   ncol, pver, top_lev, &
   do_strat_sulfate, &
   t, pmid, h2ommr, cldn, &
   dryrad, hygro, dryvol, so4dryvol, &
   dgncur_awet, troplev, &
   wetrad, wetvol, wtrvol, &
   sulfeq, wtpct, sulden, &
   specdens_1, alnsg_out, maer, &
   errmsg, errflg)

!-----------------------------------------------------------------------
!
! Purpose: Compute aerosol wet radius
!
! Method:  Kohler theory
!
! Author:  S. Ghan
!
!-----------------------------------------------------------------------

   use aerosol_properties_mod, only: aerosol_properties
   use aerosol_state_mod,      only: aerosol_state
   use wv_saturation,          only: qsat_water

   ! Arguments
   class(aerosol_properties), intent(in) :: aero_props
   class(aerosol_state),      intent(in) :: aero_state
   integer,          intent(in)  :: ncol                    ! number of columns
   integer,          intent(in)  :: pver                    ! number of vertical levels
   integer,          intent(in)  :: top_lev                 ! top level for aerosol calculations
   logical,          intent(in)  :: do_strat_sulfate        ! use stratospheric sulfate treatment
   real(kind_phys),         intent(in)  :: t(:,:)                  ! temperature (K)
   real(kind_phys),         intent(in)  :: pmid(:,:)               ! layer pressure (Pa)
   real(kind_phys),         intent(in)  :: h2ommr(:,:)             ! specific humidity (kg/kg)
   real(kind_phys),         intent(in)  :: cldn(:,:)               ! cloud fraction (0-1)
   real(kind_phys),         intent(in)  :: dryrad(:,:,:)           ! dry volume mean radius of aerosol (m)
   real(kind_phys),         intent(in)  :: hygro(:,:,:)            ! volume-weighted mean hygroscopicity (--)
   real(kind_phys),         intent(in)  :: dryvol(:,:,:)           ! dry volume of single aerosol (m3)
   real(kind_phys),         intent(in)  :: so4dryvol(:,:,:)        ! dry volume of sulfate in single aerosol (m3)
   real(kind_phys),         intent(in)  :: dgncur_awet(:,:,:)      ! wet number mode diameter, prev timestep (m)
   integer,          intent(in)  :: troplev(:)              ! tropopause level index

   real(kind_phys),         intent(out) :: wetrad(:,:,:)           ! wet radius of aerosol (m)
   real(kind_phys),         intent(out) :: wetvol(:,:,:)           ! single-particle-mean wet volume (m3)
   real(kind_phys),         intent(out) :: wtrvol(:,:,:)           ! single-particle-mean water volume (m3)
   real(kind_phys),         intent(out) :: sulfeq(:,:,:)           ! H2SO4 equilibrium mixing ratio (mol/mol)
   real(kind_phys),         intent(out) :: wtpct(:,:,:)            ! sulfate composition, weight % H2SO4
   real(kind_phys),         intent(out) :: sulden(:,:,:)           ! sulfate aerosol mass density (g/cm3)
   real(kind_phys),         intent(out) :: specdens_1(:)           ! first-species density per mode (kg/m3)
   real(kind_phys),         intent(out) :: alnsg_out(:)            ! log(sigma_g) per mode
   real(kind_phys),         intent(out) :: maer(:,:,:)             ! accumulated mode mass (kg/kg)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! local variables

   integer  :: i, k, l, m
   integer  :: nmodes, nspec
   real(kind_phys) :: so4specdens, specdens, sigmag
   character(len=32) :: spectype
   real(kind_phys) :: dmean, qh2so4_equilib, wtpct_mode, sulden_mode
   real(kind_phys) :: hystfac

   real(kind_phys), allocatable :: rhcrystal(:), rhdeliques(:)
   real(kind_phys) :: rh(ncol, pver)
   real(kind_phys) :: es(ncol), qs(ncol)

   real(kind_phys), pointer :: raer(:,:)
   !-----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   nmodes = aero_props%nbins()

   ! Initialize defaults
   wtpct(:,:,:)  = 75._kind_phys
   sulden(:,:,:) = 1.923_kind_phys
   sulfeq(:,:,:) = 0._kind_phys
   maer(:,:,:)   = 0._kind_phys

   allocate(rhcrystal(nmodes), rhdeliques(nmodes))

   so4specdens = 0._kind_phys

   ! Query mode properties and accumulate mode mass
   do m = 1, nmodes

      sigmag = exp(aero_props%alogsig(m))
      rhcrystal(m) = aero_props%rhcrystal(m)
      rhdeliques(m) = aero_props%rhdeliques(m)
      alnsg_out(m) = log(sigmag)

      ! get mode info
      nspec = aero_props%nspecies(m)

      do l = 1, nspec

         ! accumulate the aerosol masses of each mode
         call aero_state%get_ambient_mmr(species_ndx=l, bin_ndx=m, mmr=raer)
         maer(:ncol,:,m) = maer(:ncol,:,m) + raer(:ncol,:)

         ! get species interstitial mixing ratio ('a')
         call aero_props%get(m, l, density=specdens, &
                                     spectype=spectype)

         if (do_strat_sulfate .and. (trim(spectype).eq.'sulfate')) then
            so4specdens=specdens
         end if

         if (l == 1) then
            ! save off these values to be used as defaults
            specdens_1(m) = specdens
         end if

      end do

      ! Compute stratospheric sulfate equilibrium
      if (do_strat_sulfate) then
         do k = top_lev, pver
            do i = 1, ncol
               dmean = dgncur_awet(i,k,m)*exp(1.5_kind_phys*alnsg_out(m)**2)
               call calc_h2so4_equilib_mixrat( t(i,k), pmid(i,k), h2ommr(i,k), dmean, &
                                               qh2so4_equilib, wtpct_mode, sulden_mode, &
                                               errmsg, errflg )
               if (errflg /= 0) then
                  deallocate(rhcrystal, rhdeliques)
                  return
               end if
               sulfeq(i,k,m)  = qh2so4_equilib
               wtpct(i,k,m)   = wtpct_mode
               sulden(i,k,m)  = sulden_mode
            end do    ! i = 1, ncol
         end do    ! k = top_lev, pver
      end if

   end do    ! m = 1, nmodes

   ! Compute relative humidity
   do k = top_lev, pver
      call qsat_water(t(1:ncol,k), pmid(1:ncol,k), es(1:ncol), qs(1:ncol), ncol)
      do i = 1, ncol
         if (qs(i) > h2ommr(i,k)) then
            rh(i,k) = h2ommr(i,k)/qs(i)
         else
            rh(i,k) = 0.98_kind_phys
         endif
         rh(i,k) = max(rh(i,k), 0.0_kind_phys)
         rh(i,k) = min(rh(i,k), 0.98_kind_phys)
         if (cldn(i,k) .lt. 1.0_kind_phys) then
            rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_kind_phys - cldn(i,k))  ! clear portion
         end if
         rh(i,k) = max(rh(i,k), 0.0_kind_phys)
      end do
   end do

   ! Compute wet radius using Kohler theory
   do m = 1, nmodes

      hystfac = 1.0_kind_phys / max(1.0e-5_kind_phys, (rhdeliques(m) - rhcrystal(m)))

      do k = top_lev, pver
         do i = 1, ncol

            if ( do_strat_sulfate .and. (k<troplev(i))) then
               wetvol(i,k,m) = dryvol(i,k,m)-so4dryvol(i,k,m)
               wetvol(i,k,m) = wetvol(i,k,m)+so4dryvol(i,k,m)*so4specdens/sulden(i,k,m)/wtpct(i,k,m)/10._kind_phys
               wetvol(i,k,m) = max(wetvol(i,k,m), dryvol(i,k,m))
               wetrad(i,k,m) = (wetvol(i,k,m)/pi43)**third
               wetrad(i,k,m) = max(wetrad(i,k,m), dryrad(i,k,m))
               wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m)
               wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_kind_phys)
            else
              ! compute wet radius for each mode
              call modal_aero_kohler(dryrad(i:i,k,m), hygro(i:i,k,m), rh(i:i,k), wetrad(i:i,k,m), 1)

              wetrad(i,k,m) = max(wetrad(i,k,m), dryrad(i,k,m))
              wetvol(i,k,m) = pi43*wetrad(i,k,m)**3
              wetvol(i,k,m) = max(wetvol(i,k,m), dryvol(i,k,m))
              wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m)
              wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_kind_phys)

              ! apply simple treatment of deliquesence/crystallization hysteresis
              if (rh(i,k) < rhcrystal(m)) then
                 wetrad(i,k,m) = dryrad(i,k,m)
                 wetvol(i,k,m) = dryvol(i,k,m)
                 wtrvol(i,k,m) = 0.0_kind_phys
              else if (rh(i,k) < rhdeliques(m)) then
                 wtrvol(i,k,m) = wtrvol(i,k,m)*hystfac*(rh(i,k) - rhcrystal(m))
                 wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_kind_phys)
                 wetvol(i,k,m) = dryvol(i,k,m) + wtrvol(i,k,m)
                 wetrad(i,k,m) = (wetvol(i,k,m)/pi43)**third
              end if
            end if

         end do  ! columns
      end do     ! levels

   end do ! modes

   deallocate(rhcrystal, rhdeliques)

end subroutine modal_aero_wateruptake_sub

!-----------------------------------------------------------------------
subroutine modal_aero_kohler(rdry_in, hygro, s, rwet_out, im)
   ! calculates equlibrium radius r of haze droplets as function of
   ! dry particle mass and relative humidity s using kohler solution
   ! given in pruppacher and klett (eqn 6-35)

   ! for multiple aerosol types, assumes an internal mixture of aerosols

   ! arguments
   integer :: im         ! number of grid points to be processed
   real(kind_phys) :: rdry_in(:)    ! aerosol dry radius (m)
   real(kind_phys) :: hygro(:)      ! aerosol volume-mean hygroscopicity (--)
   real(kind_phys) :: s(:)          ! relative humidity (1 = saturated)
   real(kind_phys) :: rwet_out(:)   ! aerosol wet radius (m)

   ! local variables
   integer, parameter :: imax=200
   integer :: i, n, nsol

   real(kind_phys) :: a, b
   real(kind_phys) :: p40(imax),p41(imax),p42(imax),p43(imax) ! coefficients of polynomial
   real(kind_phys) :: p30(imax),p31(imax),p32(imax) ! coefficients of polynomial
   real(kind_phys) :: p
   real(kind_phys) :: r3, r4
   real(kind_phys) :: r(im)         ! wet radius (microns)
   real(kind_phys) :: rdry(imax)    ! radius of dry particle (microns)
   real(kind_phys) :: ss            ! relative humidity (1 = saturated)
   real(kind_phys) :: slog(imax)    ! log relative humidity
   real(kind_phys) :: vol(imax)     ! total volume of particle (microns**3)
   real(kind_phys) :: xi, xr

   complex(kind_phys) :: cx4(4,imax),cx3(3,imax)

   real(kind_phys), parameter :: eps = 1.e-4_kind_phys
   real(kind_phys), parameter :: mw = 18._kind_phys
   real(kind_phys), parameter :: pi_local = 3.14159_kind_phys
   real(kind_phys), parameter :: rhow = 1._kind_phys
   real(kind_phys), parameter :: surften = 76._kind_phys
   real(kind_phys), parameter :: tair = 273._kind_phys
   real(kind_phys), parameter :: ugascon = 8.3e7_kind_phys


!     effect of organics on surface tension is neglected
   a=2.e4_kind_phys*mw*surften/(ugascon*tair*rhow)

   do i=1,im
        rdry(i) = rdry_in(i)*1.0e6_kind_phys   ! convert (m) to (microns)
        vol(i) = rdry(i)**3          ! vol is r**3, not volume
        b = vol(i)*hygro(i)

!          quartic
        ss=min(s(i),1._kind_phys-eps)
        ss=max(ss,1.e-10_kind_phys)
        slog(i)=log(ss)
        p43(i)=-a/slog(i)
        p42(i)=0._kind_phys
        p41(i)=b/slog(i)-vol(i)
        p40(i)=a*vol(i)/slog(i)
!          cubic for rh=1
        p32(i)=0._kind_phys
        p31(i)=-b/a
        p30(i)=-vol(i)
   end do


    do 100 i=1,im

!       if(vol(i).le.1.e-20)then
     if(vol(i).le.1.e-12_kind_phys)then
        r(i)=rdry(i)
        go to 100
     endif

     p=abs(p31(i))/(rdry(i)*rdry(i))
     if(p.lt.eps)then
!          approximate solution for small particles
        r(i)=rdry(i)*(1._kind_phys+p*third/(1._kind_phys-slog(i)*rdry(i)/a))
     else
        call makoh_quartic(cx4(1,i),p43(i),p42(i),p41(i),p40(i),1)
!          find smallest real(kind_phys) solution
        r(i)=1000._kind_phys*rdry(i)
        nsol=0
        do n=1,4
           xr=real(cx4(n,i))
           xi=aimag(cx4(n,i))
           if(abs(xi).gt.abs(xr)*eps) cycle
           if(xr.gt.r(i)) cycle
           if(xr.lt.rdry(i)*(1._kind_phys-eps)) cycle
           if(xr.ne.xr) cycle
           r(i)=xr
           nsol=n
        end do
        if(nsol.eq.0)then
           write(*,*)   &
            'ccm kohlerc - no real(kind_phys) solution found (quartic)'
           write(*,*)'roots =', (cx4(n,i),n=1,4)
           write(*,*)'p0-p3 =', p40(i), p41(i), p42(i), p43(i)
           write(*,*)'rh=',s(i)
           write(*,*)'setting radius to dry radius=',rdry(i)
           r(i)=rdry(i)
        endif
     endif

     if(s(i).gt.1._kind_phys-eps)then
!          save quartic solution at s=1-eps
        r4=r(i)
!          cubic for rh=1
        p=abs(p31(i))/(rdry(i)*rdry(i))
        if(p.lt.eps)then
           r(i)=rdry(i)*(1._kind_phys+p*third)
        else
           call makoh_cubic(cx3,p32,p31,p30,im)
!             find smallest real(kind_phys) solution
           r(i)=1000._kind_phys*rdry(i)
           nsol=0
           do n=1,3
              xr=real(cx3(n,i))
              xi=aimag(cx3(n,i))
              if(abs(xi).gt.abs(xr)*eps) cycle
              if(xr.gt.r(i)) cycle
              if(xr.lt.rdry(i)*(1._kind_phys-eps)) cycle
              if(xr.ne.xr) cycle
              r(i)=xr
              nsol=n
           end do
           if(nsol.eq.0)then
              write(*,*)   &
               'ccm kohlerc - no real(kind_phys) solution found (cubic)'
              write(*,*)'roots =', (cx3(n,i),n=1,3)
              write(*,*)'p0-p2 =', p30(i), p31(i), p32(i)
              write(*,*)'rh=',s(i)
              write(*,*)'setting radius to dry radius=',rdry(i)
              r(i)=rdry(i)
           endif
        endif
        r3=r(i)
!          now interpolate between quartic, cubic solutions
        r(i)=(r4*(1._kind_phys-s(i))+r3*(s(i)-1._kind_phys+eps))/eps
     endif

100 continue

   ! bound and convert from microns to m
   do i=1,im
      r(i) = min(r(i),30._kind_phys) ! upper bound based on 1 day lifetime
      rwet_out(i) = r(i)*1.e-6_kind_phys
   end do

end subroutine modal_aero_kohler


!-----------------------------------------------------------------------
subroutine makoh_cubic( cx, p2, p1, p0, im )
!
!     solves  x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(kind_phys) :: p0(imx), p1(imx), p2(imx)
      complex(kind_phys) :: cx(3,imx)

      integer :: i
      real(kind_phys) :: q(imx), r(imx), sqrt3
      complex(kind_phys) :: ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

      real(kind_phys), parameter :: eps = 1.e-20_kind_phys

      ci=cmplx(0._kind_phys,1._kind_phys,kind_phys)
      sqrt3=sqrt(3._kind_phys)
      cw=0.5_kind_phys*(-1+ci*sqrt3)
      cwsq=0.5_kind_phys*(-1-ci*sqrt3)

      do i=1,im
      if(p1(i).eq.0._kind_phys)then
!        completely insoluble particle
         cx(1,i)=(-p0(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
      else
         q(i)=p1(i)/3._kind_phys
         r(i)=p0(i)/2._kind_phys
         crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
         crad(i)=sqrt(crad(i))

         cy(i)=r(i)-crad(i)
         if (abs(cy(i)).gt.eps) cy(i)=cy(i)**third
         cq=q(i)
         cz(i)=-cq/cy(i)

         cx(1,i)=-cy(i)-cz(i)
         cx(2,i)=-cw*cy(i)-cwsq*cz(i)
         cx(3,i)=-cwsq*cy(i)-cw*cz(i)
      endif
      enddo

end subroutine makoh_cubic


!-----------------------------------------------------------------------
subroutine makoh_quartic( cx, p3, p2, p1, p0, im )

!     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2, p3 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(kind_phys) :: p0(imx), p1(imx), p2(imx), p3(imx)
      complex(kind_phys) :: cx(4,imx)

      integer :: i
      real(kind_phys) :: q(imx), r(imx)
      complex(kind_phys) :: cb(imx), cb0(imx), cb1(imx),   &
                     crad(imx), cy(imx), czero


      czero=cmplx(0.0_kind_phys,0.0_kind_phys,kind_phys)

      do 10 i=1,im

      q(i)=-p2(i)*p2(i)/36._kind_phys+(p3(i)*p1(i)-4*p0(i))/12._kind_phys
      r(i)=-(p2(i)/6)**3+p2(i)*(p3(i)*p1(i)-4*p0(i))/48._kind_phys   &
       +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16

      crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
      crad(i)=sqrt(crad(i))

      cb(i)=r(i)-crad(i)
      if(cb(i).eq.czero)then
!        insoluble particle
         cx(1,i)=(-p1(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
         cx(4,i)=cx(1,i)
      else
         cb(i)=cb(i)**third

         cy(i)=-cb(i)+q(i)/cb(i)+p2(i)/6

         cb0(i)=sqrt(cy(i)*cy(i)-p0(i))
         cb1(i)=(p3(i)*cy(i)-p1(i))/(2*cb0(i))

         cb(i)=p3(i)/2+cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)+cb0(i))
         crad(i)=sqrt(crad(i))
         cx(1,i)=(-cb(i)+crad(i))/2._kind_phys
         cx(2,i)=(-cb(i)-crad(i))/2._kind_phys

         cb(i)=p3(i)/2-cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)-cb0(i))
         crad(i)=sqrt(crad(i))
         cx(3,i)=(-cb(i)+crad(i))/2._kind_phys
         cx(4,i)=(-cb(i)-crad(i))/2._kind_phys
      endif
   10 continue

end subroutine makoh_quartic

!----------------------------------------------------------------------
subroutine calc_h2so4_equilib_mixrat( temp, pres, qh2o, dmean, &
                                      qh2so4_equilib, wtpct, sulden, &
                                      errmsg, errflg )

      real(kind_phys), intent(in)  :: temp            ! temperature (K)
      real(kind_phys), intent(in)  :: pres            ! pressure (Pa)
      real(kind_phys), intent(in)  :: qh2o            ! water vapor specific humidity (kg/kg)
      real(kind_phys), intent(in)  :: dmean           ! mean diameter of particles in each mode
      real(kind_phys), intent(out) :: qh2so4_equilib  ! h2so4 saturation mixing ratios over the particles (mol/mol)
      real(kind_phys), intent(out) :: wtpct           ! sulfate composition, weight % H2SO4
      real(kind_phys), intent(out) :: sulden          ! sulfate aerosol mass density (g/cm3)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local declarations
      real(kind_phys)            :: qh2o_kelvin ! water vapor specific humidity adjusted for Kelvin effect (kg/kg)
      real(kind_phys)            :: wtpct_flat  ! sulfate composition over a flat surface, weight % H2SO4
      real(kind_phys)            :: fk1, fk4, fk4_1, fk4_2
      real(kind_phys)            :: factor_kulm                     ! Kulmala correction terms
      real(kind_phys)            :: en, t, sig1, sig2, frac, surf_tens, surf_tens_mode, dsigma_dwt
      real(kind_phys)            :: den1, den2, sulfate_density, drho_dwt
      real(kind_phys)            :: akelvin, expon, akas, rkelvinH2O, rkelvinH2O_a, rkelvinH2O_b
      real(kind_phys)            :: sulfequil, r
      real(kind_phys), parameter :: t0_kulm     = 340._kind_phys           !  T0 set in the low end of the Ayers measurement range (338-445K)
      real(kind_phys), parameter :: t_crit_kulm = 905._kind_phys           !  Critical temperature = 1.5 * Boiling point
      real(kind_phys), parameter :: fk0         = -10156._kind_phys / t0_kulm + 16.259_kind_phys   !  Log(Kulmala correction factor)
      real(kind_phys), parameter :: fk2         = 1._kind_phys / t0_kulm
      real(kind_phys), parameter :: fk3         = 0.38_kind_phys / (t_crit_kulm - t0_kulm)
      real(kind_phys), parameter :: RGAS = 8.31430e+07_kind_phys           ! ideal gas constant (erg/mol/K)
      real(kind_phys), parameter :: wtmol_h2so4 =  98.078479_kind_phys     ! molecular weight of sulphuric acid
      real(kind_phys), parameter :: wtmol_h2o   =  18.015280_kind_phys     ! molecular weight of water vapor
      real(kind_phys), parameter :: wtmol_air   =  28.9644_kind_phys       ! molecular weight of dry air
      real(kind_phys)            :: gv_wt_abcd_en(6,95), gvh2ovp(95)
      real(kind_phys)            :: dnwtp(46), dnc0(46), dnc1(46)
      real(kind_phys)            :: stwtp(15), stc0(15), stc1(15)
      integer             :: i, k


      data stwtp/0._kind_phys, 23.8141_kind_phys, 38.0279_kind_phys, 40.6856_kind_phys, 45.335_kind_phys, 52.9305_kind_phys, &
         56.2735_kind_phys, 59.8557_kind_phys, 66.2364_kind_phys, 73.103_kind_phys, 79.432_kind_phys, 85.9195_kind_phys, &
         91.7444_kind_phys, 97.6687_kind_phys, 100._kind_phys/

      data stc0/117.564_kind_phys, 103.303_kind_phys, 101.796_kind_phys, 100.42_kind_phys, 98.4993_kind_phys, 91.8866_kind_phys, &
         88.3033_kind_phys, 86.5546_kind_phys, 84.471_kind_phys, 81.2939_kind_phys, 79.3556_kind_phys, 75.608_kind_phys, &
         70.0777_kind_phys, 63.7412_kind_phys, 61.4591_kind_phys /

      data stc1/-0.153641_kind_phys, -0.0982007_kind_phys, -0.0872379_kind_phys, -0.0818509_kind_phys,           &
         -0.0746702_kind_phys, -0.0522399_kind_phys, -0.0407773_kind_phys, -0.0357946_kind_phys, -0.0317062_kind_phys,   &
         -0.025825_kind_phys, -0.0267212_kind_phys, -0.0269204_kind_phys, -0.0276187_kind_phys, -0.0302094_kind_phys,    &
         -0.0303081_kind_phys /


      data dnwtp / 0._kind_phys, 1._kind_phys, 5._kind_phys, 10._kind_phys, 20._kind_phys, 25._kind_phys, 30._kind_phys, 35._kind_phys, 40._kind_phys,  &
         41._kind_phys, 45._kind_phys, 50._kind_phys, 53._kind_phys, 55._kind_phys, 56._kind_phys, 60._kind_phys, 65._kind_phys, 66._kind_phys, 70._kind_phys, &
         72._kind_phys, 73._kind_phys, 74._kind_phys, 75._kind_phys, 76._kind_phys, 78._kind_phys, 79._kind_phys, 80._kind_phys, 81._kind_phys, 82._kind_phys, &
         83._kind_phys, 84._kind_phys, 85._kind_phys, 86._kind_phys, 87._kind_phys, 88._kind_phys, 89._kind_phys, 90._kind_phys, 91._kind_phys, 92._kind_phys, &
         93._kind_phys, 94._kind_phys, 95._kind_phys, 96._kind_phys, 97._kind_phys, 98._kind_phys, 100._kind_phys /

      data dnc0 / 1._kind_phys, 1.13185_kind_phys, 1.17171_kind_phys, 1.22164_kind_phys, 1.3219_kind_phys, 1.37209_kind_phys,         &
         1.42185_kind_phys, 1.4705_kind_phys, 1.51767_kind_phys, 1.52731_kind_phys, 1.56584_kind_phys, 1.61834_kind_phys, 1.65191_kind_phys, &
         1.6752_kind_phys, 1.68708_kind_phys, 1.7356_kind_phys, 1.7997_kind_phys, 1.81271_kind_phys, 1.86696_kind_phys, 1.89491_kind_phys,   &
         1.9092_kind_phys, 1.92395_kind_phys, 1.93904_kind_phys, 1.95438_kind_phys, 1.98574_kind_phys, 2.00151_kind_phys, 2.01703_kind_phys, &
         2.03234_kind_phys, 2.04716_kind_phys, 2.06082_kind_phys, 2.07363_kind_phys, 2.08461_kind_phys, 2.09386_kind_phys, 2.10143_kind_phys,&
         2.10764_kind_phys, 2.11283_kind_phys, 2.11671_kind_phys, 2.11938_kind_phys, 2.12125_kind_phys, 2.1219_kind_phys, 2.12723_kind_phys, &
         2.12654_kind_phys, 2.12621_kind_phys, 2.12561_kind_phys, 2.12494_kind_phys, 2.12093_kind_phys /

      data dnc1 / 0._kind_phys,  -0.000435022_kind_phys, -0.000479481_kind_phys, -0.000531558_kind_phys, -0.000622448_kind_phys, &
         -0.000660866_kind_phys, -0.000693492_kind_phys, -0.000718251_kind_phys, -0.000732869_kind_phys, -0.000735755_kind_phys, &
         -0.000744294_kind_phys, -0.000761493_kind_phys, -0.000774238_kind_phys, -0.00078392_kind_phys, -0.000788939_kind_phys,  &
         -0.00080946_kind_phys, -0.000839848_kind_phys, -0.000845825_kind_phys, -0.000874337_kind_phys, -0.000890074_kind_phys,  &
         -0.00089873_kind_phys, -0.000908778_kind_phys, -0.000920012_kind_phys, -0.000932184_kind_phys, -0.000959514_kind_phys,  &
         -0.000974043_kind_phys, -0.000988264_kind_phys, -0.00100258_kind_phys, -0.00101634_kind_phys, -0.00102762_kind_phys,    &
         -0.00103757_kind_phys, -0.00104337_kind_phys, -0.00104563_kind_phys, -0.00104458_kind_phys, -0.00104144_kind_phys,      &
         -0.00103719_kind_phys, -0.00103089_kind_phys, -0.00102262_kind_phys, -0.00101355_kind_phys, -0.00100249_kind_phys,      &
         -0.00100934_kind_phys, -0.000998299_kind_phys, -0.000990961_kind_phys, -0.000985845_kind_phys, -0.000984529_kind_phys,  &
         -0.000989315_kind_phys /

      ! Saturation vapor pressure of sulfuric acid
      !
      ! Limit extrapolation at extreme temperatures
      t=min(max(temp,140._kind_phys),450._kind_phys)

      errmsg = ''
      errflg = 0

      !!  Calculate the weight % H2SO4 composition of sulfate
      call calc_h2so4_wtpct(t, pres, qh2o, wtpct_flat, errmsg, errflg)
      if (errflg /= 0) return

      !!  Calculate surface tension (erg/cm2) of sulfate of
      !!  different compositions as a linear function of temperature.
      i = 2 ! minimum wt% is 29.517
      do while (wtpct_flat.gt.stwtp(i))
       i = i + 1
      end do
      sig1 = stc0(i-1) + stc1(i-1) * t
      sig2 = stc0(i)   + stc1(i) * t
      ! calculate derivative needed later for kelvin factor for h2o
      dsigma_dwt = (sig2-sig1) / (stwtp(i)-stwtp(i-1))
      surf_tens = sig1 + dsigma_dwt*(wtpct_flat-stwtp(i))

      !!  Calculate density (g/cm3) of sulfate of
      !!  different compositions as a linear function of temperature.
      i = 6 ! minimum wt% is 29.517
      do while (wtpct_flat .gt. dnwtp(i))
        i = i + 1
      end do
      den1 = dnc0(i-1) + dnc1(i-1) * t
      den2 = dnc0(i)   + dnc1(i) * t
      ! Calculate derivative needed later for Kelvin factor for H2O
      drho_dwt = (den2-den1) / (dnwtp(i)-dnwtp(i-1))
      sulfate_density = den1 + drho_dwt * (wtpct_flat-dnwtp(i-1))

      r=dmean*100._kind_phys/2._kind_phys ! calcuate mode radius (cm) from diameter (m)

      ! Adjust for Kelvin effect for water
      rkelvinH2O_b = 1 + wtpct_flat * drho_dwt / &
         sulfate_density - 3._kind_phys * wtpct_flat * dsigma_dwt / (2._kind_phys*surf_tens)

      rkelvinH2O_a = 2._kind_phys * wtmol_h2so4 * surf_tens / &
           (sulfate_density * RGAS * t * r)

      rkelvinH2O = exp (rkelvinH2O_a*rkelvinH2O_b)

      qh2o_kelvin = qh2o/rkelvinH2O
      call calc_h2so4_wtpct(t, pres, qh2o_kelvin, wtpct, errmsg, errflg)
      if (errflg /= 0) return


      wtpct=max(wtpct,wtpct_flat)

      ! Parameterized fit to Giauque's (1959) enthalpies v. wt %:
      en = 4.184_kind_phys * (23624.8_kind_phys - 1.14208e8_kind_phys / ((wtpct - 105.318_kind_phys)**2 + 4798.69_kind_phys))
      en = max(en, 0.0_kind_phys)

      !!  Calculate surface tension (erg/cm2) of sulfate of
      !!  different compositions as a linear function of temperature.
      i = 2 ! minimum wt% is 29.517
      do while (wtpct.gt.stwtp(i))
       i=i+1
      end do
      sig1=stc0(i-1)+stc1(i-1)*t
      sig2=stc0(i)+stc1(i)*t
      frac=(stwtp(i)-wtpct)/(stwtp(i)-stwtp(i-1))
      surf_tens_mode=sig1*frac+sig2*(1.0_kind_phys-frac)

      !!  Calculate density (g/cm3) of sulfate of
      !!  different compositions as a linear function of temperature.
      i = 6 ! minimum wt% is 29.517
      do while (wtpct .gt. dnwtp(i))
        i=i+1
      end do
      den1=dnc0(i-1)+dnc1(i-1)*t
      den2=dnc0(i)+dnc1(i)*t
      frac=(dnwtp(i)-wtpct)/(dnwtp(i)-dnwtp(i-1))
      sulden=den1*frac+den2*(1.0_kind_phys-frac)

      ! Ayers' (1980) fit to sulfuric acid equilibrium vapor pressure:
      ! (Remember this is the log)
      ! SULFEQ=-10156/Temp+16.259-En/(8.3143*Temp)
      !
      ! Kulmala correction (J. CHEM. PHYS. V.93, No.1, 1 July 1990)
      fk1   = -1._kind_phys / t
      fk4_1 = log(t0_kulm / t)
      fk4_2 = t0_kulm / t
      fk4   = 1.0_kind_phys + fk4_1 - fk4_2
      factor_kulm = fk1 + fk2 + fk3 * fk4

      ! This is for pure H2SO4
      sulfequil = fk0 + 10156._kind_phys * factor_kulm

      ! Adjust for WTPCT composition:
      sulfequil = sulfequil - en / (8.3143_kind_phys * t)

      ! Take the exponential:
      sulfequil = exp(sulfequil)

      ! Convert atmospheres ==> Pa
      sulfequil = sulfequil * 1.01325e5_kind_phys

      ! Convert Pa ==> mol/mol
      sulfequil = sulfequil / pres

      ! Calculate Kelvin curvature factor for H2SO4 interactively with temperature:
      ! (g/mol)*(erg/cm2)/(K * g/cm3 * erg/mol/K) = cm
      akelvin = 2._kind_phys * wtmol_h2so4 * surf_tens_mode / (t * sulden * RGAS)

      expon = akelvin / r  ! divide by mode radius (cm)
      expon = max(-100._kind_phys, expon)
      expon = min(100._kind_phys, expon)
      akas = exp( expon )
      qh2so4_equilib = sulfequil * akas ! reduce H2SO4 equilibrium mixing ratio by Kelvin curvature factor

end subroutine calc_h2so4_equilib_mixrat

!----------------------------------------------------------------------
subroutine calc_h2so4_wtpct( temp, pres, qh2o, wtpct, errmsg, errflg )

  !!  This function calculates the weight % H2SO4 composition of
  !!  sulfate aerosol, using Tabazadeh et. al. (GRL, 1931, 1997).
  !!  Rated for T=185-260K, activity=0.01-1.0
  !!
  !!  Argument list input:
  !!    temp = temperature (K)
  !!    pres = atmospheric pressure (Pa)
  !!    qh2o = water specific humidity (kg/kg)
  !!
  !!  Output:
  !!    wtpct = weight % H2SO4 in H2O/H2SO4 particle (0-100)
  !!
  !! @author Mike Mills
  !! @ version October 2013

      use wv_saturation, only: qsat_water

      real(kind_phys), intent(in)  :: temp  ! temperature (K)
      real(kind_phys), intent(in)  :: pres  ! pressure (Pa)
      real(kind_phys), intent(in)  :: qh2o  ! water vapor specific humidity (kg/kg)
      real(kind_phys), intent(out) :: wtpct ! sulfate weight % H2SO4 composition
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local declarations
      real(kind_phys)            :: atab1,btab1,ctab1,dtab1,atab2,btab2,ctab2,dtab2
      real(kind_phys)            :: contl, conth, contt, conwtp
      real(kind_phys)            :: activ
      real(kind_phys)            :: es ! saturation vapor pressure over water (Pa) (dummy)
      real(kind_phys)            :: qs ! saturation specific humidity over water (kg/kg)

      errmsg = ''
      errflg = 0

      ! calculate saturation specific humidity over pure water, qs (kg/kg)
      call qsat_water(temp, pres, es, qs)

      !  Activity = water specific humidity (kg/kg) / equilibrium water (kg/kg)
      activ = qh2o/qs

      if (activ.lt.0.05_kind_phys) then
        activ = max(activ,1.e-6_kind_phys)    ! restrict minimum activity
        atab1 	= 12.37208932_kind_phys
        btab1 	= -0.16125516114_kind_phys
        ctab1 	= -30.490657554_kind_phys
        dtab1 	= -2.1133114241_kind_phys
        atab2 	= 13.455394705_kind_phys
        btab2 	= -0.1921312255_kind_phys
        ctab2 	= -34.285174607_kind_phys
        dtab2 	= -1.7620073078_kind_phys
      elseif (activ.ge.0.05_kind_phys.and.activ.le.0.85_kind_phys) then
        atab1 	= 11.820654354_kind_phys
        btab1 	= -0.20786404244_kind_phys
        ctab1 	= -4.807306373_kind_phys
        dtab1 	= -5.1727540348_kind_phys
        atab2 	= 12.891938068_kind_phys
        btab2 	= -0.23233847708_kind_phys
        ctab2 	= -6.4261237757_kind_phys
        dtab2 	= -4.9005471319_kind_phys
      elseif (activ.gt.0.85_kind_phys) then
        activ = min(activ,1._kind_phys)      ! restrict maximum activity
        atab1 	= -180.06541028_kind_phys
        btab1 	= -0.38601102592_kind_phys
        ctab1 	= -93.317846778_kind_phys
        dtab1 	= 273.88132245_kind_phys
        atab2 	= -176.95814097_kind_phys
        btab2 	= -0.36257048154_kind_phys
        ctab2 	= -90.469744201_kind_phys
        dtab2 	= 267.45509988_kind_phys
      else
        write(*,*) 'calc_h2so4_wtpct: invalid activity: activ,qh2o,qs,temp,pres=',activ,qh2o,qs,temp,pres
        errmsg = 'calc_h2so4_wtpct: invalid activity'
        errflg = 1
        return
      endif

      contl = atab1*(activ**btab1)+ctab1*activ+dtab1
      conth = atab2*(activ**btab2)+ctab2*activ+dtab2

      contt = contl + (conth-contl) * ((temp -190._kind_phys)/70._kind_phys)
      conwtp = (contt*98._kind_phys) + 1000._kind_phys

      wtpct = (100._kind_phys*contt*98._kind_phys)/conwtp
      wtpct = min(max(wtpct,25._kind_phys),100._kind_phys) ! restrict between 1 and 100 %

end subroutine calc_h2so4_wtpct

!----------------------------------------------------------------------

subroutine modal_aero_wateruptake_diag( &
   aero_props, aero_state, &
   ncol, nlev, top_lev, &
   pi, rhoh2o, &
   t, pmid, h2ommr, cldn, &
   bin_idx, dgnumwet, qaerwat, &
   errmsg, errflg)

!-----------------------------------------------------------------------
!
! Recompute wet number mode diameter and aerosol water for a DIAGNOSTIC
! radiation list, returning the slices for one mode. Composes the portable
! diagnostic-list size calculation (modal_aero_calcsize_diag_run +
! modal_aero_calcdry_run) with the wet radius calculation
! (modal_aero_wateruptake_sub) and the wet-diameter/water/density
! post-processing of the climate-list driver.
!
! The stratospheric sulfate treatment is not supported for diagnostic
! lists (matching the climate-list driver, which aborts in that case),
! so do_strat_sulfate is hardwired false and troplev is unused.
!
! This routine is the target of the water-uptake-diagnostic procedure
! pointer registered with modal_aerosol_state; it must keep the plain
! abstract-interface signature so hosts without the modal aerosol
! schemes never reference this module.
!
!-----------------------------------------------------------------------

   use aerosol_properties_mod, only: aerosol_properties
   use aerosol_state_mod,      only: aerosol_state
   use modal_aero_calcsize,    only: modal_aero_calcsize_diag_run, modal_aero_calcdry_run

   ! Arguments
   class(aerosol_properties), intent(in) :: aero_props
   class(aerosol_state),      intent(in) :: aero_state
   integer,          intent(in)  :: ncol              ! number of columns
   integer,          intent(in)  :: nlev              ! number of vertical levels
   integer,          intent(in)  :: top_lev           ! top level for aerosol calculations
   real(kind_phys),         intent(in)  :: pi                ! pi
   real(kind_phys),         intent(in)  :: rhoh2o            ! density of liquid water (kg/m3)
   real(kind_phys),         intent(in)  :: t(:,:)            ! temperature (K)
   real(kind_phys),         intent(in)  :: pmid(:,:)         ! layer pressure (Pa)
   real(kind_phys),         intent(in)  :: h2ommr(:,:)       ! specific humidity (kg/kg)
   real(kind_phys),         intent(in)  :: cldn(:,:)         ! layer cloud fraction (0-1)
   integer,          intent(in)  :: bin_idx           ! mode index of the returned slices
   real(kind_phys),         intent(out) :: dgnumwet(:,:)     ! wet number mode diameter of mode bin_idx (m)
   real(kind_phys),         intent(out) :: qaerwat(:,:)      ! aerosol water of mode bin_idx (g/g)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! local variables
   integer :: i, k, m
   integer :: nmodes
   integer :: istat

   integer :: troplev(ncol)

   real(kind_phys), allocatable :: dgncur_a(:,:,:)    ! dry number mode diameter (m)
   real(kind_phys), allocatable :: dgncur_awet(:,:,:) ! wet number mode diameter (m)
   real(kind_phys), allocatable :: qaerwat_m(:,:,:)   ! aerosol water (g/g)
   real(kind_phys), allocatable :: wetdens(:,:,:)     ! wet aerosol density (kg/m3)
   real(kind_phys), allocatable :: hygro(:,:,:)       ! volume-weighted mean hygroscopicity (--)
   real(kind_phys), allocatable :: naer(:,:,:)        ! aerosol number MR (bounded!) (#/kg-air)
   real(kind_phys), allocatable :: dryvol(:,:,:)      ! single-particle-mean dry volume (m3)
   real(kind_phys), allocatable :: so4dryvol(:,:,:)   ! single-particle-mean so4 dry volume (m3)
   real(kind_phys), allocatable :: drymass(:,:,:)     ! single-particle-mean dry mass  (kg)
   real(kind_phys), allocatable :: dryrad(:,:,:)      ! dry volume mean radius of aerosol (m)

   real(kind_phys), allocatable :: wetrad(:,:,:)      ! wet radius of aerosol (m)
   real(kind_phys), allocatable :: wetvol(:,:,:)      ! single-particle-mean wet volume (m3)
   real(kind_phys), allocatable :: wtrvol(:,:,:)      ! single-particle-mean water volume in wet aerosol (m3)

   real(kind_phys), allocatable :: sulfeq(:,:,:)      ! H2SO4 equilibrium mixing ratios over particles (mol/mol)
   real(kind_phys), allocatable :: wtpct(:,:,:)       ! sulfate aerosol composition, weight % H2SO4
   real(kind_phys), allocatable :: sulden(:,:,:)      ! sulfate aerosol mass density (g/cm3)

   real(kind_phys), allocatable :: specdens_1(:)
   real(kind_phys), allocatable :: alnsg(:)
   real(kind_phys), allocatable :: maer(:,:,:)        ! accumulated aerosol mode MRs

   !-----------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   nmodes = aero_props%nbins()

   allocate( &
      dgncur_a(ncol,nlev,nmodes),    dgncur_awet(ncol,nlev,nmodes), &
      qaerwat_m(ncol,nlev,nmodes),   wetdens(ncol,nlev,nmodes),     &
      hygro(ncol,nlev,nmodes),       dryvol(ncol,nlev,nmodes),      &
      dryrad(ncol,nlev,nmodes),      drymass(ncol,nlev,nmodes),     &
      so4dryvol(ncol,nlev,nmodes),   naer(ncol,nlev,nmodes),        &
      wetrad(ncol,nlev,nmodes),      wetvol(ncol,nlev,nmodes),      &
      wtrvol(ncol,nlev,nmodes),      wtpct(ncol,nlev,nmodes),       &
      sulden(ncol,nlev,nmodes),      sulfeq(ncol,nlev,nmodes),      &
      specdens_1(nmodes),            alnsg(nmodes),                 &
      maer(ncol,nlev,nmodes),        stat=istat)
   if (istat > 0) then
      errmsg = 'modal_aero_wateruptake_diag: unable to allocate work arrays'
      errflg = 1
      return
   end if

   ! dry size distribution parameters of the diagnostic list
   call modal_aero_calcsize_diag_run( &
      aero_props = aero_props,       &
      aero_state = aero_state,       &
      ncol       = ncol,             &
      pver       = nlev,             &
      top_lev    = top_lev,          &
      pi         = pi,               &
      dgncur_a   = dgncur_a,         &
      errmsg     = errmsg,           &
      errflg     = errflg)
   if (errflg /= 0) return

   ! Zero output fields (_run writes top_lev:nlev)
   hygro(:,:,:)     = 0._kind_phys
   dryvol(:,:,:)    = 0._kind_phys
   dryrad(:,:,:)    = 0._kind_phys
   drymass(:,:,:)   = 0._kind_phys
   so4dryvol(:,:,:) = 0._kind_phys
   naer(:,:,:)      = 0._kind_phys

   call modal_aero_calcdry_run( &
      aero_props       = aero_props,   &
      aero_state       = aero_state,   &
      ncol             = ncol,         &
      pver             = nlev,         &
      top_lev          = top_lev,      &
      do_strat_sulfate = .false.,      &
      pi               = pi,           &
      dgncur_a         = dgncur_a,     &
      hygro            = hygro,        &
      dryvol           = dryvol,       &
      dryrad           = dryrad,       &
      drymass          = drymass,      &
      so4dryvol        = so4dryvol,    &
      naer             = naer,         &
      errmsg           = errmsg,       &
      errflg           = errflg)
   if (errflg /= 0) return

   ! Zero work arrays (_sub only writes top_lev:nlev).
   ! dgncur_awet is an intent(in) of _sub read only under do_strat_sulfate;
   ! zero it here for definedness before the post-processing fills it.
   wetrad(:,:,:)      = 0._kind_phys
   wetvol(:,:,:)      = 0._kind_phys
   wtrvol(:,:,:)      = 0._kind_phys
   sulfeq(:,:,:)      = 0._kind_phys
   wtpct(:,:,:)       = 0._kind_phys
   sulden(:,:,:)      = 0._kind_phys
   maer(:,:,:)        = 0._kind_phys
   dgncur_awet(:,:,:) = 0._kind_phys
   troplev(:)         = 0

   call modal_aero_wateruptake_sub( &
      aero_props       = aero_props,   &
      aero_state       = aero_state,   &
      ncol             = ncol,         &
      pver             = nlev,         &
      top_lev          = top_lev,      &
      do_strat_sulfate = .false.,      &
      t                = t,            &
      pmid             = pmid,         &
      h2ommr           = h2ommr,       &
      cldn             = cldn,         &
      dryrad           = dryrad,       &
      hygro            = hygro,        &
      dryvol           = dryvol,       &
      so4dryvol        = so4dryvol,    &
      dgncur_awet      = dgncur_awet,  &
      troplev          = troplev,      &
      wetrad           = wetrad,       &
      wetvol           = wetvol,       &
      wtrvol           = wtrvol,       &
      sulfeq           = sulfeq,       &
      wtpct            = wtpct,        &
      sulden           = sulden,       &
      specdens_1       = specdens_1,   &
      alnsg_out        = alnsg,        &
      maer             = maer,         &
      errmsg           = errmsg,       &
      errflg           = errflg)
   if (errflg /= 0) return

   ! Post-processing: wet density, qaerwat, dgncur_awet update
   qaerwat_m = 0.0_kind_phys

   do m = 1, nmodes

      do k = top_lev, nlev
         do i = 1, ncol

            dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m))
            qaerwat_m(i,k,m)   = rhoh2o*naer(i,k,m)*wtrvol(i,k,m)

            ! compute aerosol wet density (kg/m3)
            if (wetvol(i,k,m) > 1.0e-30_kind_phys) then
               wetdens(i,k,m) = (drymass(i,k,m) + rhoh2o*wtrvol(i,k,m))/wetvol(i,k,m)
            else
               wetdens(i,k,m) = specdens_1(m)
            end if
         end do
      end do

   end do    ! modes

   dgnumwet(:ncol,:nlev) = dgncur_awet(:ncol,:nlev,bin_idx)
   qaerwat (:ncol,:nlev) =   qaerwat_m(:ncol,:nlev,bin_idx)

   deallocate( &
      dgncur_a, dgncur_awet, qaerwat_m, wetdens, hygro, dryvol, dryrad, &
      drymass, so4dryvol, naer, wetrad, wetvol, wtrvol, wtpct, sulden,  &
      sulfeq, specdens_1, alnsg, maer)

end subroutine modal_aero_wateruptake_diag

!----------------------------------------------------------------------

end module modal_aero_wateruptake
