module aero_activate

! Portable (CCPP-ready) Abdul-Razzak & Ghan aerosol activation kernel
! (activate_aerosol), extracted from the CAM ndrop module.  The derived
! constants aten and sqrt(pi) are computed once by aero_activate_init from host
! physical constants; activate_aerosol receives the remaining host physical
! constants as arguments and returns errmsg/errflg instead of aborting.  The
! polymorphic aerosol_properties abstraction is deliberately host-portable.
! CAM interface / callers: ndrop, aero_convproc.

use ccpp_kinds, only: kind_phys
use wv_saturation,    only: qsat
use shr_spfn_mod,     only: erf => shr_spfn_erf

use aerosol_properties_mod, only: aerosol_properties

implicit none
private

public :: aero_activate_init, activate_aerosol

! mathematical constants
real(kind_phys), parameter :: zero     = 0._kind_phys
real(kind_phys), parameter :: third    = 1._kind_phys/3._kind_phys
real(kind_phys), parameter :: twothird = 2._kind_phys*third
real(kind_phys), parameter :: sixth    = 1._kind_phys/6._kind_phys
real(kind_phys), parameter :: sq2      = sqrt(2._kind_phys)
real(kind_phys), parameter :: tmelt    = 273._kind_phys

! derived constants (set by aero_activate_init)
real(kind_phys) :: aten
real(kind_phys) :: sqpi

!===============================================================================
contains
!===============================================================================

subroutine aero_activate_init(mwh2o, r_universal, rhoh2o, pi)

   ! Compute the derived activation constants from host physical constants.
   ! surften (surface tension of water) is a fixed property of the activation
   ! parameterization, so it is set here rather than threaded from the host.

   real(kind_phys), intent(in) :: mwh2o        ! molecular weight of water (kg/kmol)
   real(kind_phys), intent(in) :: r_universal  ! universal gas constant (J/K/kmol)
   real(kind_phys), intent(in) :: rhoh2o       ! density of liquid water (kg/m3)
   real(kind_phys), intent(in) :: pi           ! pi

   real(kind_phys) :: surften                  ! surface tension of water w/respect to air (N/m)

   surften = 0.076_kind_phys
   aten = 2._kind_phys*mwh2o*surften/(r_universal*tmelt*rhoh2o)
   sqpi = sqrt(pi)

end subroutine aero_activate_init

!===============================================================================

subroutine activate_aerosol(wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,  &
   na, nbins, volume, hygro, aero_props, &
   fn, fm, fluxn, fluxm, flux_fullact, &
   pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, errmsg, errflg, &
   smax_prescribed, in_cloud_in, smax_f)

   !      calculates number, surface, and mass fraction of aerosols activated as CCN
   !      calculates flux of cloud droplets, surface area, and aerosol mass into cloud
   !      assumes an internal mixture within each of up to nbin multiple aerosol bins
   !      a gaussiam spectrum of updrafts can be treated.

   !      mks units

   !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   !      input

   real(kind_phys), intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
   real(kind_phys), intent(in) :: sigw          ! subgrid standard deviation of vertical vel (m/s)
   real(kind_phys), intent(in) :: wdiab         ! diabatic vertical velocity (0 if adiabatic)
   real(kind_phys), intent(in) :: wminf         ! minimum updraft velocity for integration (m/s)
   real(kind_phys), intent(in) :: wmaxf         ! maximum updraft velocity for integration (m/s)
   real(kind_phys), intent(in) :: tair          ! air temperature (K)
   real(kind_phys), intent(in) :: rhoair        ! air density (kg/m3)
   real(kind_phys), intent(in) :: na(:)      ! aerosol number concentration (/m3)
   integer,  intent(in) :: nbins      ! number of aerosol bins
   real(kind_phys), intent(in) :: volume(:)  ! aerosol volume concentration (m3/m3)
   real(kind_phys), intent(in) :: hygro(:)   ! hygroscopicity of aerosol mode

   class(aerosol_properties), intent(in) :: aero_props

   !      output

   real(kind_phys), intent(out) :: fn(:)      ! number fraction of aerosols activated
   real(kind_phys), intent(out) :: fm(:)      ! mass fraction of aerosols activated
   real(kind_phys), intent(out) :: fluxn(:)   ! flux of activated aerosol number fraction into cloud (cm/s)
   real(kind_phys), intent(out) :: fluxm(:)   ! flux of activated aerosol mass fraction into cloud (cm/s)
   real(kind_phys), intent(out) :: flux_fullact   ! flux of activated aerosol fraction assuming 100% activation (cm/s)
   !    rce-comment
   !    used for consistency check -- this should match (ekd(k)*zs(k))
   !    also, fluxm/flux_fullact gives fraction of aerosol mass flux
   !       that is activated

   !      host physical constants (from physconst) + error handling
   real(kind_phys), intent(in) :: pi           ! pi
   real(kind_phys), intent(in) :: rhoh2o       ! density of liquid water (kg/m3)
   real(kind_phys), intent(in) :: rh2o         ! water vapor gas constant (J/K/kg)
   real(kind_phys), intent(in) :: gravit       ! gravitational acceleration (m/s2)
   real(kind_phys), intent(in) :: latvap       ! latent heat of vaporization (J/kg)
   real(kind_phys), intent(in) :: cpair        ! specific heat of dry air (J/K/kg)
   real(kind_phys), intent(in) :: rair         ! dry air gas constant (J/K/kg)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   !      optional
   real(kind_phys), optional, intent(in) :: smax_prescribed  ! prescribed max. supersaturation for secondary activation
   logical,  optional, intent(in) :: in_cloud_in      ! switch to modify calculations when above cloud base
   real(kind_phys), optional, intent(in) :: smax_f           ! droplet and rain size distr factor in the smax calculation
                                                      ! used when in_cloud=.true.

   !      local

   integer, parameter :: nx=200
   real(kind_phys) :: integ,integf
   real(kind_phys), parameter :: p0 = 1013.25e2_kind_phys    ! reference pressure (Pa)
   real(kind_phys) :: pres ! pressure (Pa)
   real(kind_phys) :: diff0,conduct0
   real(kind_phys) :: es ! saturation vapor pressure
   real(kind_phys) :: qs ! water vapor saturation mixing ratio
   real(kind_phys) :: dqsdt ! change in qs with temperature
   real(kind_phys) :: g ! thermodynamic function (m2/s)
   real(kind_phys) :: zeta(nbins), eta(nbins)
   real(kind_phys) :: alpha
   real(kind_phys) :: gamma
   real(kind_phys) :: beta
   real(kind_phys) :: sqrtg
   real(kind_phys) :: amcube(nbins) ! cube of dry bin radius (m)
   real(kind_phys) :: smc(nbins) ! critical supersaturation for number bin radius
   real(kind_phys) :: sumflx_fullact
   real(kind_phys) :: sumflxn(nbins)
   real(kind_phys) :: sumflxm(nbins)
   real(kind_phys) :: sumfn(nbins)
   real(kind_phys) :: sumfm(nbins)
   real(kind_phys) :: fnold(nbins)   ! number fraction activated
   real(kind_phys) :: fmold(nbins)   ! mass fraction activated
   real(kind_phys) :: wold,gold
   real(kind_phys) :: wmin,wmax,w,dw,dwmax,dwmin,wnuc,dwnew,wb
   real(kind_phys) :: dfmin,dfmax,fnew,fold,fnmin,fnbar,fmbar
   real(kind_phys) :: alw,sqrtalw
   real(kind_phys) :: smax
   real(kind_phys) :: z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
   real(kind_phys) :: etafactor1,etafactor2(nbins),etafactor2max
   real(kind_phys) :: grow
   character(len=*), parameter :: subname='activate_aerosol'

   logical :: in_cloud
   integer :: m,n
   !      numerical integration parameters
   real(kind_phys), parameter :: eps=0.3_kind_phys,fmax=0.99_kind_phys,sds=3._kind_phys

   real(kind_phys), parameter :: namin=1.e6_kind_phys   ! minimum aerosol number concentration (/m3)

   errmsg = ''
   errflg = 0

   if (present(in_cloud_in)) then
      if (.not. present(smax_f)) then
         errmsg = subname//' error: smax_f must be supplied when in_cloud is used'
         errflg = 1
         return
      end if
      in_cloud = in_cloud_in
   else
      in_cloud = .false.
   end if

   fn(:)=0._kind_phys
   fm(:)=0._kind_phys
   fluxn(:)=0._kind_phys
   fluxm(:)=0._kind_phys
   flux_fullact=0._kind_phys

   if(nbins==1.and.na(1)<1.e-20_kind_phys)return

   if(sigw<=1.e-5_kind_phys.and.wbar<=0._kind_phys)return

   if (present(smax_prescribed)) then
      if (smax_prescribed <= 0.0_kind_phys) return
   end if

   pres=rair*rhoair*tair
   diff0=0.211e-4_kind_phys*(p0/pres)*(tair/tmelt)**1.94_kind_phys
   conduct0=(5.69_kind_phys+0.017_kind_phys*(tair-tmelt))*4.186e2_kind_phys*1.e-5_kind_phys ! convert to J/m/s/deg
   call qsat(tair, pres, es, qs)
   dqsdt=latvap/(rh2o*tair*tair)*qs
   alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1._kind_phys/(rair*tair))
   gamma=(1.0_kind_phys+latvap/cpair*dqsdt)/(rhoair*qs)
   etafactor2max=1.e10_kind_phys/(alpha*wmaxf)**1.5_kind_phys ! this should make eta big if na is very small.

   grow  = 1._kind_phys/(rhoh2o/(diff0*rhoair*qs)  &
           + latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair) - 1._kind_phys))
   sqrtg = sqrt(grow)
   beta  = 2._kind_phys*pi*rhoh2o*grow*gamma

   do m=1,nbins

      if(volume(m)>1.e-39_kind_phys.and.na(m)>1.e-39_kind_phys)then
         ! number mode radius (m)
         amcube(m)=aero_props%amcube(m, volume(m),na(m))
         ! growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
         ! should depend on mean radius of mode to account for gas kinetic effects
         ! see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
         ! for approriate size to use for effective diffusivity.
         etafactor2(m)=1._kind_phys/(na(m)*beta*sqrtg)
         if(hygro(m)>1.e-10_kind_phys)then
            smc(m)=2._kind_phys*aten*sqrt(aten/(27._kind_phys*hygro(m)*amcube(m))) ! only if variable size dist
         else
            smc(m)=100._kind_phys
         end if
      else
         smc(m)=1._kind_phys
         etafactor2(m)=etafactor2max ! this should make eta big if na is very small.
      end if

   end do

   if(sigw>1.e-5_kind_phys)then ! spectrum of updrafts

      wmax=min(wmaxf,wbar+sds*sigw)
      wmin=max(wminf,-wdiab)
      wmin=max(wmin,wbar-sds*sigw)
      w=wmin
      dwmax=eps*sigw
      dw=dwmax
      dfmax=0.2_kind_phys
      dfmin=0.1_kind_phys
      if (wmax <= w) return
      do m=1,nbins
         sumflxn(m)=0._kind_phys
         sumfn(m)=0._kind_phys
         fnold(m)=0._kind_phys
         sumflxm(m)=0._kind_phys
         sumfm(m)=0._kind_phys
         fmold(m)=0._kind_phys
      end do
      sumflx_fullact=0._kind_phys

      fold=0._kind_phys
      wold=0._kind_phys
      gold=0._kind_phys

      dwmin = min(dwmax, 0.01_kind_phys)
      do n = 1, nx

100      wnuc=w+wdiab
         !           write(iulog,*)'wnuc=',wnuc
         alw=alpha*wnuc
         sqrtalw=sqrt(alw)
         etafactor1=alw*sqrtalw

         do m=1,nbins
            eta(m)=etafactor1*etafactor2(m)
            zeta(m)=twothird*sqrtalw*aten/sqrtg
         end do

         if (present(smax_prescribed)) then
            smax = smax_prescribed
         else
            smax = aero_props%maxsat(zeta,eta,smc)
         end if

         call aero_props%actfracs(nbins, smc(nbins), smax, fnew, fm(nbins))

         dwnew = dw
         if(fnew-fold>dfmax.and.n>1)then
            !              reduce updraft increment for greater accuracy in integration
            if (dw > 1.01_kind_phys*dwmin) then
               dw=0.7_kind_phys*dw
               dw=max(dw,dwmin)
               w=wold+dw
               goto 100
            else
               dwnew = dwmin
            end if
         end if

         if(fnew-fold<dfmin)then
            !              increase updraft increment to accelerate integration
            dwnew=min(1.5_kind_phys*dw,dwmax)
         end if
         fold=fnew

         z=(w-wbar)/(sigw*sq2)
         g=exp(-z*z)
         fnmin=1._kind_phys

         do m=1,nbins
            !              modal
            call aero_props%actfracs(m, smc(m), smax, fn(m), fm(m))
            fnmin=min(fn(m),fnmin)
            !               integration is second order accurate
            !               assumes linear variation of f*g with w
            fnbar=(fn(m)*g+fnold(m)*gold)
            fmbar=(fm(m)*g+fmold(m)*gold)
            wb=(w+wold)
            if(w>0._kind_phys)then
               sumflxn(m)=sumflxn(m)+sixth*(wb*fnbar           &
                  +(fn(m)*g*w+fnold(m)*gold*wold))*dw
               sumflxm(m)=sumflxm(m)+sixth*(wb*fmbar           &
                  +(fm(m)*g*w+fmold(m)*gold*wold))*dw
            end if
            sumfn(m)=sumfn(m)+0.5_kind_phys*fnbar*dw
            fnold(m)=fn(m)
            sumfm(m)=sumfm(m)+0.5_kind_phys*fmbar*dw
            fmold(m)=fm(m)
         end do
         !           same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact &
            + sixth*(wb*(g+gold) + (g*w+gold*wold))*dw
         gold=g
         wold=w
         dw=dwnew
         if (n > 1 .and. (w > wmax .or. fnmin > fmax)) exit
         w=w+dw
         if (n == nx) then
            errmsg = subname//' -- do loop is too short in activate'
            errflg = 1
            return
         end if

      end do

      if(w<wmaxf)then

         !            contribution from all updrafts stronger than wmax
         !            assuming constant f (close to fmax)
         wnuc=w+wdiab

         z1=(w-wbar)/(sigw*sq2)
         z2=(wmaxf-wbar)/(sigw*sq2)
         g=exp(-z1*z1)
         integ=sigw*0.5_kind_phys*sq2*sqpi*(erf(z2)-erf(z1))
         !            consider only upward flow into cloud base when estimating flux
         wf1=max(w,zero)
         zf1=(wf1-wbar)/(sigw*sq2)
         gf1=exp(-zf1*zf1)
         wf2=max(wmaxf,zero)
         zf2=(wf2-wbar)/(sigw*sq2)
         gf2=exp(-zf2*zf2)
         gf=(gf1-gf2)
         integf=wbar*sigw*0.5_kind_phys*sq2*sqpi*(erf(zf2)-erf(zf1))+sigw*sigw*gf

         do m=1,nbins
            sumflxn(m)=sumflxn(m)+integf*fn(m)
            sumfn(m)=sumfn(m)+fn(m)*integ
            sumflxm(m)=sumflxm(m)+integf*fm(m)
            sumfm(m)=sumfm(m)+fm(m)*integ
         end do
         !           same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact + integf
         !            sumg=sumg+integ
      end if


      do m=1,nbins
         fn(m)=sumfn(m)/(sq2*sqpi*sigw)
         !            fn(m)=sumfn(m)/(sumg)
         if(fn(m)>1.01_kind_phys)then
            errmsg = 'activate -- fn > 1'
            errflg = 1
            return
         end if
         fluxn(m)=sumflxn(m)/(sq2*sqpi*sigw)
         fm(m)=sumfm(m)/(sq2*sqpi*sigw)
         !            fm(m)=sumfm(m)/(sumg)
         fluxm(m)=sumflxm(m)/(sq2*sqpi*sigw)
      end do
      !        same form as fluxm
      flux_fullact = sumflx_fullact/(sq2*sqpi*sigw)

   else

      !        single updraft
      wnuc=wbar+wdiab

      if(wnuc>0._kind_phys)then

         w=wbar

         if(in_cloud) then

            if (smax_f > 0._kind_phys) then
               smax = alpha*w/(2.0_kind_phys*pi*rhoh2o*grow*gamma*smax_f)
            else
               smax = 1.e-20_kind_phys
            end if

         else ! at cloud base
            alw        = alpha*wnuc
            sqrtalw    = sqrt(alw)
            etafactor1 = alw*sqrtalw

            do m = 1, nbins
               eta(m)  = etafactor1*etafactor2(m)
               zeta(m) = twothird*sqrtalw*aten/sqrtg
            end do
            if (present(smax_prescribed)) then
               smax = smax_prescribed
            else
               smax = aero_props%maxsat(zeta,eta,smc)
            end if
         end if

         do m=1,nbins

            call aero_props%actfracs(m, smc(m), smax, fn(m), fm(m))

            if(wbar>0._kind_phys)then
               fluxn(m)=fn(m)*w
               fluxm(m)=fm(m)*w
            end if
         end do
         flux_fullact = w
      end if

   end if

end subroutine activate_aerosol

end module aero_activate
