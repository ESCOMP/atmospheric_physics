module nucleate_ice

!-------------------------------------------------------------------------------
! Purpose:
!  A parameterization of ice nucleation.
!
!  *** This module is intended to be a "portable" code layer.  Ideally it should
!  *** not contain any use association of modules that belong to the model framework.
!
!
! Method:
!  The current method is based on Liu & Penner (2005) & Liu et al. (2007)
!  It related the ice nucleation with the aerosol number, temperature and the
!  updraft velocity. It includes homogeneous freezing of sulfate & immersion
!  freezing on mineral dust (soot disabled) in cirrus clouds, and
!  Meyers et al. (1992) deposition nucleation in mixed-phase clouds
!
!  The effect of preexisting ice crystals on ice nucleation in cirrus clouds is included,
!  and also consider the sub-grid variability of temperature in cirrus clouds,
!  following X. Shi et al. ACP (2014).
!
!  Ice nucleation in mixed-phase clouds now uses classical nucleation theory (CNT),
!  follows Y. Wang et al. ACP (2014), Hoose et al. (2010).
!
! References:
!  Liu & Penner, 2005: https://doi.org/10.1127/0941-2948/2005/0059
!  Liu et al., 2007:   https://doi.org/10.1175/JCLI4264.1
!  Shi et al., 2014:   https://doi.org/10.5194/acp-15-1503-2015
!  Wang et al., 2014:  https://doi.org/10.5194/acp-14-10411-2014
!  Hoose et al., 2010: https://doi.org/10.1175/2010JAS3425.1
!  Meyers et al., 1992: https://doi.org/10.1175/1520-0450(1992)031<0708:NPINPI>2.0.CO;2
!
! Authors:
!  Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
!  Xiangjun Shi & Xiaohong Liu, 01/2014.
!
!  With help from C. C. Chen and B. Eaton (2014)
!-------------------------------------------------------------------------------

use ccpp_kinds,     only: kind_phys

implicit none
private

public :: nucleati_init
public :: nucleati

logical  :: use_preexisting_ice
logical  :: use_hetfrz_classnuc
logical  :: use_incloud_nuc
integer  :: iulog
real(kind_phys) :: pi
real(kind_phys) :: mincld

real(kind_phys), parameter :: Shet   = 1.3_kind_phys     ! het freezing threshold
real(kind_phys), parameter :: rhoice = 0.5e3_kind_phys   ! assumed cloud ice density (kg/m3); Wpice is not sensitive to rhoice
real(kind_phys), parameter :: minweff= 0.001_kind_phys   ! minimum effective vertical velocity for ice nucleation (m/s)
! Gamma-function values Gamma(n) for the assumed exponential (mu=0) ice size
! distribution. Only gamma4 is currently used (in the lami calculation below).
real(kind_phys), parameter :: gamma1=1.0_kind_phys       ! Gamma(1) = 0!
real(kind_phys), parameter :: gamma2=1.0_kind_phys       ! Gamma(2) = 1!
real(kind_phys), parameter :: gamma3=2.0_kind_phys       ! Gamma(3) = 2!
real(kind_phys), parameter :: gamma4=6.0_kind_phys       ! Gamma(4) = 3!

real(kind_phys) :: ci

!===============================================================================
contains
!===============================================================================

subroutine nucleati_init( &
   use_preexisting_ice_in, use_hetfrz_classnuc_in, use_incloud_nuc_in, iulog_in, pi_in, &
   mincld_in)

   logical,         intent(in) :: use_preexisting_ice_in
   logical,         intent(in) :: use_hetfrz_classnuc_in
   logical,         intent(in) :: use_incloud_nuc_in
   integer,         intent(in) :: iulog_in
   real(kind_phys), intent(in) :: pi_in
   real(kind_phys), intent(in) :: mincld_in

   use_preexisting_ice = use_preexisting_ice_in
   use_hetfrz_classnuc = use_hetfrz_classnuc_in
   use_incloud_nuc     = use_incloud_nuc_in
   iulog               = iulog_in
   pi                  = pi_in
   mincld              = mincld_in

   ci = rhoice*pi/6._kind_phys

end subroutine nucleati_init

!===============================================================================

subroutine nucleati(  &
   wbar, tair, pmid, relhum, cldn,      &
   qc, qi, ni_in, rhoair,               &
   so4_num, dst_num, soot_num, subgrid, &
   nuci, onihf, oniimm, onidep, onimey, &
   wpice, weff, fhom, regm, &
   oso4_num, odst_num, osoot_num, &
   call_frm_zm_in, add_preexisting_ice_in)

   use wv_saturation, only: svp_water, svp_ice

   ! Input Arguments
   real(kind_phys), intent(in) :: wbar        ! grid cell mean vertical velocity (m/s)
   real(kind_phys), intent(in) :: tair        ! temperature (K)
   real(kind_phys), intent(in) :: pmid        ! pressure at layer midpoints (pa)
   real(kind_phys), intent(in) :: relhum      ! relative humidity with respective to liquid
   real(kind_phys), intent(in) :: cldn        ! new value of cloud fraction    (fraction)
   real(kind_phys), intent(in) :: qc          ! liquid water mixing ratio (kg/kg)
   real(kind_phys), intent(in) :: qi          ! grid-mean preexisting cloud ice mass mixing ratio (kg/kg)
   real(kind_phys), intent(in) :: ni_in       ! grid-mean preexisting cloud ice number conc (#/kg)
   real(kind_phys), intent(in) :: rhoair      ! air density (kg/m3)
   real(kind_phys), intent(in) :: so4_num     ! so4 aerosol number (#/cm^3)
   real(kind_phys), intent(in) :: dst_num     ! total dust aerosol number (#/cm^3)
   real(kind_phys), intent(in) :: soot_num    ! soot (hydrophilic) aerosol number (#/cm^3)
   real(kind_phys), intent(in) :: subgrid     ! subgrid saturation scaling factor

   ! Output Arguments
   real(kind_phys), intent(out) :: nuci       ! ice number nucleated (#/kg)
   real(kind_phys), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4
   real(kind_phys), intent(out) :: oniimm     ! nucleated number from immersion freezing
   real(kind_phys), intent(out) :: onidep     ! nucleated number from deposition nucleation
   real(kind_phys), intent(out) :: onimey     ! nucleated number from deposition nucleation  (meyers: mixed phase)
   real(kind_phys), intent(out) :: wpice      ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
   real(kind_phys), intent(out) :: weff       ! effective Vertical velocity for ice nucleation (m/s); weff=wbar-wpice
   real(kind_phys), intent(out) :: fhom       ! how much fraction of cloud can reach Shom
   real(kind_phys), intent(out) :: regm       ! nucleation regime indiator
   real(kind_phys), intent(out) :: oso4_num   ! so4 aerosol number (#/cm^3)
   real(kind_phys), intent(out) :: odst_num   ! total dust aerosol number (#/cm^3)
   real(kind_phys), intent(out) :: osoot_num  ! soot (hydrophilic) aerosol number (#/cm^3)

   ! Optional Arguments
   logical,  intent(in), optional :: call_frm_zm_in ! true if called from ZM convection scheme
   logical,  intent(in), optional :: add_preexisting_ice_in ! only false if called with pumas_v1.21+

   ! Local workspace
   real(kind_phys) :: nihf                      ! nucleated number from homogeneous freezing of so4
   real(kind_phys) :: niimm                     ! nucleated number from immersion freezing
   real(kind_phys) :: nidep                     ! nucleated number from deposition nucleation
   real(kind_phys) :: nimey                     ! nucleated number from deposition nucleation (meyers)
   real(kind_phys) :: n1, ni                    ! nucleated number
   real(kind_phys) :: tc, A, B                  ! work variable
   real(kind_phys) :: esl, esi, deles           ! work variable
   real(kind_phys) :: wbar1, wbar2

   ! used in SUBROUTINE Vpreice
   real(kind_phys) :: Ni_preice        ! cloud ice number conc (1/m3)
   real(kind_phys) :: lami,Ri_preice   ! mean cloud ice radius (m)
   real(kind_phys) :: Shom             ! initial ice saturation ratio; if <1, use hom threshold Si
   real(kind_phys) :: detaT,RHimean    ! temperature standard deviation, mean cloudy RHi
   real(kind_phys) :: wpicehet   ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at shet

   real(kind_phys) :: weffhet    ! effective Vertical velocity for ice nucleation (m/s)  weff=wbar-wpicehet

   logical  :: call_frm_zm, add_preexisting_ice
   !-------------------------------------------------------------------------------

   nuci = 0._kind_phys

   RHimean = relhum*svp_water(tair)/svp_ice(tair)*subgrid

   ! temp variables that depend on use_preexisting_ice
   wbar1 = wbar
   wbar2 = wbar

   ! If not using prexisting ice, the homogeneous freezing happens in the
   ! entire gridbox.
   fhom = 1._kind_phys

   if (present(call_frm_zm_in)) then
     call_frm_zm = call_frm_zm_in
   else
     call_frm_zm = .false.
   end if

   if (present(add_preexisting_ice_in)) then
     add_preexisting_ice = add_preexisting_ice_in
   else
     add_preexisting_ice = .true.
   end if

   if (use_preexisting_ice .and. (.not. call_frm_zm)) then

      Ni_preice = ni_in*rhoair                    ! (convert from #/kg -> #/m3)
      Ni_preice = Ni_preice / max(mincld,cldn)   ! in-cloud ice number density

      if (Ni_preice > 10.0_kind_phys .and. qi > 1.e-10_kind_phys) then    ! > 0.01/L = 10/m3
         Shom = -1.5_kind_phys   ! if Shom<1 , Shom will be recalculated in SUBROUTINE Vpreice, according to Ren & McKenzie, 2005, https://doi.org/10.1256/qj.04.126
         lami = (gamma4*ci*ni_in/qi)**(1._kind_phys/3._kind_phys)
         Ri_preice = 0.5_kind_phys/lami                  ! radius
         Ri_preice = max(Ri_preice, 1e-8_kind_phys)       ! >0.01micron
         call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shom, wpice)
         call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shet, wpicehet)
      else
         wpice    = 0.0_kind_phys
         wpicehet = 0.0_kind_phys
      end if

      weff     = max(wbar-wpice, minweff)
      wpice    = min(wpice, wbar)
      weffhet  = max(wbar-wpicehet,minweff)
      wpicehet = min(wpicehet, wbar)

      wbar1 = weff
      wbar2 = weffhet

      detaT   = wbar/0.23_kind_phys
      if (use_incloud_nuc) then
        call frachom(tair, 1._kind_phys, detaT, fhom)
      else
        call frachom(tair, RHimean, detaT, fhom)
      end if
   end if

   ni = 0._kind_phys
   tc = tair - 273.15_kind_phys

   ! initialize
   niimm = 0._kind_phys
   nidep = 0._kind_phys
   nihf  = 0._kind_phys
   deles = 0._kind_phys
   esi   = 0._kind_phys
   regm  = 0._kind_phys

   oso4_num  = 0._kind_phys
   odst_num  = 0._kind_phys
   osoot_num = 0._kind_phys

   if ((so4_num >= 1.0e-10_kind_phys .or. (soot_num+dst_num) >= 1.0e-10_kind_phys) .and. cldn > 0._kind_phys) then

      if (RHimean >= 1.2_kind_phys) then

         if ( ((tc <= 0.0_kind_phys).and.(tc >= -37.0_kind_phys).and.(qc < 1.e-12_kind_phys)).or.(tc <= -37.0_kind_phys)) then

            if ( (soot_num+dst_num) > 0._kind_phys)   then
               A = -1.4938_kind_phys * log(soot_num+dst_num) + 12.884_kind_phys
               B = -10.41_kind_phys  * log(soot_num+dst_num) - 67.69_kind_phys
               regm = A * log(wbar1) + B
            end if

            ! heterogeneous nucleation only
            if (tc > regm .or. so4_num < 1.0e-10_kind_phys) then

               if(tc < -40._kind_phys .and. wbar1 > 1._kind_phys .and. so4_num >= 1.0e-10_kind_phys) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                  call homogeneous_freezing(tc,wbar1,relhum*subgrid,so4_num,nihf)
                  niimm=0._kind_phys
                  nidep=0._kind_phys

                  ! If some homogeneous nucleation happened, assume all of the heterogeneous
                  ! and coarse mode sulfate particles nucleated.
                  if (nihf > 1e-3_kind_phys) then ! hom occur,  add preexisting ice
                     niimm     = dst_num + soot_num       ! assuming dst_num freeze firstly
                     odst_num  = dst_num
                     osoot_num = soot_num

                     oso4_num  = nihf
                  end if

                  nihf      = nihf * fhom
                  oso4_num  = oso4_num * fhom

                  n1        = nihf + niimm
               else

                  call heterogeneous_freezing(tc,wbar2,soot_num+dst_num,niimm,nidep)

                  nihf = 0._kind_phys
                  n1   = niimm + nidep

                  if ( (soot_num+dst_num) > 0._kind_phys)   then
                     osoot_num = soot_num * (niimm + nidep) / (soot_num + dst_num)
                     odst_num  = dst_num  * (niimm + nidep) / (soot_num + dst_num)
                  end if

               end if

            ! homogeneous nucleation only
            else if (tc < regm-5._kind_phys .or. (soot_num+dst_num) < 1.0e-10_kind_phys) then

               call homogeneous_freezing(tc,wbar1,relhum*subgrid,so4_num,nihf)
               niimm=0._kind_phys
               nidep=0._kind_phys

               ! If some homogeneous nucleation happened, assume all of the
               ! heterogeneous and coarse mode sulfate particles nucleated.
               if (nihf > 1e-3_kind_phys) then !  hom occur,  add preexisting ice
                  niimm     = dst_num + soot_num       ! assuming dst_num freeze firstly
                  odst_num  = dst_num
                  osoot_num = soot_num

                  oso4_num  = nihf
               end if

               nihf      = nihf * fhom
               oso4_num  = oso4_num * fhom

               n1        = nihf + niimm

            ! transition between homogeneous and heterogeneous: interpolate in-between
            else

               if (tc < -40._kind_phys .and. wbar1 > 1._kind_phys) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                  call homogeneous_freezing(tc, wbar1, relhum*subgrid, so4_num, nihf)
                  niimm = 0._kind_phys
                  nidep = 0._kind_phys

                  ! If some homogeneous nucleation happened, assume all of the
                  ! heterogeneous and coarse mode sulfate particles nucleated.
                  if (nihf > 1e-3_kind_phys) then ! hom occur,  add preexisting ice
                     niimm     = dst_num + soot_num       ! assuming dst_num freeze firstly
                     odst_num  = dst_num
                     osoot_num = soot_num

                     oso4_num  = nihf
                  end if

                  nihf      = nihf * fhom
                  oso4_num  = oso4_num * fhom

                  n1        = nihf + niimm

               else

                  call homogeneous_freezing(regm-5._kind_phys,wbar1,relhum*subgrid,so4_num,nihf)
                  call heterogeneous_freezing(regm,wbar2,soot_num+dst_num,niimm,nidep)

                  ! If some homogeneous nucleation happened, assume all of the
                  ! heterogeneous particles nucleated and add in a fraction of
                  ! the homogeneous freezing.
                  if (nihf > 1e-3_kind_phys) then ! hom occur,  add preexisting ice
                     oso4_num  = nihf
                  end if

                  if ( (soot_num+dst_num) > 0._kind_phys)   then
                     osoot_num = soot_num * (niimm + nidep) / (soot_num + dst_num)
                     odst_num  = dst_num  * (niimm + nidep) / (soot_num + dst_num)
                  end if

                  nihf      = nihf      * fhom * ((regm - tc) / 5._kind_phys)**2
                  oso4_num  = oso4_num  * fhom * ((regm - tc) / 5._kind_phys)**2

                  n1 = niimm + nidep + nihf

               end if
            end if

            ! Scale the rates for in-cloud number, since this is what
            ! MG is expecting to find.
            ni = n1

            ! If using prexsiting ice, and allowed to add, then add it to the total.
            if (use_preexisting_ice) then
               if (add_preexisting_ice .and. (.not. call_frm_zm)) then
                  ni = ni + Ni_preice * 1e-6_kind_phys
               end if
            end if
         end if
      end if
   end if

   ! deposition/condensation nucleation in mixed clouds (-37<T<0C) (Meyers, 1992)
   if(tc < 0._kind_phys .and. tc > -37._kind_phys .and. qc > 1.e-12_kind_phys) then
      esl = svp_water(tair)     ! over water in mixed clouds
      esi = svp_ice(tair)     ! over ice
      deles = (esl - esi)
      nimey=1.e-3_kind_phys*exp(12.96_kind_phys*deles/esi - 0.639_kind_phys)
   else
      nimey=0._kind_phys
   end if

   if (use_hetfrz_classnuc) nimey = 0._kind_phys

   nuci=ni + nimey

   if(nuci > 9999._kind_phys.or.nuci < 0._kind_phys) then
      write(iulog, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
      write(iulog, *) ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst_num,so4_num
      nuci=0._kind_phys
   end if

   nuci   = nuci*1.e+6_kind_phys/rhoair    ! change unit from #/cm3 to #/kg
   onimey = nimey*1.e+6_kind_phys/rhoair
   onidep = nidep*1.e+6_kind_phys/rhoair
   oniimm = niimm*1.e+6_kind_phys/rhoair
   onihf  = nihf*1.e+6_kind_phys/rhoair
end subroutine nucleati

!===============================================================================

pure subroutine heterogeneous_freezing(T,ww,Ns,Nis,Nid)

    real(kind_phys), intent(in)  :: T, ww, Ns
    real(kind_phys), intent(out) :: Nis, Nid

    ! Empirical fit coefficients for immersion freezing on dust, Liu & Penner (2005), p. 509
    real(kind_phys), parameter :: A11 = 0.0263_kind_phys
    real(kind_phys), parameter :: A12 = -0.0185_kind_phys
    real(kind_phys), parameter :: A21 = 2.758_kind_phys
    real(kind_phys), parameter :: A22 = 1.3221_kind_phys
    real(kind_phys), parameter :: B11 = -0.008_kind_phys
    real(kind_phys), parameter :: B12 = -0.0468_kind_phys
    real(kind_phys), parameter :: B21 = -0.2667_kind_phys
    real(kind_phys), parameter :: B22 = -1.4588_kind_phys

    real(kind_phys) :: B,C

!---------------------------------------------------------------------

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0_kind_phys    ! don't include deposition nucleation for cirrus clouds when T<-37C

end subroutine heterogeneous_freezing

!===============================================================================

pure subroutine homogeneous_freezing(T,ww,RH,Na,Ni)

      real(kind_phys), intent(in)  :: T, ww, RH, Na
      real(kind_phys), intent(out) :: Ni

      ! Empirical fit coefficients for homogeneous freezing, Liu & Penner (2005), p. 504
      real(kind_phys), parameter :: A1_fast  = 0.0231_kind_phys
      real(kind_phys), parameter :: A21_fast = -1.6387_kind_phys  !(T>-64 deg)
      real(kind_phys), parameter :: A22_fast = -6.045_kind_phys   !(T<=-64 deg)
      real(kind_phys), parameter :: B1_fast  = -0.008_kind_phys
      real(kind_phys), parameter :: B21_fast = -0.042_kind_phys   !(T>-64 deg)
      real(kind_phys), parameter :: B22_fast = -0.112_kind_phys   !(T<=-64 deg)
      real(kind_phys), parameter :: C1_fast  = 0.0739_kind_phys
      real(kind_phys), parameter :: C2_fast  = 1.2372_kind_phys
      real(kind_phys), parameter :: A1_slow  = -0.3949_kind_phys
      real(kind_phys), parameter :: A2_slow  = 1.282_kind_phys
      real(kind_phys), parameter :: B1_slow  = -0.0156_kind_phys
      real(kind_phys), parameter :: B2_slow  = 0.0111_kind_phys
      real(kind_phys), parameter :: B3_slow  = 0.0217_kind_phys
      real(kind_phys), parameter :: C1_slow  = 0.120_kind_phys
      real(kind_phys), parameter :: C2_slow  = 2.312_kind_phys

      real(kind_phys) :: A2_fast,B2_fast
      real(kind_phys) :: k1_fast,k2_fast
      real(kind_phys) :: k1_slow,k2_slow
      real(kind_phys) :: regm
      real(kind_phys) :: A,B,C
      real(kind_phys) :: RHw

!---------------------------------------------------------------------

      Ni = 0.0_kind_phys

!----------------------------
!RHw parameters
      A = 6.0e-4_kind_phys*log(ww)+6.6e-3_kind_phys
      B = 6.0e-2_kind_phys*log(ww)+1.052_kind_phys
      C = 1.68_kind_phys  *log(ww)+129.35_kind_phys
      RHw=(A*T*T+B*T+C)*0.01_kind_phys

      if((T<=-37.0_kind_phys) .and. ((RH) >= RHw)) then

        regm = 6.07_kind_phys*log(ww)-55.0_kind_phys

        if(T >= regm) then    ! fast-growth regime

          if(T > -64.0_kind_phys) then
            A2_fast=A21_fast
            B2_fast=B21_fast
          else
            A2_fast=A22_fast
            B2_fast=B22_fast
          end if

          k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
          k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

          Ni = k1_fast*Na**(k2_fast)
          Ni = min(Ni,Na)

        else       ! slow-growth regime

          k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
          k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

          Ni = k1_slow*Na**(k2_slow)
          Ni = min(Ni,Na)

        end if

      end if

end subroutine homogeneous_freezing

!===============================================================================

pure subroutine Vpreice(P_in, T_in, R_in, C_in, S_in, V_out)

   !  based on  Karcher et al. (2006), https://doi.org/10.1029/2005JD006219
   !  VERTICAL VELOCITY CALCULATED FROM DEPOSITIONAL LOSS TERM

   REAL(kind_phys), INTENT(in)  :: P_in       ! [Pa],INITIAL AIR pressure
   REAL(kind_phys), INTENT(in)  :: T_in       ! [K] ,INITIAL AIR temperature
   REAL(kind_phys), INTENT(in)  :: R_in       ! [m],INITIAL MEAN  ICE CRYSTAL NUMBER RADIUS
   REAL(kind_phys), INTENT(in)  :: C_in       ! [m-3],INITIAL TOTAL ICE CRYSTAL NUMBER DENSITY, [1/cm3]
   REAL(kind_phys), INTENT(in)  :: S_in       ! [-],INITIAL ICE SATURATION RATIO;; if <1, use hom threshold Si
   REAL(kind_phys), INTENT(out) :: V_out      ! [m/s], VERTICAL VELOCITY REDUCTION (caused by preexisting ice)

   ! SUBROUTINE parameters
   REAL(kind_phys), PARAMETER :: ALPHAc  = 0.5_kind_phys ! density of ice (g/cm3), !!!V is not related to ALPHAc
   REAL(kind_phys), PARAMETER :: FA1c    = 0.601272523_kind_phys
   REAL(kind_phys), PARAMETER :: FA2c    = 0.000342181855_kind_phys
   REAL(kind_phys), PARAMETER :: FA3c    = 1.49236645E-12_kind_phys
   REAL(kind_phys), PARAMETER :: WVP1c   = 3.6E+10_kind_phys
   REAL(kind_phys), PARAMETER :: WVP2c   = 6145.0_kind_phys
   REAL(kind_phys), PARAMETER :: FVTHc   = 11713803.0_kind_phys
   REAL(kind_phys), PARAMETER :: THOUBKc = 7.24637701E+18_kind_phys
   REAL(kind_phys), PARAMETER :: SVOLc   = 3.23E-23_kind_phys    ! SVOL=XMW/RHOICE
   REAL(kind_phys), PARAMETER :: FDc     = 249.239822_kind_phys
   REAL(kind_phys), PARAMETER :: FPIVOLc = 3.89051704E+23_kind_phys
   REAL(kind_phys) :: T,P,S,R,C
   REAL(kind_phys) :: A1,A2,A3,B1,B2
   REAL(kind_phys) :: T_1,PICE,FLUX,ALP4,CISAT,DLOSS,VICE

   T = T_in          ! K  , K
   P = P_in*1e-2_kind_phys  ! Pa , hpa

   if (s_in < 1.0_kind_phys) then
      S = 2.349_kind_phys - (T/259.0_kind_phys) ! homogeneous freezing threshold, according to Ren & McKenzie, 2005
   else
      S = S_in                    ! INPUT ICE SATURATION RATIO, -,  >1
   end if

   R     = R_in*1e2_kind_phys   ! m  => cm
   C     = C_in*1e-6_kind_phys  ! m-3 => cm-3
   T_1   = 1.0_kind_phys/ T
   PICE  = WVP1c * EXP(-(WVP2c*T_1))
   ALP4  = 0.25_kind_phys * ALPHAc
   FLUX  = ALP4 * SQRT(FVTHc*T)
   CISAT = THOUBKc * PICE * T_1
   A1    = ( FA1c * T_1 - FA2c ) * T_1
   A2    = 1.0_kind_phys/ CISAT
   A3    = FA3c * T_1 / P
   B1    = FLUX * SVOLc * CISAT * ( S-1.0_kind_phys )
   B2    = FLUX * FDc * P * T_1**1.94_kind_phys
   DLOSS = FPIVOLc * C * B1 * R**2 / ( 1.0_kind_phys+ B2 * R )
   VICE  = ( A2 + A3 * S ) * DLOSS / ( A1 * S )  ! 2006,(19)
   V_out = VICE*1e-2_kind_phys  ! cm/s => m/s

end subroutine Vpreice

pure subroutine frachom(Tmean,RHimean,detaT,fhom)
   ! How much fraction of cirrus might reach Shom
   ! base on "A cirrus cloud scheme for general circulation models",
   ! B. Karcher and U. Burkhardt 2008, https://doi.org/10.1002/qj.301

   real(kind_phys), intent(in)  :: Tmean, RHimean, detaT
   real(kind_phys), intent(out) :: fhom

   real(kind_phys), parameter :: seta = 6132.9_kind_phys  ! K
   integer,  parameter :: Nbin=200          ! (Tmean - 3*detaT, Tmean + 3*detaT)

   real(kind_phys) :: PDF_T(Nbin)    ! temperature PDF;  ! PDF_T=0  outside (Tmean-3*detaT, Tmean+3*detaT)
   real(kind_phys) :: Sbin(Nbin)     ! the fluctuations of Si that are driven by the T variations
   real(kind_phys) :: Sihom, deta
   integer  :: i

   Sihom = 2.349_kind_phys-Tmean/259.0_kind_phys   ! homogeneous freezing threshold, according to Ren & McKenzie, 2005
   fhom  = 0.0_kind_phys

   do i = Nbin, 1, -1

      deta     = (i - 0.5_kind_phys - Nbin/2)*6.0_kind_phys/Nbin   ! PDF_T=0  outside (Tmean-3*detaT, Tmean+3*detaT)
      Sbin(i)  = RHimean*exp(deta*detaT*seta/Tmean**2.0_kind_phys)
      PDF_T(i) = exp(-deta**2.0_kind_phys/2.0_kind_phys)*6.0_kind_phys/(sqrt(2.0_kind_phys*Pi)*Nbin)


      if (Sbin(i) >= Sihom) then
         fhom = fhom + PDF_T(i)
      else
         exit
      end if
   end do

   fhom = min(1.0_kind_phys, fhom/0.997_kind_phys)   ! accounting for the finite limits (-3 , 3)
end subroutine frachom

end module nucleate_ice
