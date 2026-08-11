! modal_aero_newnuc.F90

!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_newnuc --- modal aerosol new-particle nucleation
!
! !INTERFACE:
module modal_aero_newnuc

! !USES:
  use ccpp_kinds, only: kind_phys

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:
  public modal_aero_newnuc_run, modal_aero_newnuc_init

! !PUBLIC DATA MEMBERS:
! species indices in the host constituent space, set by modal_aero_newnuc_init;
! default 0 = bypass (host found no h2so4 or no aitken so4/num species)
  integer  :: l_h2so4_sv = 0, l_nh3_sv = 0, lnumait_sv = 0, lnh4ait_sv = 0, lso4ait_sv = 0

! min h2so4 vapor for nuc calcs = 4.0e-16 mol/mol-air ~= 1.0e4 molecules/cm3,
  real(kind_phys), parameter :: qh2so4_cutoff = 4.0e-16_kind_phys

  real(kind_phys) :: dens_so4a_host
  real(kind_phys) :: mw_nh4a_host, mw_so4a_host

! aitken-mode geometry (nominal size + dry-diameter limits) from host mode metadata
  real(kind_phys) :: dgnum_aitken, dgnumhi_aitken, dgnumlo_aitken

! host physical constants, set by modal_aero_newnuc_init (passed from the host
! rather than hardcoded so the values stay bit-identical with the host's)
  real(kind_phys) :: pi
  real(kind_phys) :: rgas              ! gas constant (J/K/mol)
  real(kind_phys) :: avogad            ! Avogadro's number (1/mol)
  real(kind_phys) :: mw_so4a           ! molecular weight of sulfate (g/mol)
  real(kind_phys) :: mw_nh4a           ! molecular weight of ammonium (g/mol)
  real(kind_phys) :: r_universal       ! universal gas constant (J/K/kmol)

! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   R.Easter 2007.09.14:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_newnuc_run --- ...
!
! !INTERFACE:
  subroutine modal_aero_newnuc_run( &
    ncol, pver, top_lev, &
    num_q, loffset, deltat, &
    t, pmid, pdel, &
    zm, pblh, &
    qv, cld, &
    q, gravit, &
    del_h2so4_gasprod, del_h2so4_aeruptk, &
    dqdt, dotend, qsrflx, &
    errmsg, errflg)

! !USES:
    use wv_saturation, only: qsat

! !PARAMETERS:
    integer, intent(in)  :: ncol             ! number of columns in chunk
    integer, intent(in)  :: pver             ! number of vertical levels
    integer, intent(in)  :: top_lev          ! top level for modal aerosol calculations
    integer, intent(in)  :: num_q            ! number of species in q/dqdt (= gas_pcnst)
    integer, intent(in)  :: loffset          ! offset applied to modal aero "pointers"
    real(kind_phys), intent(in) :: deltat           ! model timestep (s)

    real(kind_phys), intent(in) :: t(:, :)           ! (ncol,pver) temperature (K)
    real(kind_phys), intent(in) :: pmid(:, :)        ! (ncol,pver) pressure at model levels (Pa)
    real(kind_phys), intent(in) :: pdel(:, :)        ! (ncol,pver) pressure thickness of levels (Pa)
    real(kind_phys), intent(in) :: zm(:, :)          ! (ncol,pver) midpoint height above surface (m)
    real(kind_phys), intent(in) :: pblh(:)          ! (ncol) pbl height (m)
    real(kind_phys), intent(in) :: qv(:, :)          ! (ncol,pver) specific humidity (kg/kg)
    real(kind_phys), intent(in) :: cld(:, :)         ! (ncol,pver) stratiform cloud fraction
    real(kind_phys), intent(in) :: q(:, :, :)         ! (ncol,pver,num_q)
    ! tracer mixing ratio (TMR) array
    ! *** MUST BE mol/mol-air or #/mol-air
    real(kind_phys), intent(in) :: gravit           ! gravitational acceleration (m/s2)
    real(kind_phys), intent(in) :: del_h2so4_gasprod(:, :)
    ! (ncol,pver) h2so4 gas-phase production
    ! change over deltat (mol/mol)
    real(kind_phys), intent(in) :: del_h2so4_aeruptk(:, :)
    ! (ncol,pver) h2so4 gas-phase loss to
    ! aerosol over deltat (mol/mol)
    real(kind_phys), intent(out) :: dqdt(:, :, :)     ! (ncol,pver,num_q) TMR tendency array
    logical, intent(out) :: dotend(:)       ! (num_q) flag for doing tendency
    real(kind_phys), intent(out) :: qsrflx(:, :, :)   ! (ncol,num_q,1)
    ! process-specific column tracer tendencies
    ! 1 = nucleation (for aerocom)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

! !DESCRIPTION:
!   computes changes due to aerosol nucleation (new particle formation)
!       treats both nucleation and subsequent growth of new particles
!           to aitken mode size
!   uses the following parameterizations
!       vehkamaki et al. (2002) parameterization for binary
!           homogeneous nucleation (h2so4-h2o) plus
!       kerminen and kulmala (2002) parameterization for
!           new particle loss during growth to aitken size
!
! !REVISION HISTORY:
!   R.Easter 2007.09.14:  Adapted from MIRAGE2 code and CMAQ V4.6 code
!
!EOP
!----------------------------------------------------------------------
!BOC

!   local variables
    integer :: i, itmp, k, lun
    integer :: lnumait, lso4ait, lnh4ait
    integer :: l_h2so4, l_nh3
    integer :: ldiagveh02
    integer, parameter :: ldiag1 = -1, ldiag2 = -1, ldiag3 = -1, ldiag4 = -1
    integer, parameter :: newnuc_method_flagaa = 11
!       integer, parameter :: newnuc_method_flagaa = 12
    !  1=merikanto et al (2007) ternary   2=vehkamaki et al (2002) binary
    ! 11=merikanto ternary + first-order boundary layer
    ! 12=merikanto ternary + second-order boundary layer

    real(kind_phys) :: adjust_factor
    real(kind_phys) :: aircon
    real(kind_phys) :: cldx
    real(kind_phys) :: dens_nh4so4a
    real(kind_phys) :: dmdt_ait, dmdt_aitsv1, dmdt_aitsv2, dmdt_aitsv3
    real(kind_phys) :: dndt_ait, dndt_aitsv1, dndt_aitsv2, dndt_aitsv3
    real(kind_phys) :: dndt(ncol, pver) ! nucleation rate (#/m3/s)
    real(kind_phys) :: dnh4dt_ait, dso4dt_ait
    real(kind_phys) :: dpnuc
    real(kind_phys) :: dplom_mode(1), dphim_mode(1)
    real(kind_phys) :: ev_sat(ncol, pver)
    real(kind_phys) :: mass1p
    real(kind_phys) :: mass1p_aithi, mass1p_aitlo
    real(kind_phys) :: pdel_fac
    real(kind_phys) :: qh2so4_cur, qh2so4_avg, qh2so4_del
    real(kind_phys) :: qnh3_cur, qnh3_del, qnh4a_del
    real(kind_phys) :: qnuma_del
    real(kind_phys) :: qso4a_del
    real(kind_phys) :: qv_sat(ncol, pver)
    real(kind_phys) :: qvswtr
    real(kind_phys) :: relhum, relhumav, relhumnn
    real(kind_phys) :: tmpa, tmpb, tmpc
    real(kind_phys) :: tmp_q1, tmp_q2, tmp_q3
    real(kind_phys) :: tmp_frso4, tmp_uptkrate

    logical  :: do_nh3                   ! flag for doing nh3/nh4

    character(len=1) :: tmpch1, tmpch2, tmpch3

! begin
    errmsg = ' '
    errflg = 0

    lun = 6

!--------------------------------------------------------------------------------
!!$   if (ldiag1 > 0) then
!!$   do i = 1, ncol
!!$   if (lonndx(i) /= 37) cycle
!!$   if (latndx(i) /= 23) cycle
!!$   if (nstep > 3)       cycle
!!$   write( lun, '(/a,i7,3i5,f10.2)' )   &
!!$         '*** modal_aero_newnuc_sub -- nstep, iam, lat, lon =',   &
!!$         nstep, iam, latndx(i), lonndx(i)
!!$   end do
!!$   if (nstep > 3) call endrun( '*** modal_aero_newnuc_sub -- testing halt after step 3' )
!!$!  if (ncol /= -999888777) return
!!$   end if
!--------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!   zero the tendency outputs up front: they are intent(out) and the caller
!   applies/outputs them unconditionally, including on the bypass path below
    dotend(:) = .false.
    dqdt(1:ncol, :, :) = 0.0_kind_phys
    qsrflx(1:ncol, :, :) = 0.0_kind_phys
    dndt(1:ncol, :) = 0.0_kind_phys

    l_h2so4 = l_h2so4_sv - loffset
    l_nh3 = l_nh3_sv - loffset
    lnumait = lnumait_sv - loffset
    lnh4ait = lnh4ait_sv - loffset
    lso4ait = lso4ait_sv - loffset

!   skip if no aitken mode OR if no h2so4 species
    if ((l_h2so4 <= 0) .or. (lso4ait <= 0) .or. (lnumait <= 0)) return

!   set dotend
    dotend(lnumait) = .true.
    dotend(lso4ait) = .true.
    dotend(l_h2so4) = .true.

    if ((l_nh3 > 0) .and. (l_nh3 <= num_q) .and. &
        (lnh4ait > 0) .and. (lnh4ait <= num_q)) then
      do_nh3 = .true.
      dotend(lnh4ait) = .true.
      dotend(l_nh3) = .true.
    else
      do_nh3 = .false.
    end if

!   dry-diameter limits for "grown" new particles
    dplom_mode(1) = exp(0.67_kind_phys*log(dgnumlo_aitken) &
                        + 0.33_kind_phys*log(dgnum_aitken))
    dphim_mode(1) = dgnumhi_aitken

!   mass1p_... = mass (kg) of so4 & nh4 in a single particle of diameter ...
!                (assuming same dry density for so4 & nh4)
!       mass1p_aitlo - dp = dplom_mode(1)
!       mass1p_aithi - dp = dphim_mode(1)
    tmpa = dens_so4a_host*pi/6.0_kind_phys
    mass1p_aitlo = tmpa*(dplom_mode(1)**3)
    mass1p_aithi = tmpa*(dphim_mode(1)**3)

!   compute qv_sat = saturation specific humidity
    do k = 1, pver
      call qsat(t(1:ncol, k), pmid(1:ncol, k), ev_sat(1:ncol, k), qv_sat(1:ncol, k), ncol)
    end do
!
!   loop over levels and columns to calc the renaming
!
    main_k: do k = top_lev, pver
    main_i: do i = 1, ncol

!   skip if completely cloudy,
!   because all h2so4 vapor should be cloud-borne
      if (cld(i, k) >= 0.99_kind_phys) cycle main_i

!   qh2so4_cur = current qh2so4, after aeruptk
      qh2so4_cur = q(i, k, l_h2so4)
!   skip if h2so4 vapor < qh2so4_cutoff
      if (qh2so4_cur <= qh2so4_cutoff) cycle main_i

      tmpa = max(0.0_kind_phys, del_h2so4_gasprod(i, k))
      tmp_q3 = qh2so4_cur
!   tmp_q2 = qh2so4 before aeruptk
!   (note tmp_q3, tmp_q2 both >= 0.0)
      tmp_q2 = tmp_q3 + max(0.0_kind_phys, -del_h2so4_aeruptk(i, k))

!   *** temporary -- in order to get more nucleation
!       qh2so4_cur = qh2so4_cur*1.0e1
!       tmp_q3 = tmp_q3*1.0e1
!       tmp_q2 = tmp_q2*1.0e1
!       tmpa   = tmpa  *1.0e1

!   tmpb = log( tmp_q2/tmp_q3 ) BUT with some checks added
!       tmp_uptkrate = tmpb/deltat
      if (tmp_q2 <= tmp_q3) then
        tmpb = 0.0_kind_phys
      else
        tmpc = tmp_q2*exp(-20.0_kind_phys)
        if (tmp_q3 <= tmpc) then
          tmp_q3 = tmpc
          tmpb = 20.0_kind_phys
        else
          tmpb = log(tmp_q2/tmp_q3)
        end if
      end if
!   d[ln(qh2so4)]/dt (1/s) from uptake (condensation) to aerosol
      tmp_uptkrate = tmpb/deltat

!   qh2so4_avg = estimated average qh2so4
!   when production & loss are done simultaneously
      if (tmpb <= 0.1_kind_phys) then
        qh2so4_avg = tmp_q3*(1.0_kind_phys + 0.5_kind_phys*tmpb) - 0.5_kind_phys*tmpa
      else
        tmpc = tmpa/tmpb
        qh2so4_avg = (tmp_q3 - tmpc)*((exp(tmpb) - 1.0_kind_phys)/tmpb) + tmpc
      end if
      if (qh2so4_avg <= qh2so4_cutoff) cycle main_i

      if (do_nh3) then
        qnh3_cur = max(0.0_kind_phys, q(i, k, l_nh3))
      else
        qnh3_cur = 0.0_kind_phys
      end if

!   relhumav = grid average RH
      qvswtr = qv_sat(i, k)
      qvswtr = max(qvswtr, 1.0e-20_kind_phys)
      relhumav = qv(i, k)/qvswtr
      relhumav = max(0.0_kind_phys, min(1.0_kind_phys, relhumav))
!   relhum = non-cloudy area RH
      cldx = max(0.0_kind_phys, cld(i, k))
      relhum = (relhumav - cldx)/(1.0_kind_phys - cldx)
      relhum = max(0.0_kind_phys, min(1.0_kind_phys, relhum))
!   limit RH to between 0.1% and 99%
      relhumnn = relhum
      relhumnn = max(0.01_kind_phys, min(0.99_kind_phys, relhumnn))

!   aircon = air concentration (mol-air/m3)
      aircon = 1.0e3_kind_phys*pmid(i, k)/(r_universal*t(i, k))

!   call ... routine to get nucleation rates
      ldiagveh02 = -1
!!$     if (ldiag2 > 0) then
!!$     if ((lonndx(i) == 37) .and. (latndx(i) == 23)) then
!!$     if ((k >= 24) .or. (mod(k,4) == 0)) then
!!$         ldiagveh02 = +1
!!$            write(lun,'(/a,i8,3i4,f8.2,1p,4e10.2)')   &
!!$             'veh02 call - nstep,lat,lon,k; tk,rh,p,cair',   &
!!$             nstep, latndx(i), lonndx(i), k,   &
!!$             t(i,k), relhumnn, pmid(k,k), aircon
!!$     end if
!!$     end if
!!$     end if
      call mer07_veh02_nuc_mosaic_1box( &
        newnuc_method_flagaa, &
        deltat, t(i, k), relhumnn, pmid(i, k), &
        zm(i, k), pblh(i), &
        qh2so4_cur, qh2so4_avg, qnh3_cur, tmp_uptkrate, &
        mw_so4a_host, &
        1, 1, dplom_mode, dphim_mode, &
        itmp, qnuma_del, qso4a_del, qnh4a_del, &
        qh2so4_del, qnh3_del, dens_nh4so4a, ldiagveh02)
!          qh2so4_del, qnh3_del, dens_nh4so4a )
!----------------------------------------------------------------------
!       subr mer07_veh02_nuc_mosaic_1box(   &
!          newnuc_method_flagaa,   &
!          dtnuc, temp_in, rh_in, press_in,   &
!          qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
!          nsize, maxd_asize, dplom_sect, dphim_sect,   &
!          isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
!          qh2so4_del, qnh3_del, dens_nh4so4a )
!
!! subr arguments (in)
!        real(kind_phys), intent(in) :: dtnuc             ! nucleation time step (s)
!        real(kind_phys), intent(in) :: temp_in           ! temperature, in k
!        real(kind_phys), intent(in) :: rh_in             ! relative humidity, as fraction
!        real(kind_phys), intent(in) :: press_in          ! air pressure (pa)
!
!        real(kind_phys), intent(in) :: qh2so4_cur, qh2so4_avg
!                                                  ! gas h2so4 mixing ratios (mol/mol-air)
!        real(kind_phys), intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
!             ! qxxx_cur = current value (after gas chem and condensation)
!             ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
!        real(kind_phys), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)

!
!        integer, intent(in) :: nsize                    ! number of aerosol size bins
!        integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
!        real(kind_phys), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
!        real(kind_phys), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)
!
!! subr arguments (out)
!        integer, intent(out) :: isize_nuc         ! size bin into which new particles go
!        real(kind_phys), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
!        real(kind_phys), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
!        real(kind_phys), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
!        real(kind_phys), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
!        real(kind_phys), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
!                                                  ! aerosol changes are > 0; gas changes are < 0
!        real(kind_phys), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)
!----------------------------------------------------------------------

!   convert qnuma_del from (#/mol-air) to (#/kmol-air)
      qnuma_del = qnuma_del*1.0e3_kind_phys
!   number nuc rate (#/kmol-air/s) from number nuc amt
      dndt_ait = qnuma_del/deltat
!   fraction of mass nuc going to so4
      tmpa = qso4a_del*mw_so4a_host
      tmpb = tmpa + qnh4a_del*mw_nh4a_host
      tmp_frso4 = max(tmpa, 1.0e-35_kind_phys)/max(tmpb, 1.0e-35_kind_phys)
!   mass nuc rate (kg/kmol-air/s or g/mol...) hhfrom mass nuc amts
      dmdt_ait = max(0.0_kind_phys, (tmpb/deltat))

      dndt_aitsv1 = dndt_ait
      dmdt_aitsv1 = dmdt_ait
      dndt_aitsv2 = 0.0_kind_phys
      dmdt_aitsv2 = 0.0_kind_phys
      dndt_aitsv3 = 0.0_kind_phys
      dmdt_aitsv3 = 0.0_kind_phys
      tmpch1 = ' '
      tmpch2 = ' '

      if (dndt_ait < 1.0e2_kind_phys) then
!   ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
        dndt_ait = 0.0_kind_phys
        dmdt_ait = 0.0_kind_phys
        tmpch1 = 'A'

      else
        dndt_aitsv2 = dndt_ait
        dmdt_aitsv2 = dmdt_ait
        tmpch1 = 'B'

!   mirage2 code checked for complete h2so4 depletion here,
!   but this is now done in mer07_veh02_nuc_mosaic_1box
        mass1p = dmdt_ait/dndt_ait
        dndt_aitsv3 = dndt_ait
        dmdt_aitsv3 = dmdt_ait

!   apply particle size constraints
        if (mass1p < mass1p_aitlo) then
!   reduce dndt to increase new particle size
          dndt_ait = dmdt_ait/mass1p_aitlo
          tmpch1 = 'C'
        else if (mass1p > mass1p_aithi) then
!   reduce dmdt to decrease new particle size
          dmdt_ait = dndt_ait*mass1p_aithi
          tmpch1 = 'E'
        end if
      end if

! *** apply adjustment factor to avoid unrealistically high
!     aitken number concentrations in mid and upper troposphere
!       adjust_factor = 0.5
!       dndt_ait = dndt_ait * adjust_factor
!       dmdt_ait = dmdt_ait * adjust_factor

!   set tendencies
      pdel_fac = pdel(i, k)/gravit

!   dso4dt_ait, dnh4dt_ait are (kmol/kmol-air/s)
      dso4dt_ait = dmdt_ait*tmp_frso4/mw_so4a_host
      dnh4dt_ait = dmdt_ait*(1.0_kind_phys - tmp_frso4)/mw_nh4a_host

!   dqdt tendencies are returned to the caller, which applies q = q + dqdt*deltat
      dqdt(i, k, l_h2so4) = -dso4dt_ait*(1.0_kind_phys - cldx)
      qsrflx(i, l_h2so4, 1) = qsrflx(i, l_h2so4, 1) + dqdt(i, k, l_h2so4)*pdel_fac

      dqdt(i, k, lso4ait) = dso4dt_ait*(1.0_kind_phys - cldx)
      qsrflx(i, lso4ait, 1) = qsrflx(i, lso4ait, 1) + dqdt(i, k, lso4ait)*pdel_fac
      if (lnumait > 0) then
        dqdt(i, k, lnumait) = dndt_ait*(1.0_kind_phys - cldx)
!   dndt is (#/m3/s), dqdt(:,:,lnumait) is (#/kmol-air/s), aircon is (mol-air/m3)
        dndt(i, k) = dqdt(i, k, lnumait)*aircon*1.0e-3_kind_phys
        qsrflx(i, lnumait, 1) = qsrflx(i, lnumait, 1) &
                                + dqdt(i, k, lnumait)*pdel_fac
      end if

      if ((do_nh3) .and. (dnh4dt_ait > 0.0_kind_phys)) then
        dqdt(i, k, l_nh3) = -dnh4dt_ait*(1.0_kind_phys - cldx)
        qsrflx(i, l_nh3, 1) = qsrflx(i, l_nh3, 1) + dqdt(i, k, l_nh3)*pdel_fac

        dqdt(i, k, lnh4ait) = dnh4dt_ait*(1.0_kind_phys - cldx)
        qsrflx(i, lnh4ait, 1) = qsrflx(i, lnh4ait, 1) + dqdt(i, k, lnh4ait)*pdel_fac
      end if

!!   temporary diagnostic
!        if (ldiag3 > 0) then
!        if ((dndt_ait /= 0.0_kind_phys) .or. (dmdt_ait /= 0.0_kind_phys)) then
!           write(lun,'(3a,1x,i7,3i5,1p,5e12.4)')   &
!              'newnucxx', tmpch1, tmpch2, nstep, lchnk, i, k,   &
!              dndt_ait, dmdt_ait, cldx
!!          call endrun( 'modal_aero_newnuc_sub' )
!        end if
!        end if

!   diagnostic output start ----------------------------------------
!!$     if (ldiag4 > 0) then
!!$     if ((lonndx(i) == 37) .and. (latndx(i) == 23)) then
!!$     if ((k >= 24) .or. (mod(k,4) == 0)) then
!!$        write(lun,97010) nstep, latndx(i), lonndx(i), k, t(i,k), aircon
!!$        write(lun,97020) 'pmid, pdel                   ',   &
!!$                pmid(i,k), pdel(i,k)
!!$        write(lun,97030) 'qv,qvsw, cld, rh_av, rh_clr  ',   &
!!$                qv(i,k), qvswtr, cldx, relhumav, relhum
!!$        write(lun,97020) 'h2so4_cur, _pre, _av, nh3_cur',   &
!!$             qh2so4_cur, tmp_q2, qh2so4_avg, qnh3_cur
!!$        write(lun,97020) 'del_h2so4_gasprod, _aeruptk  ',   &
!!$             del_h2so4_gasprod(i,k), del_h2so4_aeruptk(i,k),   &
!!$             tmp_uptkrate*3600.0_kind_phys
!!$        write(lun,97020) ' '
!!$        write(lun,97050) 'tmpch1, tmpch2               ', tmpch1, tmpch2
!!$        write(lun,97020) 'dndt_, dmdt_aitsv1           ',   &
!!$                              dndt_aitsv1, dmdt_aitsv1
!!$        write(lun,97020) 'dndt_, dmdt_aitsv2           ',   &
!!$                              dndt_aitsv2, dmdt_aitsv2
!!$        write(lun,97020) 'dndt_, dmdt_aitsv3           ',   &
!!$                              dndt_aitsv3, dmdt_aitsv3
!!$        write(lun,97020) 'dndt_, dmdt_ait              ',   &
!!$                              dndt_ait, dmdt_ait
!!$        write(lun,97020) 'dso4dt_, dnh4dt_ait          ',   &
!!$                              dso4dt_ait, dnh4dt_ait
!!$        write(lun,97020) 'qso4a_del, qh2so4_del        ',   &
!!$                              qso4a_del, qh2so4_del
!!$        write(lun,97020) 'qnh4a_del, qnh3_del          ',   &
!!$                              qnh4a_del, qnh3_del
!!$        write(lun,97020) 'dqdt(h2so4), (nh3)           ',   &
!!$              dqdt(i,k,l_h2so4), dqdt(i,k,l_nh3)
!!$        write(lun,97020) 'dqdt(so4a), (nh4a), (numa)   ',   &
!!$              dqdt(i,k,lso4ait), dqdt(i,k,lnh4ait), dqdt(i,k,lnumait)
!!$
!!$     dpnuc = 0.0_kind_phys
!!$     if (dndt_aitsv1 > 1.0e-5_kind_phys) dpnuc = (6.0_kind_phys*dmdt_aitsv1/   &
!!$                     (pi*dens_so4a_host*dndt_aitsv1))**0.3333333_kind_phys
!!$        if (dpnuc > 0.0_kind_phys) then
!!$        write(lun,97020) 'dpnuc,      dp_aitlo, _aithi ',   &
!!$                      dpnuc, dplom_mode(1), dphim_mode(1)
!!$        write(lun,97020) 'mass1p, mass1p_aitlo, _aithi ',   &
!!$                      mass1p, mass1p_aitlo, mass1p_aithi
!!$        end if
!!$
!!$ 97010  format( / 'NEWNUC nstep,lat,lon,k,tk,cair', i8, 3i4, f8.2, 1pe12.4 )
!!$ 97020  format( a, 1p, 6e12.4 )
!!$ 97030  format( a, 1p, 2e12.4, 0p, 5f10.6 )
!!$ 97040  format( 29x, 1p, 6e12.4 )
!!$ 97050  format( a, 2(3x,a) )
!!$        end if
!!$        end if
!!$        end if
!   diagnostic output end   ------------------------------------------

    end do main_i
    end do main_k

!   history file column-tendency fields (adv_mass/mwdry scaling + outfld of
!   qsrflx) are done by the caller, which owns the host constituent metadata

!EOC
  end subroutine modal_aero_newnuc_run

!----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine mer07_veh02_nuc_mosaic_1box( &
    newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in, &
    zm_in, pblh_in, &
    qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate, &
    mw_so4a_host, &
    nsize, maxd_asize, dplom_sect, dphim_sect, &
    isize_nuc, qnuma_del, qso4a_del, qnh4a_del, &
    qh2so4_del, qnh3_del, dens_nh4so4a, ldiagaa)
!          qh2so4_del, qnh3_del, dens_nh4so4a )
!   (rgas, avogad, mw_so4a, mw_nh4a are module-level host constants
!    set by modal_aero_newnuc_init)
!.......................................................................
!
! calculates new particle production from homogeneous nucleation
!    over timestep dtnuc, using nucleation rates from either
!    merikanto et al. (2007) h2so4-nh3-h2o ternary parameterization
!    vehkamaki et al. (2002) h2so4-h2o binary parameterization
!
! the new particles are "grown" to the lower-bound size of the host code's
!    smallest size bin.  (this "growth" is somewhat ad hoc, and would not be
!    necessary if the host code's size bins extend down to ~1 nm.)
!
!    if the h2so4 and nh3 mass mixing ratios (mixrats) of the grown new
!    particles exceed the current gas mixrats, the new particle production
!    is reduced so that the new particle mass mixrats match the gas mixrats.
!
!    the correction of kerminen and kulmala (2002) is applied to account
!    for loss of the new particles by coagulation as they are
!    growing to the "host code mininum size"
!
! revision history
!    coded by rc easter, pnnl, xx-apr-2007
!
! key routines called: subr ternary_nuc_napari
!
! references:
!    merikanto, j., i. napari, h. vehkamaki, t. anttila,
!     and m. kulmala, 2007, new parameterization of
!     sulfuric acid-ammonia-water ternary nucleation
!     rates at tropospheric conditions,
!       j. geophys. res., 112, d15207, doi:10.1029/2006jd0027977
!
!    vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!       c. timmreck, m. noppel and a. laaksonen, 2002,
!       an improved parameterization for sulfuric acid-water nucleation
!       rates for tropospheric and stratospheric conditions,
!       j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
!
!    kerminen, v., and m. kulmala, 2002,
!       analytical formulae connecting the "real" and the "apparent"
!       nucleation rate and the nuclei number concentration
!       for atmospheric nucleation events
!
!.......................................................................

! subr arguments (in)
    real(kind_phys), intent(in) :: dtnuc             ! nucleation time step (s)
    real(kind_phys), intent(in) :: temp_in           ! temperature, in k
    real(kind_phys), intent(in) :: rh_in             ! relative humidity, as fraction
    real(kind_phys), intent(in) :: press_in          ! air pressure (pa)
    real(kind_phys), intent(in) :: zm_in             ! layer midpoint height (m)
    real(kind_phys), intent(in) :: pblh_in           ! pbl height (m)

    real(kind_phys), intent(in) :: qh2so4_cur, qh2so4_avg
    ! gas h2so4 mixing ratios (mol/mol-air)
    real(kind_phys), intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
    ! qxxx_cur = current value (after gas chem and condensation)
    ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
    real(kind_phys), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)
    real(kind_phys), intent(in) :: mw_so4a_host      ! mw of so4 aerosol in host code (g/mol)

    integer, intent(in) :: newnuc_method_flagaa     ! 1=merikanto et al (2007) ternary
    ! 2=vehkamaki et al (2002) binary
    integer, intent(in) :: nsize                    ! number of aerosol size bins
    integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
    real(kind_phys), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
    real(kind_phys), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)
    integer, intent(in) :: ldiagaa

! subr arguments (out)
    integer, intent(out) :: isize_nuc         ! size bin into which new particles go
    real(kind_phys), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
    real(kind_phys), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
    real(kind_phys), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
    real(kind_phys), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
    real(kind_phys), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
    ! aerosol changes are > 0; gas changes are < 0
    real(kind_phys), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)

! subr arguments (out) passed via common block
!    these are used to duplicate the outputs of yang zhang's original test driver
!    they are not really needed in wrf-chem
    real(kind_phys) :: ratenuclt        ! j = ternary nucleation rate from napari param. (cm-3 s-1)
    real(kind_phys) :: rateloge         ! ln (j)
    real(kind_phys) :: cnum_h2so4       ! number of h2so4 molecules in the critical nucleus
    real(kind_phys) :: cnum_nh3         ! number of nh3   molecules in the critical nucleus
    real(kind_phys) :: cnum_tot         ! total number of molecules in the critical nucleus
    real(kind_phys) :: radius_cluster   ! the radius of cluster (nm)

! local variables
    integer :: i
    integer :: igrow
    integer, save :: icase = 0, icase_reldiffmax = 0
!       integer, parameter :: ldiagaa = -1
    integer :: lun
    integer :: newnuc_method_flagaa2

    real(kind_phys), parameter :: onethird = 1.0_kind_phys/3.0_kind_phys

    real(kind_phys), parameter :: accom_coef_h2so4 = 0.65_kind_phys   ! accomodation coef for h2so4 conden

! dry densities (kg/m3) molecular weights of aerosol
! ammsulf, ammbisulf, and sulfacid (from mosaic  dens_electrolyte values)
!       real(kind_phys), parameter :: dens_ammsulf   = 1.769e3
!       real(kind_phys), parameter :: dens_ammbisulf = 1.78e3
!       real(kind_phys), parameter :: dens_sulfacid  = 1.841e3
! use following to match cam3 modal_aero densities
    real(kind_phys), parameter :: dens_ammsulf = 1.770e3_kind_phys
    real(kind_phys), parameter :: dens_ammbisulf = 1.770e3_kind_phys
    real(kind_phys), parameter :: dens_sulfacid = 1.770e3_kind_phys

! molecular weights (g/mol) of aerosol ammsulf, ammbisulf, and sulfacid
!    for ammbisulf and sulfacid, use 114 & 96 here rather than 115 & 98
!    because we don't keep track of aerosol hion mass
    real(kind_phys), parameter :: mw_ammsulf = 132.0_kind_phys
    real(kind_phys), parameter :: mw_ammbisulf = 114.0_kind_phys
    real(kind_phys), parameter :: mw_sulfacid = 96.0_kind_phys

    real(kind_phys), save :: reldiffmax = 0.0_kind_phys

    real(kind_phys) :: cair                     ! dry-air molar density (mol/m3)
    real(kind_phys) :: cs_prime_kk              ! kk2002 "cs_prime" parameter (1/m2)
    real(kind_phys) :: cs_kk                    ! kk2002 "cs" parameter (1/s)
    real(kind_phys) :: dens_part                ! "grown" single-particle dry density (kg/m3)
    real(kind_phys) :: dfin_kk, dnuc_kk         ! kk2002 final/initial new particle wet diameter (nm)
    real(kind_phys) :: dpdry_clus               ! critical cluster diameter (m)
    real(kind_phys) :: dpdry_part               ! "grown" single-particle dry diameter (m)
    real(kind_phys) :: tmpa, tmpb, tmpc, tmpe, tmpq
    real(kind_phys) :: tmpa1, tmpb1
    real(kind_phys) :: tmp_m1, tmp_m2, tmp_m3, tmp_n1, tmp_n2, tmp_n3
    real(kind_phys) :: tmp_spd                  ! h2so4 vapor molecular speed (m/s)
    real(kind_phys) :: factor_kk
    real(kind_phys) :: fogas, foso4a, fonh4a, fonuma
    real(kind_phys) :: freduce                  ! reduction factor applied to nucleation rate
    ! due to limited availability of h2so4 & nh3 gases
    real(kind_phys) :: freducea, freduceb
    real(kind_phys) :: gamma_kk                 ! kk2002 "gamma" parameter (nm2*m2/h)
    real(kind_phys) :: gr_kk                    ! kk2002 "gr" parameter (nm/h)
    real(kind_phys) :: kgaero_per_moleso4a      ! (kg dry aerosol)/(mol aerosol so4)
    real(kind_phys) :: mass_part                ! "grown" single-particle dry mass (kg)
    real(kind_phys) :: molenh4a_per_moleso4a    ! (mol aerosol nh4)/(mol aerosol so4)
    real(kind_phys) :: nh3ppt, nh3ppt_bb        ! actual and bounded nh3 (ppt)
    real(kind_phys) :: nu_kk                    ! kk2002 "nu" parameter (nm)
    real(kind_phys) :: qmolnh4a_del_max         ! max production of aerosol nh4 over dtnuc (mol/mol-air)
    real(kind_phys) :: qmolso4a_del_max         ! max production of aerosol so4 over dtnuc (mol/mol-air)
    real(kind_phys) :: ratenuclt_bb             ! nucleation rate (#/m3/s)
    real(kind_phys) :: ratenuclt_kk             ! nucleation rate after kk2002 adjustment (#/m3/s)
    real(kind_phys) :: rh_bb                    ! bounded value of rh_in
    real(kind_phys) :: so4vol_in                ! concentration of h2so4 for nucl. calc., molecules cm-3
    real(kind_phys) :: so4vol_bb                ! bounded value of so4vol_in
    real(kind_phys) :: temp_bb                  ! bounded value of temp_in
    real(kind_phys) :: voldry_clus              ! critical-cluster dry volume (m3)
    real(kind_phys) :: voldry_part              ! "grown" single-particle dry volume (m3)
    real(kind_phys) :: wetvol_dryvol            ! grown particle (wet-volume)/(dry-volume)
    real(kind_phys) :: wet_volfrac_so4a         ! grown particle (dry-volume-from-so4)/(wet-volume)

!
! if h2so4 vapor < qh2so4_cutoff
! exit with new particle formation = 0
!
    isize_nuc = 1
    qnuma_del = 0.0_kind_phys
    qso4a_del = 0.0_kind_phys
    qnh4a_del = 0.0_kind_phys
    qh2so4_del = 0.0_kind_phys
    qnh3_del = 0.0_kind_phys
!       if (qh2so4_avg .le. qh2so4_cutoff) return   ! this no longer needed
!       if (qh2so4_cur .le. qh2so4_cutoff) return   ! this no longer needed

    if ((newnuc_method_flagaa /= 1) .and. &
        (newnuc_method_flagaa /= 2) .and. &
        (newnuc_method_flagaa /= 11) .and. &
        (newnuc_method_flagaa /= 12)) return

!
! make call to parameterization routine
!

! calc h2so4 in molecules/cm3 and nh3 in ppt
    cair = press_in/(temp_in*rgas)
    so4vol_in = qh2so4_avg*cair*avogad*1.0e-6_kind_phys
    nh3ppt = qnh3_cur*1.0e12_kind_phys
    ratenuclt = 1.0e-38_kind_phys
    rateloge = log(ratenuclt)

    if ((newnuc_method_flagaa /= 2) .and. &
        (nh3ppt >= 0.1_kind_phys)) then
! make call to merikanto ternary parameterization routine
! (when nh3ppt < 0.1, use binary param instead)

      if (so4vol_in >= 5.0e4_kind_phys) then
        temp_bb = max(235.0_kind_phys, min(295.0_kind_phys, temp_in))
        rh_bb = max(0.05_kind_phys, min(0.95_kind_phys, rh_in))
        so4vol_bb = max(5.0e4_kind_phys, min(1.0e9_kind_phys, so4vol_in))
        nh3ppt_bb = max(0.1_kind_phys, min(1.0e3_kind_phys, nh3ppt))
        call ternary_nuc_merik2007( &
          temp_bb, rh_bb, so4vol_bb, nh3ppt_bb, &
          rateloge, &
          cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster)
      end if
      newnuc_method_flagaa2 = 1

    else
! make call to vehkamaki binary parameterization routine

      if (so4vol_in >= 1.0e4_kind_phys) then
        temp_bb = max(230.15_kind_phys, min(305.15_kind_phys, temp_in))
        rh_bb = max(1.0e-4_kind_phys, min(1.0_kind_phys, rh_in))
        so4vol_bb = max(1.0e4_kind_phys, min(1.0e11_kind_phys, so4vol_in))
        call binary_nuc_vehk2002( &
          temp_bb, rh_bb, so4vol_bb, &
          ratenuclt, rateloge, &
          cnum_h2so4, cnum_tot, radius_cluster)
      end if
      cnum_nh3 = 0.0_kind_phys
      newnuc_method_flagaa2 = 2

    end if

! do boundary layer nuc
    if ((newnuc_method_flagaa == 11) .or. &
        (newnuc_method_flagaa == 12)) then
      if (zm_in <= max(pblh_in, 100.0_kind_phys)) then
        so4vol_bb = so4vol_in
        call pbl_nuc_wang2008(so4vol_bb, &
                              newnuc_method_flagaa, newnuc_method_flagaa2, &
                              ratenuclt, rateloge, &
                              cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster)
      end if
    end if

! if nucleation rate is less than 1e-6 #/m3/s ~= 0.1 #/cm3/day,
! exit with new particle formation = 0
    if (rateloge <= -13.82_kind_phys) return
!       if (ratenuclt .le. 1.0e-6) return
    ratenuclt = exp(rateloge)
    ratenuclt_bb = ratenuclt*1.0e6_kind_phys

! wet/dry volume ratio - use simple kohler approx for ammsulf/ammbisulf
    tmpa = max(0.10_kind_phys, min(0.95_kind_phys, rh_in))
    wetvol_dryvol = 1.0_kind_phys - 0.56_kind_phys/log(tmpa)

! determine size bin into which the new particles go
! (probably it will always be bin #1, but ...)
    voldry_clus = (max(cnum_h2so4, 1.0_kind_phys)*mw_so4a + cnum_nh3*mw_nh4a)/ &
                  (1.0e3_kind_phys*dens_sulfacid*avogad)
! correction when host code sulfate is really ammonium bisulfate/sulfate
    voldry_clus = voldry_clus*(mw_so4a_host/mw_so4a)
    dpdry_clus = (voldry_clus*6.0_kind_phys/pi)**onethird

    isize_nuc = 1
    dpdry_part = dplom_sect(1)
    if (dpdry_clus <= dplom_sect(1)) then
      igrow = 1   ! need to clusters to larger size
    else if (dpdry_clus >= dphim_sect(nsize)) then
      igrow = 0
      isize_nuc = nsize
      dpdry_part = dphim_sect(nsize)
    else
      igrow = 0
      do i = 1, nsize
        if (dpdry_clus < dphim_sect(i)) then
          isize_nuc = i
          dpdry_part = dpdry_clus
          dpdry_part = min(dpdry_part, dphim_sect(i))
          dpdry_part = max(dpdry_part, dplom_sect(i))
          exit
        end if
      end do
    end if
    voldry_part = (pi/6.0_kind_phys)*(dpdry_part**3)

!
! determine composition and density of the "grown particles"
! the grown particles are assumed to be liquid
!    (since critical clusters contain water)
!    so any (nh4/so4) molar ratio between 0 and 2 is allowed
! assume that the grown particles will have
!    (nh4/so4 molar ratio) = min( 2, (nh3/h2so4 gas molar ratio) )
!
    if (igrow <= 0) then
! no "growing" so pure sulfuric acid
      tmp_n1 = 0.0_kind_phys
      tmp_n2 = 0.0_kind_phys
      tmp_n3 = 1.0_kind_phys
    else if (qnh3_cur >= qh2so4_cur) then
! combination of ammonium sulfate and ammonium bisulfate
! tmp_n1 & tmp_n2 = mole fractions of the ammsulf & ammbisulf
      tmp_n1 = (qnh3_cur/qh2so4_cur) - 1.0_kind_phys
      tmp_n1 = max(0.0_kind_phys, min(1.0_kind_phys, tmp_n1))
      tmp_n2 = 1.0_kind_phys - tmp_n1
      tmp_n3 = 0.0_kind_phys
    else
! combination of ammonium bisulfate and sulfuric acid
! tmp_n2 & tmp_n3 = mole fractions of the ammbisulf & sulfacid
      tmp_n1 = 0.0_kind_phys
      tmp_n2 = (qnh3_cur/qh2so4_cur)
      tmp_n2 = max(0.0_kind_phys, min(1.0_kind_phys, tmp_n2))
      tmp_n3 = 1.0_kind_phys - tmp_n2
    end if

    tmp_m1 = tmp_n1*mw_ammsulf
    tmp_m2 = tmp_n2*mw_ammbisulf
    tmp_m3 = tmp_n3*mw_sulfacid
    dens_part = (tmp_m1 + tmp_m2 + tmp_m3)/ &
                ((tmp_m1/dens_ammsulf) + (tmp_m2/dens_ammbisulf) &
                 + (tmp_m3/dens_sulfacid))
    dens_nh4so4a = dens_part
    mass_part = voldry_part*dens_part
! (mol aerosol nh4)/(mol aerosol so4)
    molenh4a_per_moleso4a = 2.0_kind_phys*tmp_n1 + tmp_n2
! (kg dry aerosol)/(mol aerosol so4)
    kgaero_per_moleso4a = 1.0e-3_kind_phys*(tmp_m1 + tmp_m2 + tmp_m3)
! correction when host code sulfate is really ammonium bisulfate/sulfate
    kgaero_per_moleso4a = kgaero_per_moleso4a*(mw_so4a_host/mw_so4a)

! fraction of wet volume due to so4a
    tmpb = 1.0_kind_phys + molenh4a_per_moleso4a*17.0_kind_phys/98.0_kind_phys
    wet_volfrac_so4a = 1.0_kind_phys/(wetvol_dryvol*tmpb)

!
! calc kerminen & kulmala (2002) correction
!
    if (igrow <= 0) then
      factor_kk = 1.0_kind_phys

    else
! "gr" parameter (nm/h) = condensation growth rate of new particles
! use kk2002 eqn 21 for h2so4 uptake, and correct for nh3 & h2o uptake
      tmp_spd = 14.7_kind_phys*sqrt(temp_in)   ! h2so4 molecular speed (m/s)
      gr_kk = 3.0e-9_kind_phys*tmp_spd*mw_sulfacid*so4vol_in/ &
              (dens_part*wet_volfrac_so4a)

! "gamma" parameter (nm2/m2/h)
! use kk2002 eqn 22
!
! dfin_kk = wet diam (nm) of grown particle having dry dia = dpdry_part (m)
      dfin_kk = 1.0e9_kind_phys*dpdry_part*(wetvol_dryvol**onethird)
! dnuc_kk = wet diam (nm) of cluster
      dnuc_kk = 2.0_kind_phys*radius_cluster
      dnuc_kk = max(dnuc_kk, 1.0_kind_phys)
! neglect (dmean/150)**0.048 factor,
! which should be very close to 1.0 because of small exponent
      gamma_kk = 0.23_kind_phys*(dnuc_kk)**0.2_kind_phys &
                 *(dfin_kk/3.0_kind_phys)**0.075_kind_phys &
                 *(dens_part*1.0e-3_kind_phys)**(-0.33_kind_phys) &
                 *(temp_in/293.0_kind_phys)**(-0.75_kind_phys)

! "cs_prime parameter" (1/m2)
! instead kk2002 eqn 3, use
!     cs_prime ~= tmpa / (4*pi*tmpb * h2so4_accom_coef)
! where
!     tmpa = -d(ln(h2so4))/dt by conden to particles   (1/h units)
!     tmpb = h2so4 vapor diffusivity (m2/h units)
! this approx is generally within a few percent of the cs_prime
!     calculated directly from eqn 2,
!     which is acceptable, given overall uncertainties
! tmpa = -d(ln(h2so4))/dt by conden to particles   (1/h units)
      tmpa = h2so4_uptkrate*3600.0_kind_phys
      tmpa1 = tmpa
      tmpa = max(tmpa, 0.0_kind_phys)
! tmpb = h2so4 gas diffusivity (m2/s, then m2/h)
      tmpb = 6.7037e-6_kind_phys*(temp_in**0.75_kind_phys)/cair
      tmpb1 = tmpb         ! m2/s
      tmpb = tmpb*3600.0_kind_phys   ! m2/h
      cs_prime_kk = tmpa/(4.0_kind_phys*pi*tmpb*accom_coef_h2so4)
      cs_kk = cs_prime_kk*4.0_kind_phys*pi*tmpb1

! "nu" parameter (nm) -- kk2002 eqn 11
      nu_kk = gamma_kk*cs_prime_kk/gr_kk
! nucleation rate adjustment factor (--) -- kk2002 eqn 13
      factor_kk = exp((nu_kk/dfin_kk) - (nu_kk/dnuc_kk))

    end if
    ratenuclt_kk = ratenuclt_bb*factor_kk

! max production of aerosol dry mass (kg-aero/m3-air)
    tmpa = max(0.0_kind_phys, (ratenuclt_kk*dtnuc*mass_part))
! max production of aerosol so4 (mol-so4a/mol-air)
    tmpe = tmpa/(kgaero_per_moleso4a*cair)
! max production of aerosol so4 (mol/mol-air)
! based on ratenuclt_kk and mass_part
    qmolso4a_del_max = tmpe

! check if max production exceeds available h2so4 vapor
    freducea = 1.0_kind_phys
    if (qmolso4a_del_max > qh2so4_cur) then
      freducea = qh2so4_cur/qmolso4a_del_max
    end if

! check if max production exceeds available nh3 vapor
    freduceb = 1.0_kind_phys
    if (molenh4a_per_moleso4a >= 1.0e-10_kind_phys) then
! max production of aerosol nh4 (ppm) based on ratenuclt_kk and mass_part
      qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
      if (qmolnh4a_del_max > qnh3_cur) then
        freduceb = qnh3_cur/qmolnh4a_del_max
      end if
    end if
    freduce = min(freducea, freduceb)

! if adjusted nucleation rate is less than 1e-12 #/m3/s ~= 0.1 #/cm3/day,
! exit with new particle formation = 0
    if (freduce*ratenuclt_kk <= 1.0e-12_kind_phys) return

! note:  suppose that at this point, freduce < 1.0 (no gas-available
!    constraints) and molenh4a_per_moleso4a < 2.0
! if the gas-available constraints is do to h2so4 availability,
!    then it would be possible to condense "additional" nh3 and have
!    (nh3/h2so4 gas molar ratio) < (nh4/so4 aerosol molar ratio) <= 2
! one could do some additional calculations of
!    dens_part & molenh4a_per_moleso4a to realize this
! however, the particle "growing" is a crude approximate way to get
!    the new particles to the host code's minimum particle size,
! are such refinements worth the effort?

! changes to h2so4 & nh3 gas (in mol/mol-air), limited by amounts available
    tmpa = 0.9999_kind_phys
    qh2so4_del = min(tmpa*qh2so4_cur, freduce*qmolso4a_del_max)
    qnh3_del = min(tmpa*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a)
    qh2so4_del = -qh2so4_del
    qnh3_del = -qnh3_del

! changes to so4 & nh4 aerosol (in mol/mol-air)
    qso4a_del = -qh2so4_del
    qnh4a_del = -qnh3_del
! change to aerosol number (in #/mol-air)
    qnuma_del = 1.0e-3_kind_phys*(qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part

! do the following (tmpa, tmpb, tmpc) calculations as a check
! max production of aerosol number (#/mol-air)
    tmpa = max(0.0_kind_phys, (ratenuclt_kk*dtnuc/cair))
! adjusted production of aerosol number (#/mol-air)
    tmpb = tmpa*freduce
! relative difference from qnuma_del
    tmpc = (tmpb - qnuma_del)/max(tmpb, qnuma_del, 1.0e-35_kind_phys)

!
! diagnostic output to fort.41
! (this should be commented-out or deleted in the wrf-chem version)
!
    if (ldiagaa <= 0) return

    icase = icase + 1
    if (abs(tmpc) > abs(reldiffmax)) then
      reldiffmax = tmpc
      icase_reldiffmax = icase
    end if
!       do lun = 41, 51, 10
    do lun = 6, 6
!          write(lun,'(/)')
      write (lun, '(a,2i9,1p,e10.2)') &
        'vehkam bin-nuc icase, icase_rdmax =', &
        icase, icase_reldiffmax, reldiffmax
      if (freduceb < freducea) then
        if (abs(freducea - freduceb) > &
            3.0e-7_kind_phys*max(freduceb, freducea)) write (lun, '(a,1p,2e15.7)') &
          'freducea, b =', freducea, freduceb
      end if
    end do

! output factors so that output matches that of ternucl03
!       fogas  = 1.0e6                     ! convert mol/mol-air to ppm
!       foso4a = 1.0e9*mw_so4a/mw_air      ! convert mol-so4a/mol-air to ug/kg-air
!       fonh4a = 1.0e9*mw_nh4a/mw_air      ! convert mol-nh4a/mol-air to ug/kg-air
!       fonuma = 1.0e3/mw_air              ! convert #/mol-air to #/kg-air
    fogas = 1.0_kind_phys
    foso4a = 1.0_kind_phys
    fonh4a = 1.0_kind_phys
    fonuma = 1.0_kind_phys

!       do lun = 41, 51, 10
    do lun = 6, 6

      write (lun, '(a,2i5)') 'newnuc_method_flagaa/aa2', &
        newnuc_method_flagaa, newnuc_method_flagaa2

      write (lun, 9210)
      write (lun, 9201) temp_in, rh_in, &
        ratenuclt, 2.0_kind_phys*radius_cluster*1.0e-7_kind_phys, dpdry_part*1.0e2_kind_phys, &
        voldry_part*1.0e6_kind_phys, float(igrow)
      write (lun, 9215)
      write (lun, 9201) &
        qh2so4_avg*fogas, 0.0_kind_phys, &
        qh2so4_cur*fogas, qnh3_cur*fogas, &
        qh2so4_del*fogas, qnh3_del*fogas, &
        qso4a_del*foso4a, qnh4a_del*fonh4a

      write (lun, 9220)
      write (lun, 9201) &
        dtnuc, dens_nh4so4a*1.0e-3_kind_phys, &
        (qnh3_cur/qh2so4_cur), molenh4a_per_moleso4a, &
        qnuma_del*fonuma, tmpb*fonuma, tmpc, freduce

    end do

!       lun = 51
    lun = 6
    write (lun, 9230)
    write (lun, 9201) &
      press_in, cair*1.0e-6_kind_phys, so4vol_in, &
      wet_volfrac_so4a, wetvol_dryvol, dens_part*1.0e-3_kind_phys

    if (igrow > 0) then
      write (lun, 9240)
      write (lun, 9201) &
        tmp_spd, gr_kk, dnuc_kk, dfin_kk, &
        gamma_kk, tmpa1, tmpb1, cs_kk

      write (lun, 9250)
      write (lun, 9201) &
        cs_prime_kk, nu_kk, factor_kk, ratenuclt, &
        ratenuclt_kk*1.0e-6_kind_phys
    end if

9201 format(1p, 40e10.2)
9210 format( &
      '      temp        rh', &
      '   ratenuc  dia_clus ddry_part', &
      ' vdry_part     igrow')
9215 format( &
      '  h2so4avg  h2so4pre', &
      '  h2so4cur   nh3_cur', &
      '  h2so4del   nh3_del', &
      '  so4a_del  nh4a_del')
9220 format( &
      '     dtnuc    dens_a   nh/so g   nh/so a', &
      '  numa_del  numa_dl2   reldiff   freduce')
9230 format( &
      '  press_in      cair so4_volin', &
      ' wet_volfr wetv_dryv dens_part')
9240 format( &
      '   tmp_spd     gr_kk   dnuc_kk   dfin_kk', &
      '  gamma_kk     tmpa1     tmpb1     cs_kk')
9250 format( &
      ' cs_pri_kk     nu_kk factor_kk ratenuclt', &
      ' ratenu_kk')

  end subroutine mer07_veh02_nuc_mosaic_1box

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine pbl_nuc_wang2008(so4vol, &
                              newnuc_method_flagaa, newnuc_method_flagaa2, &
                              ratenucl, rateloge, &
                              cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster)
!
! calculates boundary nucleation nucleation rate
! using the first or second-order parameterization in
!     wang, m., and j.e. penner, 2008,
!        aerosol indirect forcing in a global model with particle nucleation,
!        atmos. chem. phys. discuss., 8, 13943-13998
!

! subr arguments (in)
    real(kind_phys), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)
    integer, intent(in)  :: newnuc_method_flagaa
    ! [11,12] value selects [first,second]-order parameterization

! subr arguments (inout)
    integer, intent(inout)  :: newnuc_method_flagaa2
    real(kind_phys), intent(inout) :: ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
    real(kind_phys), intent(inout) :: rateloge         ! log( ratenucl )

    real(kind_phys), intent(inout) :: cnum_tot         ! total number of molecules
    ! in the critical nucleus
    real(kind_phys), intent(inout) :: cnum_h2so4       ! number of h2so4 molecules
    real(kind_phys), intent(inout) :: cnum_nh3         ! number of nh3 molecules
    real(kind_phys), intent(inout) :: radius_cluster   ! the radius of cluster (nm)

! local variables
    real(kind_phys) :: tmp_diam, tmp_mass, tmp_volu
    real(kind_phys) :: tmp_rateloge, tmp_ratenucl

! executable

! nucleation rate
    if (newnuc_method_flagaa == 11) then
      tmp_ratenucl = 1.0e-6_kind_phys*so4vol
    else if (newnuc_method_flagaa == 12) then
      tmp_ratenucl = 1.0e-12_kind_phys*(so4vol**2)
    else
      return
    end if
    tmp_rateloge = log(tmp_ratenucl)

! exit if pbl nuc rate is lower than (incoming) ternary/binary rate
    if (tmp_rateloge <= rateloge) return

    rateloge = tmp_rateloge
    ratenucl = tmp_ratenucl
    newnuc_method_flagaa2 = newnuc_method_flagaa

! following wang 2002, assume fresh nuclei are 1 nm diameter
!    subsequent code will "grow" them to aitken mode size
    radius_cluster = 0.5_kind_phys

! assume fresh nuclei are pure h2so4
!    since aitken size >> initial size, the initial composition
!    has very little impact on the results
    tmp_diam = radius_cluster*2.0e-7_kind_phys   ! diameter in cm
    tmp_volu = (tmp_diam**3)*(pi/6.0_kind_phys)  ! volume in cm^3
    tmp_mass = tmp_volu*1.8_kind_phys            ! mass in g
    cnum_h2so4 = (tmp_mass/98.0_kind_phys)*6.023e23_kind_phys   ! no. of h2so4 molec assuming pure h2so4
    cnum_tot = cnum_h2so4
    cnum_nh3 = 0.0_kind_phys

  end subroutine pbl_nuc_wang2008

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine binary_nuc_vehk2002(temp, rh, so4vol, &
                                 ratenucl, rateloge, &
                                 cnum_h2so4, cnum_tot, radius_cluster)
!
! calculates binary nucleation rate and critical cluster size
! using the parameterization in
!     vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!        c. timmreck, m. noppel and a. laaksonen, 2002,
!        an improved parameterization for sulfuric acid-water nucleation
!        rates for tropospheric and stratospheric conditions,
!        j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
!

! subr arguments (in)
    real(kind_phys), intent(in) :: temp              ! temperature (k)
    real(kind_phys), intent(in) :: rh                ! relative humidity (0-1)
    real(kind_phys), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)

! subr arguments (out)
    real(kind_phys), intent(out) :: ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
    real(kind_phys), intent(out) :: rateloge         ! log( ratenucl )

    real(kind_phys), intent(out) :: cnum_h2so4       ! number of h2so4 molecules
    ! in the critical nucleus
    real(kind_phys), intent(out) :: cnum_tot         ! total number of molecules
    ! in the critical nucleus
    real(kind_phys), intent(out) :: radius_cluster   ! the radius of cluster (nm)

! local variables
    real(kind_phys) :: crit_x
    real(kind_phys) :: acoe, bcoe, ccoe, dcoe, ecoe, fcoe, gcoe, hcoe, icoe, jcoe
    real(kind_phys) :: tmpa, tmpb

! executable

! calc sulfuric acid mole fraction in critical cluster
    crit_x = 0.740997_kind_phys - 0.00266379_kind_phys*temp &
             - 0.00349998_kind_phys*log(so4vol) &
             + 0.0000504022_kind_phys*temp*log(so4vol) &
             + 0.00201048_kind_phys*log(rh) &
             - 0.000183289_kind_phys*temp*log(rh) &
             + 0.00157407_kind_phys*(log(rh))**2.0_kind_phys &
             - 0.0000179059_kind_phys*temp*(log(rh))**2.0_kind_phys &
             + 0.000184403_kind_phys*(log(rh))**3.0_kind_phys &
             - 1.50345e-6_kind_phys*temp*(log(rh))**3.0_kind_phys

! calc nucleation rate
    acoe = 0.14309_kind_phys + 2.21956_kind_phys*temp &
           - 0.0273911_kind_phys*temp**2.0_kind_phys &
           + 0.0000722811_kind_phys*temp**3.0_kind_phys + 5.91822_kind_phys/crit_x

    bcoe = 0.117489_kind_phys + 0.462532_kind_phys*temp &
           - 0.0118059_kind_phys*temp**2.0_kind_phys &
           + 0.0000404196_kind_phys*temp**3.0_kind_phys + 15.7963_kind_phys/crit_x

    ccoe = -0.215554_kind_phys - 0.0810269_kind_phys*temp &
           + 0.00143581_kind_phys*temp**2.0_kind_phys &
           - 4.7758e-6_kind_phys*temp**3.0_kind_phys &
           - 2.91297_kind_phys/crit_x

    dcoe = -3.58856_kind_phys + 0.049508_kind_phys*temp &
           - 0.00021382_kind_phys*temp**2.0_kind_phys &
           + 3.10801e-7_kind_phys*temp**3.0_kind_phys &
           - 0.0293333_kind_phys/crit_x

    ecoe = 1.14598_kind_phys - 0.600796_kind_phys*temp &
           + 0.00864245_kind_phys*temp**2.0_kind_phys &
           - 0.0000228947_kind_phys*temp**3.0_kind_phys &
           - 8.44985_kind_phys/crit_x

    fcoe = 2.15855_kind_phys + 0.0808121_kind_phys*temp &
           - 0.000407382_kind_phys*temp**2.0_kind_phys &
           - 4.01957e-7_kind_phys*temp**3.0_kind_phys &
           + 0.721326_kind_phys/crit_x

    gcoe = 1.6241_kind_phys - 0.0160106_kind_phys*temp &
           + 0.0000377124_kind_phys*temp**2.0_kind_phys &
           + 3.21794e-8_kind_phys*temp**3.0_kind_phys &
           - 0.0113255_kind_phys/crit_x

    hcoe = 9.71682_kind_phys - 0.115048_kind_phys*temp &
           + 0.000157098_kind_phys*temp**2.0_kind_phys &
           + 4.00914e-7_kind_phys*temp**3.0_kind_phys &
           + 0.71186_kind_phys/crit_x

    icoe = -1.05611_kind_phys + 0.00903378_kind_phys*temp &
           - 0.0000198417_kind_phys*temp**2.0_kind_phys &
           + 2.46048e-8_kind_phys*temp**3.0_kind_phys &
           - 0.0579087_kind_phys/crit_x

    jcoe = -0.148712_kind_phys + 0.00283508_kind_phys*temp &
           - 9.24619e-6_kind_phys*temp**2.0_kind_phys &
           + 5.00427e-9_kind_phys*temp**3.0_kind_phys &
           - 0.0127081_kind_phys/crit_x

    tmpa = ( &
           acoe &
           + bcoe*log(rh) &
           + ccoe*(log(rh))**2.0_kind_phys &
           + dcoe*(log(rh))**3.0_kind_phys &
           + ecoe*log(so4vol) &
           + fcoe*(log(rh))*(log(so4vol)) &
           + gcoe*((log(rh))**2.0_kind_phys) &
           *(log(so4vol)) &
           + hcoe*(log(so4vol))**2.0_kind_phys &
           + icoe*log(rh) &
           *((log(so4vol))**2.0_kind_phys) &
           + jcoe*(log(so4vol))**3.0_kind_phys &
           )
    rateloge = tmpa
    tmpa = min(tmpa, log(1.0e38_kind_phys))
    ratenucl = exp(tmpa)
!       write(*,*) 'tmpa, ratenucl =', tmpa, ratenucl

! calc number of molecules in critical cluster
    acoe = -0.00295413_kind_phys - 0.0976834_kind_phys*temp &
           + 0.00102485_kind_phys*temp**2.0_kind_phys &
           - 2.18646e-6_kind_phys*temp**3.0_kind_phys - 0.101717_kind_phys/crit_x

    bcoe = -0.00205064_kind_phys - 0.00758504_kind_phys*temp &
           + 0.000192654_kind_phys*temp**2.0_kind_phys &
           - 6.7043e-7_kind_phys*temp**3.0_kind_phys - 0.255774_kind_phys/crit_x

    ccoe = +0.00322308_kind_phys + 0.000852637_kind_phys*temp &
           - 0.0000154757_kind_phys*temp**2.0_kind_phys &
           + 5.66661e-8_kind_phys*temp**3.0_kind_phys &
           + 0.0338444_kind_phys/crit_x

    dcoe = +0.0474323_kind_phys - 0.000625104_kind_phys*temp &
           + 2.65066e-6_kind_phys*temp**2.0_kind_phys &
           - 3.67471e-9_kind_phys*temp**3.0_kind_phys &
           - 0.000267251_kind_phys/crit_x

    ecoe = -0.0125211_kind_phys + 0.00580655_kind_phys*temp &
           - 0.000101674_kind_phys*temp**2.0_kind_phys &
           + 2.88195e-7_kind_phys*temp**3.0_kind_phys &
           + 0.0942243_kind_phys/crit_x

    fcoe = -0.038546_kind_phys - 0.000672316_kind_phys*temp &
           + 2.60288e-6_kind_phys*temp**2.0_kind_phys &
           + 1.19416e-8_kind_phys*temp**3.0_kind_phys &
           - 0.00851515_kind_phys/crit_x

    gcoe = -0.0183749_kind_phys + 0.000172072_kind_phys*temp &
           - 3.71766e-7_kind_phys*temp**2.0_kind_phys &
           - 5.14875e-10_kind_phys*temp**3.0_kind_phys &
           + 0.00026866_kind_phys/crit_x

    hcoe = -0.0619974_kind_phys + 0.000906958_kind_phys*temp &
           - 9.11728e-7_kind_phys*temp**2.0_kind_phys &
           - 5.36796e-9_kind_phys*temp**3.0_kind_phys &
           - 0.00774234_kind_phys/crit_x

    icoe = +0.0121827_kind_phys - 0.00010665_kind_phys*temp &
           + 2.5346e-7_kind_phys*temp**2.0_kind_phys &
           - 3.63519e-10_kind_phys*temp**3.0_kind_phys &
           + 0.000610065_kind_phys/crit_x

    jcoe = +0.000320184_kind_phys - 0.0000174762_kind_phys*temp &
           + 6.06504e-8_kind_phys*temp**2.0_kind_phys &
           - 1.4177e-11_kind_phys*temp**3.0_kind_phys &
           + 0.000135751_kind_phys/crit_x

    cnum_tot = exp( &
               acoe &
               + bcoe*log(rh) &
               + ccoe*(log(rh))**2.0_kind_phys &
               + dcoe*(log(rh))**3.0_kind_phys &
               + ecoe*log(so4vol) &
               + fcoe*(log(rh))*(log(so4vol)) &
               + gcoe*((log(rh))**2.0_kind_phys) &
               *(log(so4vol)) &
               + hcoe*(log(so4vol))**2.0_kind_phys &
               + icoe*log(rh) &
               *((log(so4vol))**2.0_kind_phys) &
               + jcoe*(log(so4vol))**3.0_kind_phys &
               )

    cnum_h2so4 = cnum_tot*crit_x

!   calc radius (nm) of critical cluster
    radius_cluster = exp(-1.6524245_kind_phys + 0.42316402_kind_phys*crit_x &
                         + 0.3346648_kind_phys*log(cnum_tot))

  end subroutine binary_nuc_vehk2002

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_newnuc_init(l_h2so4_in, l_nh3_in, &
                                    lnumait_in, lnh4ait_in, lso4ait_in, &
                                    mw_so4a_host_in, mw_nh4a_host_in, dens_so4a_host_in, &
                                    dgnum_aitken_in, dgnumhi_aitken_in, dgnumlo_aitken_in, &
                                    pi_in, rgas_in, avogad_in, mw_so4a_in, mw_nh4a_in, &
                                    r_universal_in, errmsg, errflg)

!-----------------------------------------------------------------------
!
! Purpose:
!    store the resolved species indices, aitken-mode so4/nh4 properties,
!       and host physical constants used by modal_aero_newnuc_run
!    index resolution and history-field registration are host
!       responsibilities (see modal_aero_newnuc_cam)
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! arguments
    integer, intent(in) :: l_h2so4_in          ! h2so4 gas index (host constituent space)
    integer, intent(in) :: l_nh3_in            ! nh3 gas index (host constituent space)
    integer, intent(in) :: lnumait_in          ! aitken number index (host constituent space)
    integer, intent(in) :: lnh4ait_in          ! aitken nh4 index (host constituent space)
    integer, intent(in) :: lso4ait_in          ! aitken so4 index (host constituent space)
    real(kind_phys), intent(in) :: mw_so4a_host_in     ! mw of so4 aerosol in host code (g/mol)
    real(kind_phys), intent(in) :: mw_nh4a_host_in     ! mw of nh4 aerosol in host code (g/mol)
    real(kind_phys), intent(in) :: dens_so4a_host_in   ! dry density of so4 aerosol in host code (kg/m3)
    real(kind_phys), intent(in) :: dgnum_aitken_in     ! aitken mode nominal dry diameter (m)
    real(kind_phys), intent(in) :: dgnumhi_aitken_in   ! aitken mode upper dry-diameter limit (m)
    real(kind_phys), intent(in) :: dgnumlo_aitken_in   ! aitken mode lower dry-diameter limit (m)
    real(kind_phys), intent(in) :: pi_in
    real(kind_phys), intent(in) :: rgas_in             ! gas constant (J/K/mol)
    real(kind_phys), intent(in) :: avogad_in           ! Avogadro's number (1/mol)
    real(kind_phys), intent(in) :: mw_so4a_in          ! molecular weight of sulfate (g/mol)
    real(kind_phys), intent(in) :: mw_nh4a_in          ! molecular weight of ammonium (g/mol)
    real(kind_phys), intent(in) :: r_universal_in      ! universal gas constant (J/K/kmol)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    !-----------------------------------------------------------------------

    errmsg = ' '
    errflg = 0

    l_h2so4_sv = l_h2so4_in
    l_nh3_sv = l_nh3_in
    lnumait_sv = lnumait_in
    lnh4ait_sv = lnh4ait_in
    lso4ait_sv = lso4ait_in

    mw_so4a_host = mw_so4a_host_in
    mw_nh4a_host = mw_nh4a_host_in
    dens_so4a_host = dens_so4a_host_in

    dgnum_aitken = dgnum_aitken_in
    dgnumhi_aitken = dgnumhi_aitken_in
    dgnumlo_aitken = dgnumlo_aitken_in

    pi = pi_in
    rgas = rgas_in
    avogad = avogad_in
    mw_so4a = mw_so4a_in
    mw_nh4a = mw_nh4a_in
    r_universal = r_universal_in

  end subroutine modal_aero_newnuc_init

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine ternary_nuc_merik2007(t, rh, c2, c3, j_log, ntot, nacid, namm, r)
!subroutine ternary_fit(          t, rh, c2, c3, j_log, ntot, nacid, namm, r )
! *************************** ternary_fit.f90 ********************************
! joonas merikanto, 2006
!
! fortran 90 subroutine that calculates the parameterized composition
! and nucleation rate of critical clusters in h2o-h2so4-nh3 vapor
!
! warning: the fit should not be used outside its limits of validity
! (limits indicated below)
!
! in:
! t:     temperature (k), limits 235-295 k
! rh:    relative humidity as fraction (eg. 0.5=50%) limits 0.05-0.95
! c2:    sulfuric acid concentration (molecules/cm3) limits 5x10^4 - 10^9 molecules/cm3
! c3:    ammonia mixing ratio (ppt) limits 0.1 - 1000 ppt
!
! out:
! j_log: logarithm of nucleation rate (1/(s cm3))
! ntot:  total number of molecules in the critical cluster
! nacid: number of sulfuric acid molecules in the critical cluster
! namm:  number of ammonia molecules in the critical cluster
! r:     radius of the critical cluster (nm)
!  ****************************************************************************

    real(kind_phys), intent(in) :: t, rh, c2, c3
    real(kind_phys), intent(out) :: j_log, ntot, nacid, namm, r
    real(kind_phys) :: j, t_onset

    t_onset = 143.6002929064716_kind_phys + 1.0178856665693992_kind_phys*rh + &
              10.196398812974294_kind_phys*log(c2) - &
              0.1849879416839113_kind_phys*log(c2)**2 - 17.161783213150173_kind_phys*log(c3) + &
              (109.92469248546053_kind_phys*log(c3))/log(c2) + &
              0.7734119613144357_kind_phys*log(c2)*log(c3) - 0.15576469879527022_kind_phys*log(c3)**2

    if (t_onset > t) then

      j_log = -12.861848898625231_kind_phys + 4.905527742256349_kind_phys*c3 - 358.2337705052991_kind_phys*rh - &
              0.05463019231872484_kind_phys*c3*t + 4.8630382337426985_kind_phys*rh*t + &
              0.00020258394697064567_kind_phys*c3*t**2 - 0.02175548069741675_kind_phys*rh*t**2 - &
              2.502406532869512e-7_kind_phys*c3*t**3 + 0.00003212869941055865_kind_phys*rh*t**3 - &
              4.39129415725234e6_kind_phys/log(c2)**2 + (56383.93843154586_kind_phys*t)/log(c2)**2 - &
              (239.835990963361_kind_phys*t**2)/log(c2)**2 + &
              (0.33765136625580167_kind_phys*t**3)/log(c2)**2 - &
              (629.7882041830943_kind_phys*rh)/(c3**3*log(c2)) + &
              (7.772806552631709_kind_phys*rh*t)/(c3**3*log(c2)) - &
              (0.031974053936299256_kind_phys*rh*t**2)/(c3**3*log(c2)) + &
              (0.00004383764128775082_kind_phys*rh*t**3)/(c3**3*log(c2)) + &
              1200.472096232311_kind_phys*log(c2) - 17.37107890065621_kind_phys*t*log(c2) + &
              0.08170681335921742_kind_phys*t**2*log(c2) - 0.00012534476159729881_kind_phys*t**3*log(c2) - &
              14.833042158178936_kind_phys*log(c2)**2 + 0.2932631303555295_kind_phys*t*log(c2)**2 - &
              0.0016497524241142845_kind_phys*t**2*log(c2)**2 + &
              2.844074805239367e-6_kind_phys*t**3*log(c2)**2 - 231375.56676032578_kind_phys*log(c3) - &
              100.21645273730675_kind_phys*rh*log(c3) + 2919.2852552424706_kind_phys*t*log(c3) + &
              0.977886555834732_kind_phys*rh*t*log(c3) - 12.286497122264588_kind_phys*t**2*log(c3) - &
              0.0030511783284506377_kind_phys*rh*t**2*log(c3) + &
              0.017249301826661612_kind_phys*t**3*log(c3) + 2.967320346100855e-6_kind_phys*rh*t**3*log(c3) + &
              (2.360931724951942e6_kind_phys*log(c3))/log(c2) - &
              (29752.130254319443_kind_phys*t*log(c3))/log(c2) + &
              (125.04965118142027_kind_phys*t**2*log(c3))/log(c2) - &
              (0.1752996881934318_kind_phys*t**3*log(c3))/log(c2) + &
              5599.912337254629_kind_phys*log(c2)*log(c3) - 70.70896612937771_kind_phys*t*log(c2)*log(c3) + &
              0.2978801613269466_kind_phys*t**2*log(c2)*log(c3) - &
              0.00041866525019504_kind_phys*t**3*log(c2)*log(c3) + 75061.15281456841_kind_phys*log(c3)**2 - &
              931.8802278173565_kind_phys*t*log(c3)**2 + 3.863266220840964_kind_phys*t**2*log(c3)**2 - &
              0.005349472062284983_kind_phys*t**3*log(c3)**2 - &
              (732006.8180571689_kind_phys*log(c3)**2)/log(c2) + &
              (9100.06398573816_kind_phys*t*log(c3)**2)/log(c2) - &
              (37.771091915932004_kind_phys*t**2*log(c3)**2)/log(c2) + &
              (0.05235455395566905_kind_phys*t**3*log(c3)**2)/log(c2) - &
              1911.0303773001353_kind_phys*log(c2)*log(c3)**2 + &
              23.6903969622286_kind_phys*t*log(c2)*log(c3)**2 - &
              0.09807872005428583_kind_phys*t**2*log(c2)*log(c3)**2 + &
              0.00013564560238552576_kind_phys*t**3*log(c2)*log(c3)**2 - &
              3180.5610833308_kind_phys*log(c3)**3 + 39.08268568672095_kind_phys*t*log(c3)**3 - &
              0.16048521066690752_kind_phys*t**2*log(c3)**3 + &
              0.00022031380023793877_kind_phys*t**3*log(c3)**3 + &
              (40751.075322248245_kind_phys*log(c3)**3)/log(c2) - &
              (501.66977622013934_kind_phys*t*log(c3)**3)/log(c2) + &
              (2.063469732254135_kind_phys*t**2*log(c3)**3)/log(c2) - &
              (0.002836873785758324_kind_phys*t**3*log(c3)**3)/log(c2) + &
              2.792313345723013_kind_phys*log(c2)**2*log(c3)**3 - &
              0.03422552111802899_kind_phys*t*log(c2)**2*log(c3)**3 + &
              0.00014019195277521142_kind_phys*t**2*log(c2)**2*log(c3)**3 - &
              1.9201227328396297e-7_kind_phys*t**3*log(c2)**2*log(c3)**3 - &
              980.923146020468_kind_phys*log(rh) + 10.054155220444462_kind_phys*t*log(rh) - &
              0.03306644502023841_kind_phys*t**2*log(rh) + 0.000034274041225891804_kind_phys*t**3*log(rh) + &
              (16597.75554295064_kind_phys*log(rh))/log(c2) - &
              (175.2365504237746_kind_phys*t*log(rh))/log(c2) + &
              (0.6033215603167458_kind_phys*t**2*log(rh))/log(c2) - &
              (0.0006731787599587544_kind_phys*t**3*log(rh))/log(c2) - &
              89.38961120336789_kind_phys*log(c3)*log(rh) + 1.153344219304926_kind_phys*t*log(c3)*log(rh) - &
              0.004954549700267233_kind_phys*t**2*log(c3)*log(rh) + &
              7.096309866238719e-6_kind_phys*t**3*log(c3)*log(rh) + &
              3.1712136610383244_kind_phys*log(c3)**3*log(rh) - &
              0.037822330602328806_kind_phys*t*log(c3)**3*log(rh) + &
              0.0001500555743561457_kind_phys*t**2*log(c3)**3*log(rh) - &
              1.9828365865570703e-7_kind_phys*t**3*log(c3)**3*log(rh)

      j = exp(j_log)

      ntot = 57.40091052369212_kind_phys - 0.2996341884645408_kind_phys*t + &
             0.0007395477768531926_kind_phys*t**2 - &
             5.090604835032423_kind_phys*log(c2) + 0.011016634044531128_kind_phys*t*log(c2) + &
             0.06750032251225707_kind_phys*log(c2)**2 - 0.8102831333223962_kind_phys*log(c3) + &
             0.015905081275952426_kind_phys*t*log(c3) - 0.2044174683159531_kind_phys*log(c2)*log(c3) + &
             0.08918159167625832_kind_phys*log(c3)**2 - 0.0004969033586666147_kind_phys*t*log(c3)**2 + &
             0.005704394549007816_kind_phys*log(c3)**3 + 3.4098703903474368_kind_phys*log(j) - &
             0.014916956508210809_kind_phys*t*log(j) + 0.08459090011666293_kind_phys*log(c3)*log(j) - &
             0.00014800625143907616_kind_phys*t*log(c3)*log(j) + 0.00503804694656905_kind_phys*log(j)**2

      r = 3.2888553966535506e-10_kind_phys - 3.374171768439839e-12_kind_phys*t + &
          1.8347359507774313e-14_kind_phys*t**2 + 2.5419844298881856e-12_kind_phys*log(c2) - &
          9.498107643050827e-14_kind_phys*t*log(c2) + 7.446266520834559e-13_kind_phys*log(c2)**2 + &
          2.4303397746137294e-11_kind_phys*log(c3) + 1.589324325956633e-14_kind_phys*t*log(c3) - &
          2.034596219775266e-12_kind_phys*log(c2)*log(c3) - 5.59303954457172e-13_kind_phys*log(c3)**2 - &
          4.889507104645867e-16_kind_phys*t*log(c3)**2 + 1.3847024107506764e-13_kind_phys*log(c3)**3 + &
          4.141077193427042e-15_kind_phys*log(j) - 2.6813110884009767e-14_kind_phys*t*log(j) + &
          1.2879071621313094e-12_kind_phys*log(c3)*log(j) - &
          3.80352446061867e-15_kind_phys*t*log(c3)*log(j) - 1.8790172502456827e-14_kind_phys*log(j)**2

      nacid = -4.7154180661803595_kind_phys + 0.13436423483953885_kind_phys*t - &
              0.00047184686478816176_kind_phys*t**2 - &
              2.564010713640308_kind_phys*log(c2) + 0.011353312899114723_kind_phys*t*log(c2) + &
              0.0010801941974317014_kind_phys*log(c2)**2 + 0.5171368624197119_kind_phys*log(c3) - &
              0.0027882479896204665_kind_phys*t*log(c3) + 0.8066971907026886_kind_phys*log(c3)**2 - &
              0.0031849094214409335_kind_phys*t*log(c3)**2 - 0.09951184152927882_kind_phys*log(c3)**3 + &
              0.00040072788891745513_kind_phys*t*log(c3)**3 + 1.3276469271073974_kind_phys*log(j) - &
              0.006167654171986281_kind_phys*t*log(j) - 0.11061390967822708_kind_phys*log(c3)*log(j) + &
              0.0004367575329273496_kind_phys*t*log(c3)*log(j) + 0.000916366357266258_kind_phys*log(j)**2

      namm = 71.20073903979772_kind_phys - 0.8409600103431923_kind_phys*t + &
             0.0024803006590334922_kind_phys*t**2 + &
             2.7798606841602607_kind_phys*log(c2) - 0.01475023348171676_kind_phys*t*log(c2) + &
             0.012264508212031405_kind_phys*log(c2)**2 - 2.009926050440182_kind_phys*log(c3) + &
             0.008689123511431527_kind_phys*t*log(c3) - 0.009141180198955415_kind_phys*log(c2)*log(c3) + &
             0.1374122553905617_kind_phys*log(c3)**2 - 0.0006253227821679215_kind_phys*t*log(c3)**2 + &
             0.00009377332742098946_kind_phys*log(c3)**3 + 0.5202974341687757_kind_phys*log(j) - &
             0.002419872323052805_kind_phys*t*log(j) + 0.07916392322884074_kind_phys*log(c3)*log(j) - &
             0.0003021586030317366_kind_phys*t*log(c3)*log(j) + 0.0046977006608603395_kind_phys*log(j)**2

    else
! nucleation rate less that 5e-6, setting j_log arbitrary small
      j_log = -300._kind_phys
    end if

  end subroutine ternary_nuc_merik2007

!----------------------------------------------------------------------
end module modal_aero_newnuc

