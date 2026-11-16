! Modal aerosol coagulation
! RCE 07.04.13:  Adapted from MIRAGE2 code
module modal_aero_coag
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_coag_init
  public :: modal_aero_coag_run

  integer, protected, public :: pair_option_acoag = 0
  ! specifies pairs of modes for which coagulation is calculated
  ! (set by modal_aero_coag_init from the host; default 0 = do no coag)
  !   1 -- [aitken-->accum]
  !   2 -- [aitken-->accum], and [pcarbon-->accum]
  !   3 -- [aitken-->accum], [pcarbon-->accum],
  !        and [aitken-->pcarbon--(aging)-->accum]
  ! other -- do no coag
  integer, parameter, public :: maxpair_acoag = 10
  integer, protected, public :: maxspec_acoag != nspec_max

  ! coagulation-pair tables (host constituent-index space), handed to
  ! modal_aero_coag_init by the host wrapper (see modal_aero_coag_cam)
  integer, protected, public :: npair_acoag = 0
  integer, protected, public :: modefrm_acoag(maxpair_acoag)
  integer, protected, public :: modetoo_acoag(maxpair_acoag)
  integer, protected, public :: modetooeff_acoag(maxpair_acoag)
  integer, protected, public :: nspecfrm_acoag(maxpair_acoag)
  integer, allocatable, protected, public :: lspecfrm_acoag(:, :)
  integer, allocatable, protected, public :: lspectoo_acoag(:, :)

  integer :: ip_aitacc, ip_aitpca, ip_pcaacc
  real(kind_phys), allocatable :: fac_m2v_aitage(:), fac_m2v_pcarbon(:)

  ! mode metadata from the host (set by modal_aero_coag_init):
  integer :: ntot_amode
  integer :: modeptr_accum, modeptr_aitken, modeptr_pcarbon
  integer, allocatable :: numptr_amode(:), mprognum_amode(:), nspec_amode(:)
  integer, allocatable :: lmassptr_amode(:, :)
  real(kind_phys), allocatable :: alnsg_amode(:), sigmag_amode(:)

  real(kind_phys) :: r_universal   ! universal gas constant (J/K/kmol)
  real(kind_phys) :: p0            ! standard pressure (Pa)
  real(kind_phys) :: tmelt         ! freezing point of water (K)
  real(kind_phys) :: boltz         ! Boltzmann constant (J/K)

contains

  subroutine modal_aero_coag_init(pair_option_acoag_in, &
                                  npair_acoag_in, modefrm_acoag_in, modetoo_acoag_in, &
                                  modetooeff_acoag_in, nspecfrm_acoag_in, &
                                  lspecfrm_acoag_in, lspectoo_acoag_in, &
                                  ip_aitacc_in, ip_aitpca_in, ip_pcaacc_in, &
                                  fac_m2v_aitage_in, fac_m2v_pcarbon_in, &
                                  nspec_max_in, ntot_amode_in, &
                                  modeptr_accum_in, modeptr_aitken_in, modeptr_pcarbon_in, &
                                  numptr_amode_in, mprognum_amode_in, nspec_amode_in, &
                                  lmassptr_amode_in, alnsg_amode_in, sigmag_amode_in, &
                                  r_universal_in, pstd_in, tmelt_in, boltz_in, &
                                  errmsg, errflg)
    !   store the resolved coagulation-pair tables, mode metadata, and host
    !   physical constants used by modal_aero_coag_run
    !   pair/species resolution and history-field registration are host
    !   responsibilities (see modal_aero_coag_cam)
    integer, intent(in) :: pair_option_acoag_in     ! pair selection (see module header)
    integer, intent(in) :: npair_acoag_in           ! number of coagulation pairs
    integer, intent(in) :: modefrm_acoag_in(:)      ! (maxpair_acoag) "from" mode of each pair
    integer, intent(in) :: modetoo_acoag_in(:)      ! (maxpair_acoag) "too" mode of each pair
    integer, intent(in) :: modetooeff_acoag_in(:)   ! (maxpair_acoag) effective "too" mode of each pair
    integer, intent(in) :: nspecfrm_acoag_in(:)     ! (maxpair_acoag) species count of each pair
    integer, intent(in) :: lspecfrm_acoag_in(:, :)   ! (nspec_max,maxpair_acoag) "from" species
    ! indices (host constituent space)
    integer, intent(in) :: lspectoo_acoag_in(:, :)   ! (nspec_max,maxpair_acoag) "too" species
    ! indices (host constituent space)
    integer, intent(in) :: ip_aitacc_in             ! pair index of [aitken-->accum]
    integer, intent(in) :: ip_aitpca_in             ! pair index of [aitken-->pcarbon]
    integer, intent(in) :: ip_pcaacc_in             ! pair index of [pcarbon-->accum]
    real(kind_phys), intent(in) :: fac_m2v_aitage_in(:)     ! (nspec_max) mixing-ratio to volume factors
    real(kind_phys), intent(in) :: fac_m2v_pcarbon_in(:)    ! (nspec_max) for aging shell/core calcs
    integer, intent(in) :: nspec_max_in             ! max number of species in a mode
    integer, intent(in) :: ntot_amode_in            ! number of aerosol modes
    integer, intent(in) :: modeptr_accum_in         ! accumulation mode index
    integer, intent(in) :: modeptr_aitken_in        ! aitken mode index
    integer, intent(in) :: modeptr_pcarbon_in       ! primary carbon mode index
    integer, intent(in) :: numptr_amode_in(:)       ! (ntot_amode) number indices (host constituent space)
    integer, intent(in) :: mprognum_amode_in(:)     ! (ntot_amode) prognostic-number flags
    integer, intent(in) :: nspec_amode_in(:)        ! (ntot_amode) species counts
    integer, intent(in) :: lmassptr_amode_in(:, :)   ! (nspec_max,ntot_amode) mass indices
    ! (host constituent space)
    real(kind_phys), intent(in) :: alnsg_amode_in(:)        ! (ntot_amode) ln(sigmag)
    real(kind_phys), intent(in) :: sigmag_amode_in(:)       ! (ntot_amode) geometric standard deviation
    real(kind_phys), intent(in) :: r_universal_in           ! universal gas constant (J/K/kmol)
    real(kind_phys), intent(in) :: pstd_in                  ! standard pressure (Pa)
    real(kind_phys), intent(in) :: tmelt_in                 ! freezing point of water (K)
    real(kind_phys), intent(in) :: boltz_in                 ! Boltzmann constant (J/K)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    pair_option_acoag = pair_option_acoag_in

    maxspec_acoag = nspec_max_in
    allocate (lspecfrm_acoag(maxspec_acoag, maxpair_acoag))
    allocate (lspectoo_acoag(maxspec_acoag, maxpair_acoag))
    allocate (fac_m2v_aitage(nspec_max_in), fac_m2v_pcarbon(nspec_max_in))

    npair_acoag = npair_acoag_in
    modefrm_acoag(:) = modefrm_acoag_in(:)
    modetoo_acoag(:) = modetoo_acoag_in(:)
    modetooeff_acoag(:) = modetooeff_acoag_in(:)
    nspecfrm_acoag(:) = nspecfrm_acoag_in(:)
    lspecfrm_acoag(:, :) = lspecfrm_acoag_in(:, :)
    lspectoo_acoag(:, :) = lspectoo_acoag_in(:, :)

    ip_aitacc = ip_aitacc_in
    ip_aitpca = ip_aitpca_in
    ip_pcaacc = ip_pcaacc_in

    fac_m2v_aitage(:) = fac_m2v_aitage_in(:)
    fac_m2v_pcarbon(:) = fac_m2v_pcarbon_in(:)

    ntot_amode = ntot_amode_in
    modeptr_accum = modeptr_accum_in
    modeptr_aitken = modeptr_aitken_in
    modeptr_pcarbon = modeptr_pcarbon_in

    allocate (numptr_amode(ntot_amode_in), mprognum_amode(ntot_amode_in), &
              nspec_amode(ntot_amode_in))
    allocate (lmassptr_amode(nspec_max_in, ntot_amode_in))
    allocate (alnsg_amode(ntot_amode_in), sigmag_amode(ntot_amode_in))
    numptr_amode(:) = numptr_amode_in(:)
    mprognum_amode(:) = mprognum_amode_in(:)
    nspec_amode(:) = nspec_amode_in(:)
    lmassptr_amode(:, :) = lmassptr_amode_in(:, :)
    alnsg_amode(:) = alnsg_amode_in(:)
    sigmag_amode(:) = sigmag_amode_in(:)

    r_universal = r_universal_in
    p0 = pstd_in
    tmelt = tmelt_in
    boltz = boltz_in

  end subroutine modal_aero_coag_init

  !   computes changes due to coagulation involving
  !       aitken mode (modeptr_aitken) with accumulation mode (modeptr_accum)
  !   this version will
  !       compute changes to mass and number, but not to surface area
  !       calculates coagulation rate coefficients using either
  !           new CMAQ V4.6 fast method
  !           older cmaq slow method (direct gauss-hermite quadrature)
  ! Authors: R. Easter
  ! RCE 07.04.15:  Adapted from MIRAGE2 code and CMAQ V4.6 code
  subroutine modal_aero_coag_run( &
    ncol, pver, top_lev, &
    loffset, nstep, &
    deltat_main, &
    t, pmid, &
    q, &
    dgncur_a, dgncur_awet, &
    wetdens_a, &
    dqdt, dotend, &
    errmsg, errflg)
    use modal_aero_gasaerexch, only: n_so4_monolayers_pcage

    integer, intent(in)  :: ncol
    integer, intent(in)  :: pver
    integer, intent(in)  :: top_lev          ! top level for modal aerosol calculations
    integer, intent(in)  :: loffset          ! offset applied to modal aero "pointers"
    integer, intent(in)  :: nstep            ! model step (for coagulation sub-cycling)

    real(kind_phys), intent(in) :: deltat_main      ! model timestep (s)

    real(kind_phys), intent(in) :: t(:, :)           ! (ncol,pver) temperature (K)
    real(kind_phys), intent(in) :: pmid(:, :)        ! (ncol,pver) pressure at model levels (Pa)

    real(kind_phys), intent(inout) :: q(:, :, :)      ! (ncol,pver,:)
    ! tracer mixing ratio (TMR) array
    ! *** MUST BE mol/mol-air or #/mol-air
    ! *** updated in place: number changes are
    !     direct assignments, and dqdt*deltat_main
    !     is NOT bit-identical to the stored change
    !     (deltatinv_main carries a 1+1e-15 guard),
    !     so this scheme cannot be tendency-return
    real(kind_phys), intent(in) :: dgncur_a(:, :, :)  ! (ncol,pver,ntot_amode)
    ! dry geo. mean dia. (m) of number distrib.
    real(kind_phys), intent(in) :: dgncur_awet(:, :, :)
    ! (ncol,pver,ntot_amode)
    ! wet geo. mean dia. (m) of number distrib.
    real(kind_phys), intent(in) :: wetdens_a(:, :, :) ! (ncol,pver,ntot_amode)
    ! density of wet aerosol (kg/m3)
    real(kind_phys), intent(out) :: dqdt(:, :, :)     ! (ncol,pver,:) TMR "dq/dt" array
    ! (diagnostic only; q is updated in place)
    logical, intent(out) :: dotend(:)       ! identifies the species that
    ! tendencies are computed for
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: i, ipair, iq
    integer :: idomode(ntot_amode), iselfcoagdone(ntot_amode)
    integer :: jfreqcoag
    integer :: k
    integer :: l, lsfrm, lstoo, lunout
    integer :: modefrm, modetoo, mait, macc, mpca
    integer ::  n, nfreqcoag

    integer, save :: nerr = 0       ! number of errors for entire run
    integer, save :: nerrmax = 9999 ! maximum number of errors before abort
    integer, parameter :: ldiag1 = -1, ldiag2 = -1, ldiag3 = -1

    logical, parameter :: fastcoag_flag = .true. ! selects coag rate-coef method

    real(kind_phys) :: aircon
    real(kind_phys) :: deltat, deltatinv_main
    real(kind_phys) :: dr_so4_monolayers_pcage
    real(kind_phys) :: dumexp, dumloss, dumprod
    real(kind_phys) :: dumsfc_frm_old, dumsfc_frm_new
    real(kind_phys) :: dum_m2v
    real(kind_phys) :: fac_volsfc_pcarbon
    real(kind_phys) :: lnsg_frm, lnsg_too
    real(kind_phys) :: sg_frm, sg_too
    real(kind_phys) :: tmpa, tmpb, tmpc, tmpf, tmpg, tmph, tmpn
    real(kind_phys) :: tmp1, tmp2
    real(kind_phys) :: vol_core, vol_shell
    real(kind_phys) :: xbetaij0, xbetaij2i, xbetaij2j, xbetaij3, &
                       xbetaii0, xbetaii2, xbetajj0, xbetajj2
    real(kind_phys) :: xferamt, xferfracvol, xferfrac_pcage, xferfrac_max
    real(kind_phys) :: xnumbconc(ntot_amode)
    real(kind_phys) :: xnumbconcavg(ntot_amode), xnumbconcnew(ntot_amode)
    real(kind_phys) :: ybetaij0(maxpair_acoag), ybetaij3(maxpair_acoag)
    real(kind_phys) :: ybetaii0(maxpair_acoag), ybetajj0(maxpair_acoag)

    ! begin
    errmsg = ' '
    errflg = 0

    !   zero the tendency outputs up front: they are intent(out) and the caller
    !   uses them unconditionally, including on the bypass paths below
    dotend(:) = .false.
    dqdt(1:ncol, :, :) = 0.0_kind_phys

    !   check if any coagulation pairs exist
    if (npair_acoag <= 0) return

    lunout = 6
    !
    !   determine if coagulation will be done on this time-step
    !   currently coagulation is done every 3 hours
    !
    !       deltat = 3600.0*3.0
    deltat = deltat_main
    nfreqcoag = max(1, nint(deltat/deltat_main))
    jfreqcoag = nfreqcoag/2
    xferfrac_max = 1.0_kind_phys - 10.0_kind_phys*epsilon(1.0_kind_phys)   ! 1-eps

    if (nfreqcoag > 1) then
      if (mod(nstep, nfreqcoag) /= jfreqcoag) return
    end if

!
!   set idomode
!
    idomode(:) = 0
    do ipair = 1, npair_acoag
      idomode(modefrm_acoag(ipair)) = 1
      idomode(modetoo_acoag(ipair)) = 1
    end do

!
!   other init
!
    macc = modeptr_accum
    mait = modeptr_aitken
    mpca = modeptr_pcarbon

    if (mpca > 0 .and. mpca <= ntot_amode) then
      ! use 1 mol (bi-)sulfate = 65 cm^3 --> 1 molecule = (4.76e-10 m)^3
      dr_so4_monolayers_pcage = n_so4_monolayers_pcage*4.76e-10_kind_phys
      fac_volsfc_pcarbon = exp(2.5_kind_phys*(alnsg_amode(mpca)**2))
    end if

!
!   loop over levels and columns to calc the coagulation
!
!   integrate coagulation changes over deltat = nfreqcoag*deltat_main
!   then compute tendencies as
!      dqdt = (q(t+deltat) - q(t))/deltat_main
!   because tendencies are applied (in physics_update) over deltat_main
!
    deltat = nfreqcoag*deltat_main
    deltatinv_main = 1.0_kind_phys/(deltat_main*(1.0_kind_phys + 1.0e-15_kind_phys))

    main_k: do k = top_lev, pver
    main_i: do i = 1, ncol

!   air molar density (kmol/m3)
      aircon = (pmid(i, k)/(r_universal*t(i, k)))

!   calculate number conc. (#/m3) for modes doing coagulation
      do n = 1, ntot_amode
        if (idomode(n) > 0) then
          xnumbconc(n) = q(i, k, numptr_amode(n) - loffset)*aircon
          xnumbconc(n) = max(0.0_kind_phys, xnumbconc(n))
        end if
        iselfcoagdone(n) = 0
      end do

      ! calculate coagulation rates for each pair
      main_ipair1: do ipair = 1, npair_acoag

        modefrm = modefrm_acoag(ipair)
        modetoo = modetoo_acoag(ipair)

        !
        ! compute coagulation rates using cmaq "fast" method
        !    (based on E. Whitby's approximation approach)
        ! here subr. arguments are all in mks unit
        !
        call getcoags_wrapper_f( &
          t(i, k), pmid(i, k), &
          dgncur_awet(i, k, modefrm), dgncur_awet(i, k, modetoo), &
          sigmag_amode(modefrm), sigmag_amode(modetoo), &
          alnsg_amode(modefrm), alnsg_amode(modetoo), &
          wetdens_a(i, k, modefrm), wetdens_a(i, k, modetoo), &
          xbetaij0, xbetaij2i, xbetaij2j, xbetaij3, &
          xbetaii0, xbetaii2, xbetajj0, xbetajj2)

        ybetaij0(ipair) = xbetaij0
        ybetaij3(ipair) = xbetaij3
        ybetaii0(ipair) = xbetaii0
        ybetajj0(ipair) = xbetajj0

      end do main_ipair1

      if ((pair_option_acoag == 1) .or. &
          (pair_option_acoag == 2)) then
        ! calculate number and mass changes for pair_option_acoag == 1,2
        main_ipair2: do ipair = 1, npair_acoag

          modefrm = modefrm_acoag(ipair)
          modetoo = modetoo_acoag(ipair)

          !   calculate number changes
          !   apply self-coagulation losses only once to a mode (when iselfcoagdone=0)
          !       first calc change to "too" mode
          !       next  calc change to "frm" mode, using average number conc of "too"
          if ((mprognum_amode(modetoo) > 0) .and. &
              (iselfcoagdone(modetoo) <= 0)) then
            iselfcoagdone(modetoo) = 1
            tmpn = xnumbconc(modetoo)
            xnumbconcnew(modetoo) = tmpn/(1.0_kind_phys + deltat*ybetajj0(ipair)*tmpn)
            xnumbconcavg(modetoo) = 0.5_kind_phys*(xnumbconcnew(modetoo) + tmpn)
            lstoo = numptr_amode(modetoo) - loffset
            q(i, k, lstoo) = xnumbconcnew(modetoo)/aircon
            dqdt(i, k, lstoo) = (xnumbconcnew(modetoo) - tmpn)*deltatinv_main/aircon
          end if

          if ((mprognum_amode(modefrm) > 0) .and. &
              (iselfcoagdone(modefrm) <= 0)) then
            iselfcoagdone(modefrm) = 1
            tmpn = xnumbconc(modefrm)
            tmpa = deltat*ybetaij0(ipair)*xnumbconcavg(modetoo)
            tmpb = deltat*ybetaii0(ipair)
            tmpc = tmpa + tmpb*tmpn
            if (abs(tmpc) < 0.01_kind_phys) then
              xnumbconcnew(modefrm) = tmpn*exp(-tmpc)
            else if (abs(tmpa) < 0.001_kind_phys) then
              xnumbconcnew(modefrm) = &
                exp(-tmpa)*tmpn/(1.0_kind_phys + tmpb*tmpn)
            else
              tmpf = tmpb*tmpn/tmpc
              tmpg = exp(-tmpa)
              tmph = tmpg*(1.0_kind_phys - tmpf)/(1.0_kind_phys - tmpg*tmpf)
              xnumbconcnew(modefrm) = tmpn*max(0.0_kind_phys, min(1.0_kind_phys, tmph))
            end if
            xnumbconcavg(modefrm) = 0.5_kind_phys*(xnumbconcnew(modefrm) + tmpn)
            lsfrm = numptr_amode(modefrm) - loffset
            q(i, k, lsfrm) = xnumbconcnew(modefrm)/aircon
            dqdt(i, k, lsfrm) = (xnumbconcnew(modefrm) - tmpn)*deltatinv_main/aircon
          end if

          !   calculate mass changes
          !     xbetaij3*xnumbconc(modetoo) = first order loss rate for modefrm volume
          !     xferfracvol = fraction of modefrm volume transferred to modetoo over deltat
          dumloss = ybetaij3(ipair)*xnumbconcavg(modetoo)
          xferfracvol = 1.0_kind_phys - exp(-dumloss*deltat)
          xferfracvol = max(0.0_kind_phys, min(xferfrac_max, xferfracvol))

          do iq = 1, nspecfrm_acoag(ipair)
            lsfrm = lspecfrm_acoag(iq, ipair) - loffset
            lstoo = lspectoo_acoag(iq, ipair) - loffset
            if (lsfrm > 0) then
              xferamt = q(i, k, lsfrm)*xferfracvol
              dqdt(i, k, lsfrm) = dqdt(i, k, lsfrm) - xferamt*deltatinv_main
              q(i, k, lsfrm) = q(i, k, lsfrm) - xferamt
              if (lstoo > 0) then
                dqdt(i, k, lstoo) = dqdt(i, k, lstoo) + xferamt*deltatinv_main
                q(i, k, lstoo) = q(i, k, lstoo) + xferamt
              end if
            end if
          end do

        end do main_ipair2

      else if (pair_option_acoag == 3) then
      !
      !   calculate number and mass changes for pair_option_acoag == 3
      !

        ! calculate number changes to accum mode
        if (mprognum_amode(macc) > 0) then
          tmpn = xnumbconc(macc)
          xnumbconcnew(macc) = tmpn/(1.0_kind_phys + deltat*ybetajj0(ip_aitacc)*tmpn)
          xnumbconcavg(macc) = 0.5_kind_phys*(xnumbconcnew(macc) + tmpn)
          lstoo = numptr_amode(macc) - loffset
          q(i, k, lstoo) = xnumbconcnew(macc)/aircon
          dqdt(i, k, lstoo) = (xnumbconcnew(macc) - tmpn)*deltatinv_main/aircon
        end if

        ! calculate number changes to primary carbon mode
        modefrm = modeptr_pcarbon
        if (mprognum_amode(mpca) > 0) then
          tmpn = xnumbconc(mpca)
          tmpa = deltat*ybetaij0(ip_pcaacc)*xnumbconcavg(macc)
          tmpb = deltat*ybetaii0(ip_pcaacc)
          tmpc = tmpa + tmpb*tmpn
          if (abs(tmpc) < 0.01_kind_phys) then
            xnumbconcnew(mpca) = tmpn*exp(-tmpc)
          else if (abs(tmpa) < 0.001_kind_phys) then
            xnumbconcnew(mpca) = &
              exp(-tmpa)*tmpn/(1.0_kind_phys + tmpb*tmpn)
          else
            tmpf = tmpb*tmpn/tmpc
            tmpg = exp(-tmpa)
            tmph = tmpg*(1.0_kind_phys - tmpf)/(1.0_kind_phys - tmpg*tmpf)
            xnumbconcnew(mpca) = tmpn*max(0.0_kind_phys, min(1.0_kind_phys, tmph))
          end if
          xnumbconcavg(mpca) = 0.5_kind_phys*(xnumbconcnew(mpca) + tmpn)
          lsfrm = numptr_amode(mpca) - loffset
          q(i, k, lsfrm) = xnumbconcnew(mpca)/aircon
          dqdt(i, k, lsfrm) = (xnumbconcnew(mpca) - tmpn)*deltatinv_main/aircon
        end if

        ! calculate number changes to aitken mode
        if (mprognum_amode(mait) > 0) then
          tmpn = xnumbconc(mait)
          tmpa = deltat*(ybetaij0(ip_aitacc)*xnumbconcavg(macc) &
                         + ybetaij0(ip_aitpca)*xnumbconcavg(mpca))
          tmpb = deltat*ybetaii0(ip_aitacc)
          tmpc = tmpa + tmpb*tmpn
          if (abs(tmpc) < 0.01_kind_phys) then
            xnumbconcnew(mait) = tmpn*exp(-tmpc)
          else if (abs(tmpa) < 0.001_kind_phys) then
            xnumbconcnew(mait) = &
              exp(-tmpa)*tmpn/(1.0_kind_phys + tmpb*tmpn)
          else
            tmpf = tmpb*tmpn/tmpc
            tmpg = exp(-tmpa)
            tmph = tmpg*(1.0_kind_phys - tmpf)/(1.0_kind_phys - tmpg*tmpf)
            xnumbconcnew(mait) = tmpn*max(0.0_kind_phys, min(1.0_kind_phys, tmph))
          end if
          xnumbconcavg(mait) = 0.5_kind_phys*(xnumbconcnew(mait) + tmpn)
          lsfrm = numptr_amode(mait) - loffset
          q(i, k, lsfrm) = xnumbconcnew(mait)/aircon
          dqdt(i, k, lsfrm) = (xnumbconcnew(mait) - tmpn)*deltatinv_main/aircon
        end if

        !   calculate mass changes from aitken-->accum direct coagulation and
        !       aitken-->pcarbon-->accum coagulation/aging
        !   also calc volume of shell material (so4 & nh4 from aitken-->pcarbon)
        dumloss = ybetaij3(ip_aitacc)*xnumbconcavg(macc) &
                  + ybetaij3(ip_aitpca)*xnumbconcavg(mpca)
        tmpa = ybetaij3(ip_aitpca)*xnumbconcavg(mpca)/max(dumloss, 1.0e-37_kind_phys)
        xferfracvol = 1.0_kind_phys - exp(-dumloss*deltat)
        xferfracvol = max(0.0_kind_phys, min(xferfrac_max, xferfracvol))
        vol_shell = 0.0_kind_phys

        ipair = ip_aitacc
        do iq = 1, nspecfrm_acoag(ipair)
          lsfrm = lspecfrm_acoag(iq, ipair) - loffset
          lstoo = lspectoo_acoag(iq, ipair) - loffset
          if (lsfrm > 0) then
            xferamt = q(i, k, lsfrm)*xferfracvol
            dqdt(i, k, lsfrm) = dqdt(i, k, lsfrm) - xferamt*deltatinv_main
            q(i, k, lsfrm) = q(i, k, lsfrm) - xferamt
            if (lstoo > 0) then
              dqdt(i, k, lstoo) = dqdt(i, k, lstoo) + xferamt*deltatinv_main
              q(i, k, lstoo) = q(i, k, lstoo) + xferamt
            end if
            vol_shell = vol_shell + xferamt*tmpa*fac_m2v_aitage(iq)
          end if
        end do

        !   now calculate aging transfer fraction for pcarbon-->accum
        !   this duplicates the code in modal_aero_gasaerexch
        vol_core = 0.0_kind_phys
        do l = 1, nspec_amode(mpca)
          vol_core = vol_core + &
                     q(i, k, lmassptr_amode(l, mpca) - loffset)*fac_m2v_pcarbon(l)
        end do
        tmp1 = vol_shell*dgncur_a(i, k, mpca)*fac_volsfc_pcarbon
        tmp2 = 6.0_kind_phys*dr_so4_monolayers_pcage*vol_core
        tmp2 = max(tmp2, 0.0_kind_phys)
        if (tmp1 >= tmp2) then
          xferfrac_pcage = xferfrac_max
        else
          xferfrac_pcage = min(tmp1/tmp2, xferfrac_max)
        end if

        !   calculate mass changes from pcarbon-->accum by direct coagulation
        !   and aging
        dumloss = ybetaij3(ip_pcaacc)*xnumbconcavg(macc)
        xferfracvol = 1.0_kind_phys - exp(-dumloss*deltat)
        xferfracvol = xferfracvol + xferfrac_pcage
        xferfracvol = max(0.0_kind_phys, min(xferfrac_max, xferfracvol))

        ipair = ip_pcaacc
        do iq = 1, nspecfrm_acoag(ipair)
          lsfrm = lspecfrm_acoag(iq, ipair) - loffset
          lstoo = lspectoo_acoag(iq, ipair) - loffset
          if (lsfrm > 0) then
            xferamt = q(i, k, lsfrm)*xferfracvol
            dqdt(i, k, lsfrm) = dqdt(i, k, lsfrm) - xferamt*deltatinv_main
            q(i, k, lsfrm) = q(i, k, lsfrm) - xferamt
            if (lstoo > 0) then
              dqdt(i, k, lstoo) = dqdt(i, k, lstoo) + xferamt*deltatinv_main
              q(i, k, lstoo) = q(i, k, lstoo) + xferamt
            end if
          end if
        end do

        lsfrm = numptr_amode(mpca) - loffset
        lstoo = numptr_amode(macc) - loffset
        if (lsfrm > 0) then
          xferamt = q(i, k, lsfrm)*xferfrac_pcage
          dqdt(i, k, lsfrm) = dqdt(i, k, lsfrm) - xferamt*deltatinv_main
          q(i, k, lsfrm) = q(i, k, lsfrm) - xferamt
          if (lstoo > 0) then
            dqdt(i, k, lstoo) = dqdt(i, k, lstoo) + xferamt*deltatinv_main
            q(i, k, lstoo) = q(i, k, lstoo) + xferamt
          end if
        end if

      else   ! (pair_option_acoag /= 1,2,3) then

        write (lunout, *) '*** modal_aero_coag_sub error'
        write (lunout, *) '    cannot do _coag_sub error pair_option_acoag =', &
          pair_option_acoag
        errmsg = 'modal_aero_coag_sub error'
        errflg = 1
        return

      end if   ! (pair_option_acoag == ...)

    end do main_i
    end do main_k

    ! set dotend's
    do ipair = 1, npair_acoag
      modefrm = modefrm_acoag(ipair)
      modetoo = modetoo_acoag(ipair)

      do iq = 1, nspecfrm_acoag(ipair)
        lsfrm = lspecfrm_acoag(iq, ipair) - loffset
        lstoo = lspectoo_acoag(iq, ipair) - loffset
        if (lsfrm > 0) dotend(lsfrm) = .true.
        if (lstoo > 0) dotend(lstoo) = .true.
      end do

      if (mprognum_amode(modefrm) > 0) then
        lsfrm = numptr_amode(modefrm) - loffset
        if (lsfrm > 0) dotend(lsfrm) = .true.
      end if
      if (mprognum_amode(modetoo) > 0) then
        lstoo = numptr_amode(modetoo) - loffset
        if (lstoo > 0) dotend(lstoo) = .true.
      end if

    end do

  end subroutine modal_aero_coag_run

  subroutine getcoags_wrapper_f( &
    airtemp, airprs, &
    dgatk, dgacc, &
    sgatk, sgacc, &
    xxlsgat, xxlsgac, &
    pdensat, pdensac, &
    betaij0, betaij2i, betaij2j, betaij3, &
    betaii0, betaii2, betajj0, betajj2)
!   (p0, tmelt, boltz are module-level host constants
!    set by modal_aero_coag_init)
!
! interface to subr. getcoags
!
! interface code adapted from subr. aeroproc of cmaq v4.6,
!     with some of the parameter values from module aero_info_ae4
!

! *** arguments

    real(kind_phys), intent(in) :: airtemp  ! air temperature [ k ]
    real(kind_phys), intent(in) :: airprs   ! air pressure in [ pa ]

    real(kind_phys), intent(in) :: dgatk    ! aitken mode geometric mean diameter [m]
    real(kind_phys), intent(in) :: dgacc    ! accumulation mode geometric mean diam [m]

    real(kind_phys), intent(in) :: sgatk    ! aitken mode geometric standard deviation
    real(kind_phys), intent(in) :: sgacc    ! accumulation mode geometric standard deviation

    real(kind_phys), intent(in) :: xxlsgat  ! natural log of geometric standard
    real(kind_phys), intent(in) :: xxlsgac  !  deviations

    real(kind_phys), intent(in) :: pdensat  ! aitken mode particle density [ kg / m**3 ]
    real(kind_phys), intent(in) :: pdensac  ! accumulation mode density [ kg / m**3 ]

    real(kind_phys), intent(out) :: betaij0, betaij2i, betaij2j, betaij3, &
                                    betaii0, betaii2, betajj0, betajj2

! *** local parameters
    real(kind_phys) :: t0  ! standard surface temperature (15 deg C) [ k ]
    real(kind_phys), parameter :: two3 = 2.0_kind_phys/3.0_kind_phys

! *** local variables
    real(kind_phys) :: amu            ! atmospheric dynamic viscosity [ kg/m s ]
    real(kind_phys) :: sqrt_temp      ! square root of ambient temperature
    real(kind_phys) :: lamda          ! mean free path [ m ]

! *** intramodal coagulation rates [ m**3/s ] ( 0th & 2nd moments )
    real(kind_phys)    :: batat(2)  ! aitken mode
    real(kind_phys)    :: bacac(2)  ! accumulation mode
! *** intermodal coagulation rates [ m**3/s ] ( 0th & 2nd moments )
    real(kind_phys)    :: batac(2)  ! aitken to accumulation
    real(kind_phys)    :: bacat(2)  ! accumulation from aitken
! *** intermodal coagulation rate [ m**3/s ] ( 3rd moment )
    real(kind_phys)    :: c3ij        ! aitken to accumulation
! *** 3rd moment intermodal transfer rate by coagulation
    real(kind_phys)    :: c30atac     ! aitken to accumulation

! *** near continnuum regime (independent of mode)
    real(kind_phys)    :: knc         ! knc = two3 * boltz *  airtemp / amu
! *** free molecular regime (depends upon modal density)
    real(kind_phys)    :: kfmat       ! kfmat = sqrt(3.0*boltz*airtemp/pdensat)
    real(kind_phys)    :: kfmac       ! kfmac = sqrt(3.0*boltz*airtemp/pdensac)
    real(kind_phys)    :: kfmatac     ! kfmatac = sqrt( 6.0 * boltz * airtemp /
    !                ( pdensat + pdensac ) )

    real(kind_phys)    :: dumacc2, dumatk2, dumatk3

    t0 = tmelt + 15._kind_phys

    sqrt_temp = sqrt(airtemp)

! *** calculate mean free path [ m ]:
!     6.6328e-8 is the sea level value given in table i.2.8
!     on page 10 of u.s. standard atmosphere 1962
    lamda = 6.6328e-8_kind_phys*p0*airtemp/(t0*airprs)

! *** calculate dynamic viscosity [ kg m**-1 s**-1 ]:
!     u.s. standard atmosphere 1962 page 14 expression
!     for dynamic viscosity is:
!     dynamic viscosity =  beta * t * sqrt(t) / ( t + s)
!     where beta = 1.458e-6 [ kg sec^-1 k**-0.5 ], s = 110.4 [ k ].
    amu = 1.458e-6_kind_phys*airtemp*sqrt_temp/(airtemp + 110.4_kind_phys)

! *** coagulation
!     calculate coagulation coefficients using a method dictated by
!     the value of fastcoag_flag.  if true, the computationally-
!     efficient getcoags routine is used.  if false, the more intensive
!     gauss-hermite numerical quadrature method is used.  see section
!     2.1 of bhave et al. (2004) for further discussion.

! *** calculate term used in equation a6 of binkowski & shankar (1995)
    knc = two3*boltz*airtemp/amu
! *** calculate terms used in equation a5 of binkowski & shankar (1995)
    kfmat = sqrt(3.0_kind_phys*boltz*airtemp/pdensat)
    kfmac = sqrt(3.0_kind_phys*boltz*airtemp/pdensac)
    kfmatac = sqrt(6.0_kind_phys*boltz*airtemp/(pdensat + pdensac))

! *** transfer of number to accumulation mode from aitken mode is zero
    bacat(1) = 0.0_kind_phys

! *** calculate intermodal and intramodal coagulation coefficients
!     for zeroth and second moments, and intermodal coagulation
!     coefficient for third moment
    call getcoags(lamda, kfmatac, kfmat, kfmac, knc, &
                  dgatk, dgacc, sgatk, sgacc, &
                  xxlsgat, xxlsgac, &
                  batat(2), batat(1), bacac(2), bacac(1), &
                  batac(2), bacat(2), batac(1), c3ij)

! convert from the "cmaq" coag rate parameters
! to the "mirage2" parameters
    dumacc2 = ((dgacc**2)*exp(2.0_kind_phys*xxlsgac*xxlsgac))
    dumatk2 = ((dgatk**2)*exp(2.0_kind_phys*xxlsgat*xxlsgat))
    dumatk3 = ((dgatk**3)*exp(4.5_kind_phys*xxlsgat*xxlsgat))

    betaii0 = max(0.0_kind_phys, batat(1))
    betajj0 = max(0.0_kind_phys, bacac(1))
    betaij0 = max(0.0_kind_phys, batac(1))
    betaij3 = max(0.0_kind_phys, c3ij/dumatk3)

    betajj2 = max(0.0_kind_phys, bacac(2)/dumacc2)
    betaii2 = max(0.0_kind_phys, batat(2)/dumatk2)
    betaij2i = max(0.0_kind_phys, batac(2)/dumatk2)
    betaij2j = max(0.0_kind_phys, bacat(2)/dumatk2)

  end subroutine getcoags_wrapper_f

! //////////////////////////////////////////////////////////////////
!  subroutine getcoags calculates the coagulation rates using a new
!     approximate algorithm for the 2nd moment.  the 0th and 3rd moments
!     are done by analytic expressions from whitby et al. (1991).  the
!     correction factors are also similar to those from whitby et al.
!     (1991), but are derived from the gauss-hermite numerical
!     quadratures used by binkowski and roselle (2003).
!
!     called from aerostep as:
!     call getcoags( lamda, kfmatac, kfmat, kfmac, knc,
!                    dgat,dgac, sgatk, sgacc, xxlsgat,xxlsgac,
!                    batat(2), batat(1), bacac(2), bacac(1),
!                    batac(2), bacat(2), batac(1), c3ij )
!     where all input and outputs are real*8
!
!  revision history:
!   fsb 08/25/03 coded by dr. francis s. binkowksi
!
!   fsb 08/25/04 added in-line documentation
!
!   rce 04/15/2007
!       code taken from cmaq v4.6 code; converted to f90;
!       added "intent" to subr arguments;
!       renamed "r4" & "kind_phys" variables to "rx4" & "rx8";
!       changed "real*N" declarations to "real(rN)" (N = 4 or 8)
!
!  references:
!   1. whitby, e. r., p. h. mcmurry, u. shankar, and f. s. binkowski,
!   modal aerosol dynamics modeling, rep. 600/3-91/020, atmospheric
!   research and exposure assessment laboratory,
!   u.s. environmental protection agency, research triangle park, n.c.,
!   (ntis pb91-161729/as), 1991
!
!   2. binkowski, f.s. an u. shankar, the regional particulate matter
!   model 1. model decsription and preliminary results, journal of
!   geophysical research, 100, d12, pp 26,191-26,209,
!   december 20, 1995.
!
!   3. binkowski, f.s. and s.j. roselle, models-3 community
!      multiscale air quality (cmaq) model aerosol component 1:
!      model description.  j. geophys. res., vol 108, no d6, 4183
!      doi:10.1029/2001jd001409, 2003.

  subroutine getcoags(lamda, kfmatac, kfmat, kfmac, knc, &
                      dgatk, dgacc, sgatk, sgacc, xxlsgat, xxlsgac, &
                      qs11, qn11, qs22, qn22, &
                      qs12, qs21, qn12, qv12)

    real(kind_phys), intent(in) ::  lamda     ! mean free path [ m ]

! *** coefficients for free molecular regime
    real(kind_phys), intent(in) ::  kfmat     ! aitken mode
    real(kind_phys), intent(in) ::  kfmac     ! accumulation mode
    real(kind_phys), intent(in) ::  kfmatac   ! aitken to accumulation mode

    real(kind_phys), intent(in) ::  knc   ! coefficient for near continnuum regime

! *** modal geometric mean diameters: [ m ]
    real(kind_phys), intent(in) :: dgatk          ! aitken mode
    real(kind_phys), intent(in) :: dgacc          ! accumulation mode

! *** modal geometric standard deviation
    real(kind_phys), intent(in) :: sgatk          ! atken mode
    real(kind_phys), intent(in) :: sgacc          ! accumulation mode

! *** natural log of modal geometric standard deviation
    real(kind_phys), intent(in) :: xxlsgat         ! aitken mode
    real(kind_phys), intent(in) :: xxlsgac         ! accumulation mode

! *** coagulation coefficients
    real(kind_phys), intent(out) :: qs11, qn11, qs22, qn22, &
                                    qs12, qs21, qn12, qv12

    integer :: ibeta, n1, n2a, n2n ! indices for correction factors

    real(kind_phys)  :: i1fm_at
    real(kind_phys)  :: i1nc_at
    real(kind_phys)  :: i1_at

    real(kind_phys)  :: i1fm_ac
    real(kind_phys)  :: i1nc_ac
    real(kind_phys)  :: i1_ac

    real(kind_phys)  :: i1fm
    real(kind_phys)  :: i1nc
    real(kind_phys)  :: i1

    real(kind_phys) :: constii

    real(kind_phys)    :: kngat, kngac
    real(kind_phys)    :: one, two, half
    parameter(one=1.0_kind_phys, two=2.0_kind_phys, half=0.5_kind_phys)
    real(kind_phys)    :: a
!       parameter( a = 2.492_kind_phys)
    parameter(a=1.246_kind_phys)
    real(kind_phys)      :: two3rds
    parameter(two3rds=2._kind_phys/3._kind_phys)

    real(kind_phys)   :: sqrttwo  !  sqrt(two)
    real(kind_phys)   :: dlgsqt2  !  1/ln( sqrt( 2 ) )

    real(kind_phys)    :: esat01         ! aitken mode exp( log^2( sigmag )/8 )
    real(kind_phys)    :: esac01         ! accumulation mode exp( log^2( sigmag )/8 )

    real(kind_phys)    :: esat04
    real(kind_phys)    :: esac04

    real(kind_phys)    :: esat05
    real(kind_phys)    :: esac05

    real(kind_phys)    :: esat08
    real(kind_phys)    :: esac08

    real(kind_phys)    :: esat09
    real(kind_phys)    :: esac09

    real(kind_phys)    :: esat16
    real(kind_phys)    :: esac16

    real(kind_phys)    :: esat20
    real(kind_phys)    :: esac20

    real(kind_phys)    :: esat24
    real(kind_phys)    :: esac24

    real(kind_phys)    :: esat25
    real(kind_phys)    :: esac25

    real(kind_phys)    :: esat36
    real(kind_phys)    :: esac36

    real(kind_phys)    :: esat49

    real(kind_phys)    :: esat64
    real(kind_phys)    :: esac64

    real(kind_phys)    :: esat100

    real(kind_phys) :: dgat2, dgac2, dgat3, dgac3
    real(kind_phys) :: sqdgat, sqdgac
    real(kind_phys) :: sqdgat5, sqdgac5
    real(kind_phys) :: sqdgat7
    real(kind_phys) :: r, r2, r3, rx4, r5, r6, rx8
    real(kind_phys) :: ri1, ri2, ri3, ri4
    real(kind_phys) :: rat
    real(kind_phys) :: coagfm0, coagnc0
    real(kind_phys) :: coagfm3, coagnc3
    real(kind_phys) :: coagfm_at, coagfm_ac
    real(kind_phys) :: coagnc_at, coagnc_ac
    real(kind_phys) :: coagatat0
    real(kind_phys) :: coagacac0
    real(kind_phys) :: coagatat2
    real(kind_phys) :: coagacac2
    real(kind_phys) :: coagatac0, coagatac3
    real(kind_phys) :: coagatac2
    real(kind_phys) :: coagacat2
    real(kind_phys) :: xm2at, xm3at, xm2ac, xm3ac

! *** correction factors for coagulation rates
    real(kind_phys), save :: bm0(10)          ! m0 intramodal fm - rpm values
    real(kind_phys), save :: bm0ij(10, 10, 10) ! m0 intermodal fm
    real(kind_phys), save :: bm3i(10, 10, 10) ! m3 intermodal fm- rpm values
    real(kind_phys), save :: bm2ii(10) ! m2 intramodal fm
    real(kind_phys), save :: bm2iitt(10) ! m2 intramodal total
    real(kind_phys), save :: bm2ij(10, 10, 10) ! m2 intermodal fm i to j
    real(kind_phys), save :: bm2ji(10, 10, 10) ! m2 total intermodal  j from i

! *** populate the arrays for the correction factors.

! rpm 0th moment correction factors for unimodal fm coagulation  rates
    data bm0/ &
      0.707106785165097_kind_phys, 0.726148960080488_kind_phys, 0.766430744110958_kind_phys, &
      0.814106389441342_kind_phys, 0.861679526483207_kind_phys, 0.903600509090092_kind_phys, &
      0.936578814219156_kind_phys, 0.960098926735545_kind_phys, 0.975646823342881_kind_phys, &
      0.985397173215326_kind_phys/

! fsb new fm correction factors for m0 intermodal coagulation

    data(bm0ij(1, 1, ibeta), ibeta=1, 10)/ &
      0.628539_kind_phys, 0.639610_kind_phys, 0.664514_kind_phys, 0.696278_kind_phys, 0.731558_kind_phys, &
      0.768211_kind_phys, 0.804480_kind_phys, 0.838830_kind_phys, 0.870024_kind_phys, 0.897248_kind_phys/
    data(bm0ij(1, 2, ibeta), ibeta=1, 10)/ &
      0.639178_kind_phys, 0.649966_kind_phys, 0.674432_kind_phys, 0.705794_kind_phys, 0.740642_kind_phys, &
      0.776751_kind_phys, 0.812323_kind_phys, 0.845827_kind_phys, 0.876076_kind_phys, 0.902324_kind_phys/
    data(bm0ij(1, 3, ibeta), ibeta=1, 10)/ &
      0.663109_kind_phys, 0.673464_kind_phys, 0.697147_kind_phys, 0.727637_kind_phys, 0.761425_kind_phys, &
      0.796155_kind_phys, 0.829978_kind_phys, 0.861419_kind_phys, 0.889424_kind_phys, 0.913417_kind_phys/
    data(bm0ij(1, 4, ibeta), ibeta=1, 10)/ &
      0.693693_kind_phys, 0.703654_kind_phys, 0.726478_kind_phys, 0.755786_kind_phys, 0.787980_kind_phys, &
      0.820626_kind_phys, 0.851898_kind_phys, 0.880459_kind_phys, 0.905465_kind_phys, 0.926552_kind_phys/
    data(bm0ij(1, 5, ibeta), ibeta=1, 10)/ &
      0.727803_kind_phys, 0.737349_kind_phys, 0.759140_kind_phys, 0.786870_kind_phys, 0.816901_kind_phys, &
      0.846813_kind_phys, 0.874906_kind_phys, 0.900060_kind_phys, 0.921679_kind_phys, 0.939614_kind_phys/
    data(bm0ij(1, 6, ibeta), ibeta=1, 10)/ &
      0.763461_kind_phys, 0.772483_kind_phys, 0.792930_kind_phys, 0.818599_kind_phys, 0.845905_kind_phys, &
      0.872550_kind_phys, 0.897051_kind_phys, 0.918552_kind_phys, 0.936701_kind_phys, 0.951528_kind_phys/
    data(bm0ij(1, 7, ibeta), ibeta=1, 10)/ &
      0.799021_kind_phys, 0.807365_kind_phys, 0.826094_kind_phys, 0.849230_kind_phys, 0.873358_kind_phys, &
      0.896406_kind_phys, 0.917161_kind_phys, 0.935031_kind_phys, 0.949868_kind_phys, 0.961828_kind_phys/
    data(bm0ij(1, 8, ibeta), ibeta=1, 10)/ &
      0.833004_kind_phys, 0.840514_kind_phys, 0.857192_kind_phys, 0.877446_kind_phys, 0.898147_kind_phys, &
      0.917518_kind_phys, 0.934627_kind_phys, 0.949106_kind_phys, 0.960958_kind_phys, 0.970403_kind_phys/
    data(bm0ij(1, 9, ibeta), ibeta=1, 10)/ &
      0.864172_kind_phys, 0.870734_kind_phys, 0.885153_kind_phys, 0.902373_kind_phys, 0.919640_kind_phys, &
      0.935494_kind_phys, 0.949257_kind_phys, 0.960733_kind_phys, 0.970016_kind_phys, 0.977346_kind_phys/
    data(bm0ij(1, 10, ibeta), ibeta=1, 10)/ &
      0.891658_kind_phys, 0.897227_kind_phys, 0.909343_kind_phys, 0.923588_kind_phys, 0.937629_kind_phys, &
      0.950307_kind_phys, 0.961151_kind_phys, 0.970082_kind_phys, 0.977236_kind_phys, 0.982844_kind_phys/
    data(bm0ij(2, 1, ibeta), ibeta=1, 10)/ &
      0.658724_kind_phys, 0.670587_kind_phys, 0.697539_kind_phys, 0.731890_kind_phys, 0.769467_kind_phys, &
      0.807391_kind_phys, 0.843410_kind_phys, 0.875847_kind_phys, 0.903700_kind_phys, 0.926645_kind_phys/
    data(bm0ij(2, 2, ibeta), ibeta=1, 10)/ &
      0.667070_kind_phys, 0.678820_kind_phys, 0.705538_kind_phys, 0.739591_kind_phys, 0.776758_kind_phys, &
      0.814118_kind_phys, 0.849415_kind_phys, 0.881020_kind_phys, 0.908006_kind_phys, 0.930121_kind_phys/
    data(bm0ij(2, 3, ibeta), ibeta=1, 10)/ &
      0.686356_kind_phys, 0.697839_kind_phys, 0.723997_kind_phys, 0.757285_kind_phys, 0.793389_kind_phys, &
      0.829313_kind_phys, 0.862835_kind_phys, 0.892459_kind_phys, 0.917432_kind_phys, 0.937663_kind_phys/
    data(bm0ij(2, 4, ibeta), ibeta=1, 10)/ &
      0.711425_kind_phys, 0.722572_kind_phys, 0.747941_kind_phys, 0.780055_kind_phys, 0.814518_kind_phys, &
      0.848315_kind_phys, 0.879335_kind_phys, 0.906290_kind_phys, 0.928658_kind_phys, 0.946526_kind_phys/
    data(bm0ij(2, 5, ibeta), ibeta=1, 10)/ &
      0.739575_kind_phys, 0.750307_kind_phys, 0.774633_kind_phys, 0.805138_kind_phys, 0.837408_kind_phys, &
      0.868504_kind_phys, 0.896517_kind_phys, 0.920421_kind_phys, 0.939932_kind_phys, 0.955299_kind_phys/
    data(bm0ij(2, 6, ibeta), ibeta=1, 10)/ &
      0.769143_kind_phys, 0.779346_kind_phys, 0.802314_kind_phys, 0.830752_kind_phys, 0.860333_kind_phys, &
      0.888300_kind_phys, 0.913014_kind_phys, 0.933727_kind_phys, 0.950370_kind_phys, 0.963306_kind_phys/
    data(bm0ij(2, 7, ibeta), ibeta=1, 10)/ &
      0.798900_kind_phys, 0.808431_kind_phys, 0.829700_kind_phys, 0.855653_kind_phys, 0.882163_kind_phys, &
      0.906749_kind_phys, 0.928075_kind_phys, 0.945654_kind_phys, 0.959579_kind_phys, 0.970280_kind_phys/
    data(bm0ij(2, 8, ibeta), ibeta=1, 10)/ &
      0.827826_kind_phys, 0.836542_kind_phys, 0.855808_kind_phys, 0.878954_kind_phys, 0.902174_kind_phys, &
      0.923316_kind_phys, 0.941345_kind_phys, 0.955989_kind_phys, 0.967450_kind_phys, 0.976174_kind_phys/
    data(bm0ij(2, 9, ibeta), ibeta=1, 10)/ &
      0.855068_kind_phys, 0.862856_kind_phys, 0.879900_kind_phys, 0.900068_kind_phys, 0.919956_kind_phys, &
      0.937764_kind_phys, 0.952725_kind_phys, 0.964726_kind_phys, 0.974027_kind_phys, 0.981053_kind_phys/
    data(bm0ij(2, 10, ibeta), ibeta=1, 10)/ &
      0.879961_kind_phys, 0.886755_kind_phys, 0.901484_kind_phys, 0.918665_kind_phys, 0.935346_kind_phys, &
      0.950065_kind_phys, 0.962277_kind_phys, 0.971974_kind_phys, 0.979432_kind_phys, 0.985033_kind_phys/
    data(bm0ij(3, 1, ibeta), ibeta=1, 10)/ &
      0.724166_kind_phys, 0.735474_kind_phys, 0.761359_kind_phys, 0.794045_kind_phys, 0.828702_kind_phys, &
      0.862061_kind_phys, 0.891995_kind_phys, 0.917385_kind_phys, 0.937959_kind_phys, 0.954036_kind_phys/
    data(bm0ij(3, 2, ibeta), ibeta=1, 10)/ &
      0.730416_kind_phys, 0.741780_kind_phys, 0.767647_kind_phys, 0.800116_kind_phys, 0.834344_kind_phys, &
      0.867093_kind_phys, 0.896302_kind_phys, 0.920934_kind_phys, 0.940790_kind_phys, 0.956237_kind_phys/
    data(bm0ij(3, 3, ibeta), ibeta=1, 10)/ &
      0.745327_kind_phys, 0.756664_kind_phys, 0.782255_kind_phys, 0.814026_kind_phys, 0.847107_kind_phys, &
      0.878339_kind_phys, 0.905820_kind_phys, 0.928699_kind_phys, 0.946931_kind_phys, 0.960977_kind_phys/
    data(bm0ij(3, 4, ibeta), ibeta=1, 10)/ &
      0.765195_kind_phys, 0.776312_kind_phys, 0.801216_kind_phys, 0.831758_kind_phys, 0.863079_kind_phys, &
      0.892159_kind_phys, 0.917319_kind_phys, 0.937939_kind_phys, 0.954145_kind_phys, 0.966486_kind_phys/
    data(bm0ij(3, 5, ibeta), ibeta=1, 10)/ &
      0.787632_kind_phys, 0.798347_kind_phys, 0.822165_kind_phys, 0.850985_kind_phys, 0.880049_kind_phys, &
      0.906544_kind_phys, 0.929062_kind_phys, 0.947218_kind_phys, 0.961288_kind_phys, 0.971878_kind_phys/
    data(bm0ij(3, 6, ibeta), ibeta=1, 10)/ &
      0.811024_kind_phys, 0.821179_kind_phys, 0.843557_kind_phys, 0.870247_kind_phys, 0.896694_kind_phys, &
      0.920365_kind_phys, 0.940131_kind_phys, 0.955821_kind_phys, 0.967820_kind_phys, 0.976753_kind_phys/
    data(bm0ij(3, 7, ibeta), ibeta=1, 10)/ &
      0.834254_kind_phys, 0.843709_kind_phys, 0.864356_kind_phys, 0.888619_kind_phys, 0.912245_kind_phys, &
      0.933019_kind_phys, 0.950084_kind_phys, 0.963438_kind_phys, 0.973530_kind_phys, 0.980973_kind_phys/
    data(bm0ij(3, 8, ibeta), ibeta=1, 10)/ &
      0.856531_kind_phys, 0.865176_kind_phys, 0.883881_kind_phys, 0.905544_kind_phys, 0.926290_kind_phys, &
      0.944236_kind_phys, 0.958762_kind_phys, 0.969988_kind_phys, 0.978386_kind_phys, 0.984530_kind_phys/
    data(bm0ij(3, 9, ibeta), ibeta=1, 10)/ &
      0.877307_kind_phys, 0.885070_kind_phys, 0.901716_kind_phys, 0.920729_kind_phys, 0.938663_kind_phys, &
      0.953951_kind_phys, 0.966169_kind_phys, 0.975512_kind_phys, 0.982442_kind_phys, 0.987477_kind_phys/
    data(bm0ij(3, 10, ibeta), ibeta=1, 10)/ &
      0.896234_kind_phys, 0.903082_kind_phys, 0.917645_kind_phys, 0.934069_kind_phys, 0.949354_kind_phys, &
      0.962222_kind_phys, 0.972396_kind_phys, 0.980107_kind_phys, 0.985788_kind_phys, 0.989894_kind_phys/
    data(bm0ij(4, 1, ibeta), ibeta=1, 10)/ &
      0.799294_kind_phys, 0.809144_kind_phys, 0.831293_kind_phys, 0.858395_kind_phys, 0.885897_kind_phys, &
      0.911031_kind_phys, 0.932406_kind_phys, 0.949642_kind_phys, 0.963001_kind_phys, 0.973062_kind_phys/
    data(bm0ij(4, 2, ibeta), ibeta=1, 10)/ &
      0.804239_kind_phys, 0.814102_kind_phys, 0.836169_kind_phys, 0.862984_kind_phys, 0.890003_kind_phys, &
      0.914535_kind_phys, 0.935274_kind_phys, 0.951910_kind_phys, 0.964748_kind_phys, 0.974381_kind_phys/
    data(bm0ij(4, 3, ibeta), ibeta=1, 10)/ &
      0.815910_kind_phys, 0.825708_kind_phys, 0.847403_kind_phys, 0.873389_kind_phys, 0.899185_kind_phys, &
      0.922275_kind_phys, 0.941543_kind_phys, 0.956826_kind_phys, 0.968507_kind_phys, 0.977204_kind_phys/
    data(bm0ij(4, 4, ibeta), ibeta=1, 10)/ &
      0.831348_kind_phys, 0.840892_kind_phys, 0.861793_kind_phys, 0.886428_kind_phys, 0.910463_kind_phys, &
      0.931614_kind_phys, 0.948993_kind_phys, 0.962593_kind_phys, 0.972872_kind_phys, 0.980456_kind_phys/
    data(bm0ij(4, 5, ibeta), ibeta=1, 10)/ &
      0.848597_kind_phys, 0.857693_kind_phys, 0.877402_kind_phys, 0.900265_kind_phys, 0.922180_kind_phys, &
      0.941134_kind_phys, 0.956464_kind_phys, 0.968298_kind_phys, 0.977143_kind_phys, 0.983611_kind_phys/
    data(bm0ij(4, 6, ibeta), ibeta=1, 10)/ &
      0.866271_kind_phys, 0.874764_kind_phys, 0.892984_kind_phys, 0.913796_kind_phys, 0.933407_kind_phys, &
      0.950088_kind_phys, 0.963380_kind_phys, 0.973512_kind_phys, 0.981006_kind_phys, 0.986440_kind_phys/
    data(bm0ij(4, 7, ibeta), ibeta=1, 10)/ &
      0.883430_kind_phys, 0.891216_kind_phys, 0.907762_kind_phys, 0.926388_kind_phys, 0.943660_kind_phys, &
      0.958127_kind_phys, 0.969499_kind_phys, 0.978070_kind_phys, 0.984351_kind_phys, 0.988872_kind_phys/
    data(bm0ij(4, 8, ibeta), ibeta=1, 10)/ &
      0.899483_kind_phys, 0.906505_kind_phys, 0.921294_kind_phys, 0.937719_kind_phys, 0.952729_kind_phys, &
      0.965131_kind_phys, 0.974762_kind_phys, 0.981950_kind_phys, 0.987175_kind_phys, 0.990912_kind_phys/
    data(bm0ij(4, 9, ibeta), ibeta=1, 10)/ &
      0.914096_kind_phys, 0.920337_kind_phys, 0.933373_kind_phys, 0.947677_kind_phys, 0.960579_kind_phys, &
      0.971111_kind_phys, 0.979206_kind_phys, 0.985196_kind_phys, 0.989520_kind_phys, 0.992597_kind_phys/
    data(bm0ij(4, 10, ibeta), ibeta=1, 10)/ &
      0.927122_kind_phys, 0.932597_kind_phys, 0.943952_kind_phys, 0.956277_kind_phys, 0.967268_kind_phys, &
      0.976147_kind_phys, 0.982912_kind_phys, 0.987882_kind_phys, 0.991450_kind_phys, 0.993976_kind_phys/
    data(bm0ij(5, 1, ibeta), ibeta=1, 10)/ &
      0.865049_kind_phys, 0.872851_kind_phys, 0.889900_kind_phys, 0.909907_kind_phys, 0.929290_kind_phys, &
      0.946205_kind_phys, 0.959991_kind_phys, 0.970706_kind_phys, 0.978764_kind_phys, 0.984692_kind_phys/
    data(bm0ij(5, 2, ibeta), ibeta=1, 10)/ &
      0.868989_kind_phys, 0.876713_kind_phys, 0.893538_kind_phys, 0.913173_kind_phys, 0.932080_kind_phys, &
      0.948484_kind_phys, 0.961785_kind_phys, 0.972080_kind_phys, 0.979796_kind_phys, 0.985457_kind_phys/
    data(bm0ij(5, 3, ibeta), ibeta=1, 10)/ &
      0.878010_kind_phys, 0.885524_kind_phys, 0.901756_kind_phys, 0.920464_kind_phys, 0.938235_kind_phys, &
      0.953461_kind_phys, 0.965672_kind_phys, 0.975037_kind_phys, 0.982005_kind_phys, 0.987085_kind_phys/
    data(bm0ij(5, 4, ibeta), ibeta=1, 10)/ &
      0.889534_kind_phys, 0.896698_kind_phys, 0.912012_kind_phys, 0.929395_kind_phys, 0.945647_kind_phys, &
      0.959366_kind_phys, 0.970227_kind_phys, 0.978469_kind_phys, 0.984547_kind_phys, 0.988950_kind_phys/
    data(bm0ij(5, 5, ibeta), ibeta=1, 10)/ &
      0.902033_kind_phys, 0.908713_kind_phys, 0.922848_kind_phys, 0.938648_kind_phys, 0.953186_kind_phys, &
      0.965278_kind_phys, 0.974729_kind_phys, 0.981824_kind_phys, 0.987013_kind_phys, 0.990746_kind_phys/
    data(bm0ij(5, 6, ibeta), ibeta=1, 10)/ &
      0.914496_kind_phys, 0.920599_kind_phys, 0.933389_kind_phys, 0.947485_kind_phys, 0.960262_kind_phys, &
      0.970743_kind_phys, 0.978839_kind_phys, 0.984858_kind_phys, 0.989225_kind_phys, 0.992348_kind_phys/
    data(bm0ij(5, 7, ibeta), ibeta=1, 10)/ &
      0.926281_kind_phys, 0.931761_kind_phys, 0.943142_kind_phys, 0.955526_kind_phys, 0.966600_kind_phys, &
      0.975573_kind_phys, 0.982431_kind_phys, 0.987485_kind_phys, 0.991128_kind_phys, 0.993718_kind_phys/
    data(bm0ij(5, 8, ibeta), ibeta=1, 10)/ &
      0.937029_kind_phys, 0.941877_kind_phys, 0.951868_kind_phys, 0.962615_kind_phys, 0.972112_kind_phys, &
      0.979723_kind_phys, 0.985488_kind_phys, 0.989705_kind_phys, 0.992725_kind_phys, 0.994863_kind_phys/
    data(bm0ij(5, 9, ibeta), ibeta=1, 10)/ &
      0.946580_kind_phys, 0.950819_kind_phys, 0.959494_kind_phys, 0.968732_kind_phys, 0.976811_kind_phys, &
      0.983226_kind_phys, 0.988047_kind_phys, 0.991550_kind_phys, 0.994047_kind_phys, 0.995806_kind_phys/
    data(bm0ij(5, 10, ibeta), ibeta=1, 10)/ &
      0.954909_kind_phys, 0.958581_kind_phys, 0.966049_kind_phys, 0.973933_kind_phys, 0.980766_kind_phys, &
      0.986149_kind_phys, 0.990166_kind_phys, 0.993070_kind_phys, 0.995130_kind_phys, 0.996577_kind_phys/
    data(bm0ij(6, 1, ibeta), ibeta=1, 10)/ &
      0.914182_kind_phys, 0.919824_kind_phys, 0.931832_kind_phys, 0.945387_kind_phys, 0.957999_kind_phys, &
      0.968606_kind_phys, 0.976982_kind_phys, 0.983331_kind_phys, 0.988013_kind_phys, 0.991407_kind_phys/
    data(bm0ij(6, 2, ibeta), ibeta=1, 10)/ &
      0.917139_kind_phys, 0.922665_kind_phys, 0.934395_kind_phys, 0.947580_kind_phys, 0.959792_kind_phys, &
      0.970017_kind_phys, 0.978062_kind_phys, 0.984138_kind_phys, 0.988609_kind_phys, 0.991843_kind_phys/
    data(bm0ij(6, 3, ibeta), ibeta=1, 10)/ &
      0.923742_kind_phys, 0.928990_kind_phys, 0.940064_kind_phys, 0.952396_kind_phys, 0.963699_kind_phys, &
      0.973070_kind_phys, 0.980381_kind_phys, 0.985866_kind_phys, 0.989878_kind_phys, 0.992768_kind_phys/
    data(bm0ij(6, 4, ibeta), ibeta=1, 10)/ &
      0.931870_kind_phys, 0.936743_kind_phys, 0.946941_kind_phys, 0.958162_kind_phys, 0.968318_kind_phys, &
      0.976640_kind_phys, 0.983069_kind_phys, 0.987853_kind_phys, 0.991330_kind_phys, 0.993822_kind_phys/
    data(bm0ij(6, 5, ibeta), ibeta=1, 10)/ &
      0.940376_kind_phys, 0.944807_kind_phys, 0.954004_kind_phys, 0.963999_kind_phys, 0.972928_kind_phys, &
      0.980162_kind_phys, 0.985695_kind_phys, 0.989779_kind_phys, 0.992729_kind_phys, 0.994833_kind_phys/
    data(bm0ij(6, 6, ibeta), ibeta=1, 10)/ &
      0.948597_kind_phys, 0.952555_kind_phys, 0.960703_kind_phys, 0.969454_kind_phys, 0.977181_kind_phys, &
      0.983373_kind_phys, 0.988067_kind_phys, 0.991507_kind_phys, 0.993977_kind_phys, 0.995730_kind_phys/
    data(bm0ij(6, 7, ibeta), ibeta=1, 10)/ &
      0.956167_kind_phys, 0.959648_kind_phys, 0.966763_kind_phys, 0.974326_kind_phys, 0.980933_kind_phys, &
      0.986177_kind_phys, 0.990121_kind_phys, 0.992993_kind_phys, 0.995045_kind_phys, 0.996495_kind_phys/
    data(bm0ij(6, 8, ibeta), ibeta=1, 10)/ &
      0.962913_kind_phys, 0.965937_kind_phys, 0.972080_kind_phys, 0.978552_kind_phys, 0.984153_kind_phys, &
      0.988563_kind_phys, 0.991857_kind_phys, 0.994242_kind_phys, 0.995938_kind_phys, 0.997133_kind_phys/
    data(bm0ij(6, 9, ibeta), ibeta=1, 10)/ &
      0.968787_kind_phys, 0.971391_kind_phys, 0.976651_kind_phys, 0.982148_kind_phys, 0.986869_kind_phys, &
      0.990560_kind_phys, 0.993301_kind_phys, 0.995275_kind_phys, 0.996675_kind_phys, 0.997657_kind_phys/
    data(bm0ij(6, 10, ibeta), ibeta=1, 10)/ &
      0.973822_kind_phys, 0.976047_kind_phys, 0.980523_kind_phys, 0.985170_kind_phys, 0.989134_kind_phys, &
      0.992215_kind_phys, 0.994491_kind_phys, 0.996124_kind_phys, 0.997277_kind_phys, 0.998085_kind_phys/
    data(bm0ij(7, 1, ibeta), ibeta=1, 10)/ &
      0.947410_kind_phys, 0.951207_kind_phys, 0.959119_kind_phys, 0.967781_kind_phys, 0.975592_kind_phys, &
      0.981981_kind_phys, 0.986915_kind_phys, 0.990590_kind_phys, 0.993266_kind_phys, 0.995187_kind_phys/
    data(bm0ij(7, 2, ibeta), ibeta=1, 10)/ &
      0.949477_kind_phys, 0.953161_kind_phys, 0.960824_kind_phys, 0.969187_kind_phys, 0.976702_kind_phys, &
      0.982831_kind_phys, 0.987550_kind_phys, 0.991057_kind_phys, 0.993606_kind_phys, 0.995434_kind_phys/
    data(bm0ij(7, 3, ibeta), ibeta=1, 10)/ &
      0.954008_kind_phys, 0.957438_kind_phys, 0.964537_kind_phys, 0.972232_kind_phys, 0.979095_kind_phys, &
      0.984653_kind_phys, 0.988907_kind_phys, 0.992053_kind_phys, 0.994330_kind_phys, 0.995958_kind_phys/
    data(bm0ij(7, 4, ibeta), ibeta=1, 10)/ &
      0.959431_kind_phys, 0.962539_kind_phys, 0.968935_kind_phys, 0.975808_kind_phys, 0.981882_kind_phys, &
      0.986759_kind_phys, 0.990466_kind_phys, 0.993190_kind_phys, 0.995153_kind_phys, 0.996552_kind_phys/
    data(bm0ij(7, 5, ibeta), ibeta=1, 10)/ &
      0.964932_kind_phys, 0.967693_kind_phys, 0.973342_kind_phys, 0.979355_kind_phys, 0.984620_kind_phys, &
      0.988812_kind_phys, 0.991974_kind_phys, 0.994285_kind_phys, 0.995943_kind_phys, 0.997119_kind_phys/
    data(bm0ij(7, 6, ibeta), ibeta=1, 10)/ &
      0.970101_kind_phys, 0.972517_kind_phys, 0.977428_kind_phys, 0.982612_kind_phys, 0.987110_kind_phys, &
      0.990663_kind_phys, 0.993326_kind_phys, 0.995261_kind_phys, 0.996644_kind_phys, 0.997621_kind_phys/
    data(bm0ij(7, 7, ibeta), ibeta=1, 10)/ &
      0.974746_kind_phys, 0.976834_kind_phys, 0.981055_kind_phys, 0.985475_kind_phys, 0.989280_kind_phys, &
      0.992265_kind_phys, 0.994488_kind_phys, 0.996097_kind_phys, 0.997241_kind_phys, 0.998048_kind_phys/
    data(bm0ij(7, 8, ibeta), ibeta=1, 10)/ &
      0.978804_kind_phys, 0.980591_kind_phys, 0.984187_kind_phys, 0.987927_kind_phys, 0.991124_kind_phys, &
      0.993617_kind_phys, 0.995464_kind_phys, 0.996795_kind_phys, 0.997739_kind_phys, 0.998403_kind_phys/
    data(bm0ij(7, 9, ibeta), ibeta=1, 10)/ &
      0.982280_kind_phys, 0.983799_kind_phys, 0.986844_kind_phys, 0.989991_kind_phys, 0.992667_kind_phys, &
      0.994742_kind_phys, 0.996273_kind_phys, 0.997372_kind_phys, 0.998149_kind_phys, 0.998695_kind_phys/
    data(bm0ij(7, 10, ibeta), ibeta=1, 10)/ &
      0.985218_kind_phys, 0.986503_kind_phys, 0.989071_kind_phys, 0.991711_kind_phys, 0.993945_kind_phys, &
      0.995669_kind_phys, 0.996937_kind_phys, 0.997844_kind_phys, 0.998484_kind_phys, 0.998932_kind_phys/
    data(bm0ij(8, 1, ibeta), ibeta=1, 10)/ &
      0.968507_kind_phys, 0.970935_kind_phys, 0.975916_kind_phys, 0.981248_kind_phys, 0.985947_kind_phys, &
      0.989716_kind_phys, 0.992580_kind_phys, 0.994689_kind_phys, 0.996210_kind_phys, 0.997297_kind_phys/
    data(bm0ij(8, 2, ibeta), ibeta=1, 10)/ &
      0.969870_kind_phys, 0.972210_kind_phys, 0.977002_kind_phys, 0.982119_kind_phys, 0.986619_kind_phys, &
      0.990219_kind_phys, 0.992951_kind_phys, 0.994958_kind_phys, 0.996405_kind_phys, 0.997437_kind_phys/
    data(bm0ij(8, 3, ibeta), ibeta=1, 10)/ &
      0.972820_kind_phys, 0.974963_kind_phys, 0.979339_kind_phys, 0.983988_kind_phys, 0.988054_kind_phys, &
      0.991292_kind_phys, 0.993738_kind_phys, 0.995529_kind_phys, 0.996817_kind_phys, 0.997734_kind_phys/
    data(bm0ij(8, 4, ibeta), ibeta=1, 10)/ &
      0.976280_kind_phys, 0.978186_kind_phys, 0.982060_kind_phys, 0.986151_kind_phys, 0.989706_kind_phys, &
      0.992520_kind_phys, 0.994636_kind_phys, 0.996179_kind_phys, 0.997284_kind_phys, 0.998069_kind_phys/
    data(bm0ij(8, 5, ibeta), ibeta=1, 10)/ &
      0.979711_kind_phys, 0.981372_kind_phys, 0.984735_kind_phys, 0.988263_kind_phys, 0.991309_kind_phys, &
      0.993706_kind_phys, 0.995499_kind_phys, 0.996801_kind_phys, 0.997730_kind_phys, 0.998389_kind_phys/
    data(bm0ij(8, 6, ibeta), ibeta=1, 10)/ &
      0.982863_kind_phys, 0.984292_kind_phys, 0.987172_kind_phys, 0.990174_kind_phys, 0.992750_kind_phys, &
      0.994766_kind_phys, 0.996266_kind_phys, 0.997352_kind_phys, 0.998125_kind_phys, 0.998670_kind_phys/
    data(bm0ij(8, 7, ibeta), ibeta=1, 10)/ &
      0.985642_kind_phys, 0.986858_kind_phys, 0.989301_kind_phys, 0.991834_kind_phys, 0.993994_kind_phys, &
      0.995676_kind_phys, 0.996923_kind_phys, 0.997822_kind_phys, 0.998460_kind_phys, 0.998910_kind_phys/
    data(bm0ij(8, 8, ibeta), ibeta=1, 10)/ &
      0.988029_kind_phys, 0.989058_kind_phys, 0.991116_kind_phys, 0.993240_kind_phys, 0.995043_kind_phys, &
      0.996440_kind_phys, 0.997472_kind_phys, 0.998214_kind_phys, 0.998739_kind_phys, 0.999108_kind_phys/
    data(bm0ij(8, 9, ibeta), ibeta=1, 10)/ &
      0.990046_kind_phys, 0.990912_kind_phys, 0.992640_kind_phys, 0.994415_kind_phys, 0.995914_kind_phys, &
      0.997073_kind_phys, 0.997925_kind_phys, 0.998536_kind_phys, 0.998968_kind_phys, 0.999271_kind_phys/
    data(bm0ij(8, 10, ibeta), ibeta=1, 10)/ &
      0.991732_kind_phys, 0.992459_kind_phys, 0.993906_kind_phys, 0.995386_kind_phys, 0.996633_kind_phys, &
      0.997592_kind_phys, 0.998296_kind_phys, 0.998799_kind_phys, 0.999154_kind_phys, 0.999403_kind_phys/
    data(bm0ij(9, 1, ibeta), ibeta=1, 10)/ &
      0.981392_kind_phys, 0.982893_kind_phys, 0.985938_kind_phys, 0.989146_kind_phys, 0.991928_kind_phys, &
      0.994129_kind_phys, 0.995783_kind_phys, 0.996991_kind_phys, 0.997857_kind_phys, 0.998473_kind_phys/
    data(bm0ij(9, 2, ibeta), ibeta=1, 10)/ &
      0.982254_kind_phys, 0.983693_kind_phys, 0.986608_kind_phys, 0.989673_kind_phys, 0.992328_kind_phys, &
      0.994424_kind_phys, 0.995998_kind_phys, 0.997146_kind_phys, 0.997969_kind_phys, 0.998553_kind_phys/
    data(bm0ij(9, 3, ibeta), ibeta=1, 10)/ &
      0.984104_kind_phys, 0.985407_kind_phys, 0.988040_kind_phys, 0.990798_kind_phys, 0.993178_kind_phys, &
      0.995052_kind_phys, 0.996454_kind_phys, 0.997474_kind_phys, 0.998204_kind_phys, 0.998722_kind_phys/
    data(bm0ij(9, 4, ibeta), ibeta=1, 10)/ &
      0.986243_kind_phys, 0.987386_kind_phys, 0.989687_kind_phys, 0.992087_kind_phys, 0.994149_kind_phys, &
      0.995765_kind_phys, 0.996971_kind_phys, 0.997846_kind_phys, 0.998470_kind_phys, 0.998913_kind_phys/
    data(bm0ij(9, 5, ibeta), ibeta=1, 10)/ &
      0.988332_kind_phys, 0.989313_kind_phys, 0.991284_kind_phys, 0.993332_kind_phys, 0.995082_kind_phys, &
      0.996449_kind_phys, 0.997465_kind_phys, 0.998200_kind_phys, 0.998723_kind_phys, 0.999093_kind_phys/
    data(bm0ij(9, 6, ibeta), ibeta=1, 10)/ &
      0.990220_kind_phys, 0.991053_kind_phys, 0.992721_kind_phys, 0.994445_kind_phys, 0.995914_kind_phys, &
      0.997056_kind_phys, 0.997902_kind_phys, 0.998513_kind_phys, 0.998947_kind_phys, 0.999253_kind_phys/
    data(bm0ij(9, 7, ibeta), ibeta=1, 10)/ &
      0.991859_kind_phys, 0.992561_kind_phys, 0.993961_kind_phys, 0.995403_kind_phys, 0.996626_kind_phys, &
      0.997574_kind_phys, 0.998274_kind_phys, 0.998778_kind_phys, 0.999136_kind_phys, 0.999387_kind_phys/
    data(bm0ij(9, 8, ibeta), ibeta=1, 10)/ &
      0.993250_kind_phys, 0.993837_kind_phys, 0.995007_kind_phys, 0.996208_kind_phys, 0.997223_kind_phys, &
      0.998007_kind_phys, 0.998584_kind_phys, 0.998999_kind_phys, 0.999293_kind_phys, 0.999499_kind_phys/
    data(bm0ij(9, 9, ibeta), ibeta=1, 10)/ &
      0.994413_kind_phys, 0.994903_kind_phys, 0.995878_kind_phys, 0.996876_kind_phys, 0.997716_kind_phys, &
      0.998363_kind_phys, 0.998839_kind_phys, 0.999180_kind_phys, 0.999421_kind_phys, 0.999591_kind_phys/
    data(bm0ij(9, 10, ibeta), ibeta=1, 10)/ &
      0.995376_kind_phys, 0.995785_kind_phys, 0.996597_kind_phys, 0.997425_kind_phys, 0.998121_kind_phys, &
      0.998655_kind_phys, 0.999048_kind_phys, 0.999328_kind_phys, 0.999526_kind_phys, 0.999665_kind_phys/
    data(bm0ij(10, 1, ibeta), ibeta=1, 10)/ &
      0.989082_kind_phys, 0.989991_kind_phys, 0.991819_kind_phys, 0.993723_kind_phys, 0.995357_kind_phys, &
      0.996637_kind_phys, 0.997592_kind_phys, 0.998286_kind_phys, 0.998781_kind_phys, 0.999132_kind_phys/
    data(bm0ij(10, 2, ibeta), ibeta=1, 10)/ &
      0.989613_kind_phys, 0.990480_kind_phys, 0.992224_kind_phys, 0.994039_kind_phys, 0.995594_kind_phys, &
      0.996810_kind_phys, 0.997717_kind_phys, 0.998375_kind_phys, 0.998845_kind_phys, 0.999178_kind_phys/
    data(bm0ij(10, 3, ibeta), ibeta=1, 10)/ &
      0.990744_kind_phys, 0.991523_kind_phys, 0.993086_kind_phys, 0.994708_kind_phys, 0.996094_kind_phys, &
      0.997176_kind_phys, 0.997981_kind_phys, 0.998564_kind_phys, 0.998980_kind_phys, 0.999274_kind_phys/
    data(bm0ij(10, 4, ibeta), ibeta=1, 10)/ &
      0.992041_kind_phys, 0.992716_kind_phys, 0.994070_kind_phys, 0.995470_kind_phys, 0.996662_kind_phys, &
      0.997591_kind_phys, 0.998280_kind_phys, 0.998778_kind_phys, 0.999133_kind_phys, 0.999383_kind_phys/
    data(bm0ij(10, 5, ibeta), ibeta=1, 10)/ &
      0.993292_kind_phys, 0.993867_kind_phys, 0.995015_kind_phys, 0.996199_kind_phys, 0.997205_kind_phys, &
      0.997985_kind_phys, 0.998564_kind_phys, 0.998981_kind_phys, 0.999277_kind_phys, 0.999487_kind_phys/
    data(bm0ij(10, 6, ibeta), ibeta=1, 10)/ &
      0.994411_kind_phys, 0.994894_kind_phys, 0.995857_kind_phys, 0.996847_kind_phys, 0.997685_kind_phys, &
      0.998334_kind_phys, 0.998814_kind_phys, 0.999159_kind_phys, 0.999404_kind_phys, 0.999577_kind_phys/
    data(bm0ij(10, 7, ibeta), ibeta=1, 10)/ &
      0.995373_kind_phys, 0.995776_kind_phys, 0.996577_kind_phys, 0.997400_kind_phys, 0.998094_kind_phys, &
      0.998630_kind_phys, 0.999026_kind_phys, 0.999310_kind_phys, 0.999512_kind_phys, 0.999654_kind_phys/
    data(bm0ij(10, 8, ibeta), ibeta=1, 10)/ &
      0.996181_kind_phys, 0.996516_kind_phys, 0.997181_kind_phys, 0.997861_kind_phys, 0.998435_kind_phys, &
      0.998877_kind_phys, 0.999202_kind_phys, 0.999435_kind_phys, 0.999601_kind_phys, 0.999717_kind_phys/
    data(bm0ij(10, 9, ibeta), ibeta=1, 10)/ &
      0.996851_kind_phys, 0.997128_kind_phys, 0.997680_kind_phys, 0.998242_kind_phys, 0.998715_kind_phys, &
      0.999079_kind_phys, 0.999346_kind_phys, 0.999538_kind_phys, 0.999673_kind_phys, 0.999769_kind_phys/
    data(bm0ij(10, 10, ibeta), ibeta=1, 10)/ &
      0.997402_kind_phys, 0.997632_kind_phys, 0.998089_kind_phys, 0.998554_kind_phys, 0.998945_kind_phys, &
      0.999244_kind_phys, 0.999464_kind_phys, 0.999622_kind_phys, 0.999733_kind_phys, 0.999811_kind_phys/

! rpm....   3rd moment nuclei mode corr. fac. for bimodal fm coag rate

    data(bm3i(1, 1, ibeta), ibeta=1, 10)/ &
      0.70708_kind_phys, 0.71681_kind_phys, 0.73821_kind_phys, 0.76477_kind_phys, 0.79350_kind_phys, 0.82265_kind_phys, 0.85090_kind_phys, 0.87717_kind_phys, &
      0.90069_kind_phys, 0.92097_kind_phys/
    data(bm3i(1, 2, ibeta), ibeta=1, 10)/ &
      0.72172_kind_phys, 0.73022_kind_phys, 0.74927_kind_phys, 0.77324_kind_phys, 0.79936_kind_phys, 0.82601_kind_phys, 0.85199_kind_phys, 0.87637_kind_phys, &
      0.89843_kind_phys, 0.91774_kind_phys/
    data(bm3i(1, 3, ibeta), ibeta=1, 10)/ &
      0.78291_kind_phys, 0.78896_kind_phys, 0.80286_kind_phys, 0.82070_kind_phys, 0.84022_kind_phys, 0.85997_kind_phys, 0.87901_kind_phys, 0.89669_kind_phys, &
      0.91258_kind_phys, 0.92647_kind_phys/
    data(bm3i(1, 4, ibeta), ibeta=1, 10)/ &
      0.87760_kind_phys, 0.88147_kind_phys, 0.89025_kind_phys, 0.90127_kind_phys, 0.91291_kind_phys, 0.92420_kind_phys, 0.93452_kind_phys, 0.94355_kind_phys, &
      0.95113_kind_phys, 0.95726_kind_phys/
    data(bm3i(1, 5, ibeta), ibeta=1, 10)/ &
      0.94988_kind_phys, 0.95184_kind_phys, 0.95612_kind_phys, 0.96122_kind_phys, 0.96628_kind_phys, 0.97085_kind_phys, 0.97467_kind_phys, 0.97763_kind_phys, &
      0.97971_kind_phys, 0.98089_kind_phys/
    data(bm3i(1, 6, ibeta), ibeta=1, 10)/ &
      0.98318_kind_phys, 0.98393_kind_phys, 0.98551_kind_phys, 0.98728_kind_phys, 0.98889_kind_phys, 0.99014_kind_phys, 0.99095_kind_phys, 0.99124_kind_phys, &
      0.99100_kind_phys, 0.99020_kind_phys/
    data(bm3i(1, 7, ibeta), ibeta=1, 10)/ &
      0.99480_kind_phys, 0.99504_kind_phys, 0.99551_kind_phys, 0.99598_kind_phys, 0.99629_kind_phys, 0.99635_kind_phys, 0.99611_kind_phys, 0.99550_kind_phys, &
      0.99450_kind_phys, 0.99306_kind_phys/
    data(bm3i(1, 8, ibeta), ibeta=1, 10)/ &
      0.99842_kind_phys, 0.99848_kind_phys, 0.99858_kind_phys, 0.99861_kind_phys, 0.99850_kind_phys, 0.99819_kind_phys, 0.99762_kind_phys, 0.99674_kind_phys, &
      0.99550_kind_phys, 0.99388_kind_phys/
    data(bm3i(1, 9, ibeta), ibeta=1, 10)/ &
      0.99951_kind_phys, 0.99951_kind_phys, 0.99949_kind_phys, 0.99939_kind_phys, 0.99915_kind_phys, 0.99872_kind_phys, 0.99805_kind_phys, 0.99709_kind_phys, &
      0.99579_kind_phys, 0.99411_kind_phys/
    data(bm3i(1, 10, ibeta), ibeta=1, 10)/ &
      0.99984_kind_phys, 0.99982_kind_phys, 0.99976_kind_phys, 0.99962_kind_phys, 0.99934_kind_phys, 0.99888_kind_phys, 0.99818_kind_phys, 0.99719_kind_phys, &
      0.99587_kind_phys, 0.99417_kind_phys/
    data(bm3i(2, 1, ibeta), ibeta=1, 10)/ &
      0.72957_kind_phys, 0.73993_kind_phys, 0.76303_kind_phys, 0.79178_kind_phys, 0.82245_kind_phys, 0.85270_kind_phys, 0.88085_kind_phys, 0.90578_kind_phys, &
      0.92691_kind_phys, 0.94415_kind_phys/
    data(bm3i(2, 2, ibeta), ibeta=1, 10)/ &
      0.72319_kind_phys, 0.73320_kind_phys, 0.75547_kind_phys, 0.78323_kind_phys, 0.81307_kind_phys, 0.84287_kind_phys, 0.87107_kind_phys, 0.89651_kind_phys, &
      0.91852_kind_phys, 0.93683_kind_phys/
    data(bm3i(2, 3, ibeta), ibeta=1, 10)/ &
      0.74413_kind_phys, 0.75205_kind_phys, 0.76998_kind_phys, 0.79269_kind_phys, 0.81746_kind_phys, 0.84258_kind_phys, 0.86685_kind_phys, 0.88938_kind_phys, &
      0.90953_kind_phys, 0.92695_kind_phys/
    data(bm3i(2, 4, ibeta), ibeta=1, 10)/ &
      0.82588_kind_phys, 0.83113_kind_phys, 0.84309_kind_phys, 0.85825_kind_phys, 0.87456_kind_phys, 0.89072_kind_phys, 0.90594_kind_phys, 0.91972_kind_phys, &
      0.93178_kind_phys, 0.94203_kind_phys/
    data(bm3i(2, 5, ibeta), ibeta=1, 10)/ &
      0.91886_kind_phys, 0.92179_kind_phys, 0.92831_kind_phys, 0.93624_kind_phys, 0.94434_kind_phys, 0.95192_kind_phys, 0.95856_kind_phys, 0.96409_kind_phys, &
      0.96845_kind_phys, 0.97164_kind_phys/
    data(bm3i(2, 6, ibeta), ibeta=1, 10)/ &
      0.97129_kind_phys, 0.97252_kind_phys, 0.97515_kind_phys, 0.97818_kind_phys, 0.98108_kind_phys, 0.98354_kind_phys, 0.98542_kind_phys, 0.98665_kind_phys, &
      0.98721_kind_phys, 0.98709_kind_phys/
    data(bm3i(2, 7, ibeta), ibeta=1, 10)/ &
      0.99104_kind_phys, 0.99145_kind_phys, 0.99230_kind_phys, 0.99320_kind_phys, 0.99394_kind_phys, 0.99439_kind_phys, 0.99448_kind_phys, 0.99416_kind_phys, &
      0.99340_kind_phys, 0.99217_kind_phys/
    data(bm3i(2, 8, ibeta), ibeta=1, 10)/ &
      0.99730_kind_phys, 0.99741_kind_phys, 0.99763_kind_phys, 0.99779_kind_phys, 0.99782_kind_phys, 0.99762_kind_phys, 0.99715_kind_phys, 0.99636_kind_phys, &
      0.99519_kind_phys, 0.99363_kind_phys/
    data(bm3i(2, 9, ibeta), ibeta=1, 10)/ &
      0.99917_kind_phys, 0.99919_kind_phys, 0.99921_kind_phys, 0.99915_kind_phys, 0.99895_kind_phys, 0.99856_kind_phys, 0.99792_kind_phys, 0.99698_kind_phys, &
      0.99570_kind_phys, 0.99404_kind_phys/
    data(bm3i(2, 10, ibeta), ibeta=1, 10)/ &
      0.99973_kind_phys, 0.99973_kind_phys, 0.99968_kind_phys, 0.99955_kind_phys, 0.99928_kind_phys, 0.99883_kind_phys, 0.99814_kind_phys, 0.99716_kind_phys, &
      0.99584_kind_phys, 0.99415_kind_phys/
    data(bm3i(3, 1, ibeta), ibeta=1, 10)/ &
      0.78358_kind_phys, 0.79304_kind_phys, 0.81445_kind_phys, 0.84105_kind_phys, 0.86873_kind_phys, 0.89491_kind_phys, 0.91805_kind_phys, 0.93743_kind_phys, &
      0.95300_kind_phys, 0.96510_kind_phys/
    data(bm3i(3, 2, ibeta), ibeta=1, 10)/ &
      0.76412_kind_phys, 0.77404_kind_phys, 0.79635_kind_phys, 0.82404_kind_phys, 0.85312_kind_phys, 0.88101_kind_phys, 0.90610_kind_phys, 0.92751_kind_phys, &
      0.94500_kind_phys, 0.95879_kind_phys/
    data(bm3i(3, 3, ibeta), ibeta=1, 10)/ &
      0.74239_kind_phys, 0.75182_kind_phys, 0.77301_kind_phys, 0.79956_kind_phys, 0.82809_kind_phys, 0.85639_kind_phys, 0.88291_kind_phys, 0.90658_kind_phys, &
      0.92683_kind_phys, 0.94350_kind_phys/
    data(bm3i(3, 4, ibeta), ibeta=1, 10)/ &
      0.78072_kind_phys, 0.78758_kind_phys, 0.80317_kind_phys, 0.82293_kind_phys, 0.84437_kind_phys, 0.86589_kind_phys, 0.88643_kind_phys, 0.90526_kind_phys, &
      0.92194_kind_phys, 0.93625_kind_phys/
    data(bm3i(3, 5, ibeta), ibeta=1, 10)/ &
      0.87627_kind_phys, 0.88044_kind_phys, 0.88981_kind_phys, 0.90142_kind_phys, 0.91357_kind_phys, 0.92524_kind_phys, 0.93585_kind_phys, 0.94510_kind_phys, &
      0.95285_kind_phys, 0.95911_kind_phys/
    data(bm3i(3, 6, ibeta), ibeta=1, 10)/ &
      0.95176_kind_phys, 0.95371_kind_phys, 0.95796_kind_phys, 0.96297_kind_phys, 0.96792_kind_phys, 0.97233_kind_phys, 0.97599_kind_phys, 0.97880_kind_phys, &
      0.98072_kind_phys, 0.98178_kind_phys/
    data(bm3i(3, 7, ibeta), ibeta=1, 10)/ &
      0.98453_kind_phys, 0.98523_kind_phys, 0.98670_kind_phys, 0.98833_kind_phys, 0.98980_kind_phys, 0.99092_kind_phys, 0.99160_kind_phys, 0.99179_kind_phys, &
      0.99145_kind_phys, 0.99058_kind_phys/
    data(bm3i(3, 8, ibeta), ibeta=1, 10)/ &
      0.99534_kind_phys, 0.99555_kind_phys, 0.99597_kind_phys, 0.99637_kind_phys, 0.99662_kind_phys, 0.99663_kind_phys, 0.99633_kind_phys, 0.99569_kind_phys, &
      0.99465_kind_phys, 0.99318_kind_phys/
    data(bm3i(3, 9, ibeta), ibeta=1, 10)/ &
      0.99859_kind_phys, 0.99864_kind_phys, 0.99872_kind_phys, 0.99873_kind_phys, 0.99860_kind_phys, 0.99827_kind_phys, 0.99768_kind_phys, 0.99679_kind_phys, &
      0.99555_kind_phys, 0.99391_kind_phys/
    data(bm3i(3, 10, ibeta), ibeta=1, 10)/ &
      0.99956_kind_phys, 0.99956_kind_phys, 0.99953_kind_phys, 0.99942_kind_phys, 0.99918_kind_phys, 0.99875_kind_phys, 0.99807_kind_phys, 0.99711_kind_phys, &
      0.99580_kind_phys, 0.99412_kind_phys/
    data(bm3i(4, 1, ibeta), ibeta=1, 10)/ &
      0.84432_kind_phys, 0.85223_kind_phys, 0.86990_kind_phys, 0.89131_kind_phys, 0.91280_kind_phys, 0.93223_kind_phys, 0.94861_kind_phys, 0.96172_kind_phys, &
      0.97185_kind_phys, 0.97945_kind_phys/
    data(bm3i(4, 2, ibeta), ibeta=1, 10)/ &
      0.82299_kind_phys, 0.83164_kind_phys, 0.85101_kind_phys, 0.87463_kind_phys, 0.89857_kind_phys, 0.92050_kind_phys, 0.93923_kind_phys, 0.95443_kind_phys, &
      0.96629_kind_phys, 0.97529_kind_phys/
    data(bm3i(4, 3, ibeta), ibeta=1, 10)/ &
      0.77870_kind_phys, 0.78840_kind_phys, 0.81011_kind_phys, 0.83690_kind_phys, 0.86477_kind_phys, 0.89124_kind_phys, 0.91476_kind_phys, 0.93460_kind_phys, &
      0.95063_kind_phys, 0.96316_kind_phys/
    data(bm3i(4, 4, ibeta), ibeta=1, 10)/ &
      0.76386_kind_phys, 0.77233_kind_phys, 0.79147_kind_phys, 0.81557_kind_phys, 0.84149_kind_phys, 0.86719_kind_phys, 0.89126_kind_phys, 0.91275_kind_phys, &
      0.93116_kind_phys, 0.94637_kind_phys/
    data(bm3i(4, 5, ibeta), ibeta=1, 10)/ &
      0.82927_kind_phys, 0.83488_kind_phys, 0.84756_kind_phys, 0.86346_kind_phys, 0.88040_kind_phys, 0.89704_kind_phys, 0.91257_kind_phys, 0.92649_kind_phys, &
      0.93857_kind_phys, 0.94874_kind_phys/
    data(bm3i(4, 6, ibeta), ibeta=1, 10)/ &
      0.92184_kind_phys, 0.92481_kind_phys, 0.93136_kind_phys, 0.93925_kind_phys, 0.94724_kind_phys, 0.95462_kind_phys, 0.96104_kind_phys, 0.96634_kind_phys, &
      0.97048_kind_phys, 0.97348_kind_phys/
    data(bm3i(4, 7, ibeta), ibeta=1, 10)/ &
      0.97341_kind_phys, 0.97457_kind_phys, 0.97706_kind_phys, 0.97991_kind_phys, 0.98260_kind_phys, 0.98485_kind_phys, 0.98654_kind_phys, 0.98760_kind_phys, &
      0.98801_kind_phys, 0.98777_kind_phys/
    data(bm3i(4, 8, ibeta), ibeta=1, 10)/ &
      0.99192_kind_phys, 0.99229_kind_phys, 0.99305_kind_phys, 0.99385_kind_phys, 0.99449_kind_phys, 0.99486_kind_phys, 0.99487_kind_phys, 0.99449_kind_phys, &
      0.99367_kind_phys, 0.99239_kind_phys/
    data(bm3i(4, 9, ibeta), ibeta=1, 10)/ &
      0.99758_kind_phys, 0.99768_kind_phys, 0.99787_kind_phys, 0.99800_kind_phys, 0.99799_kind_phys, 0.99777_kind_phys, 0.99727_kind_phys, 0.99645_kind_phys, &
      0.99527_kind_phys, 0.99369_kind_phys/
    data(bm3i(4, 10, ibeta), ibeta=1, 10)/ &
      0.99926_kind_phys, 0.99928_kind_phys, 0.99928_kind_phys, 0.99921_kind_phys, 0.99900_kind_phys, 0.99860_kind_phys, 0.99795_kind_phys, 0.99701_kind_phys, &
      0.99572_kind_phys, 0.99405_kind_phys/
    data(bm3i(5, 1, ibeta), ibeta=1, 10)/ &
      0.89577_kind_phys, 0.90190_kind_phys, 0.91522_kind_phys, 0.93076_kind_phys, 0.94575_kind_phys, 0.95876_kind_phys, 0.96932_kind_phys, 0.97751_kind_phys, &
      0.98367_kind_phys, 0.98820_kind_phys/
    data(bm3i(5, 2, ibeta), ibeta=1, 10)/ &
      0.87860_kind_phys, 0.88547_kind_phys, 0.90052_kind_phys, 0.91828_kind_phys, 0.93557_kind_phys, 0.95075_kind_phys, 0.96319_kind_phys, 0.97292_kind_phys, &
      0.98028_kind_phys, 0.98572_kind_phys/
    data(bm3i(5, 3, ibeta), ibeta=1, 10)/ &
      0.83381_kind_phys, 0.84240_kind_phys, 0.86141_kind_phys, 0.88425_kind_phys, 0.90707_kind_phys, 0.92770_kind_phys, 0.94510_kind_phys, 0.95906_kind_phys, &
      0.96986_kind_phys, 0.97798_kind_phys/
    data(bm3i(5, 4, ibeta), ibeta=1, 10)/ &
      0.78530_kind_phys, 0.79463_kind_phys, 0.81550_kind_phys, 0.84127_kind_phys, 0.86813_kind_phys, 0.89367_kind_phys, 0.91642_kind_phys, 0.93566_kind_phys, &
      0.95125_kind_phys, 0.96347_kind_phys/
    data(bm3i(5, 5, ibeta), ibeta=1, 10)/ &
      0.79614_kind_phys, 0.80332_kind_phys, 0.81957_kind_phys, 0.84001_kind_phys, 0.86190_kind_phys, 0.88351_kind_phys, 0.90368_kind_phys, 0.92169_kind_phys, &
      0.93718_kind_phys, 0.95006_kind_phys/
    data(bm3i(5, 6, ibeta), ibeta=1, 10)/ &
      0.88192_kind_phys, 0.88617_kind_phys, 0.89565_kind_phys, 0.90728_kind_phys, 0.91931_kind_phys, 0.93076_kind_phys, 0.94107_kind_phys, 0.94997_kind_phys, &
      0.95739_kind_phys, 0.96333_kind_phys/
    data(bm3i(5, 7, ibeta), ibeta=1, 10)/ &
      0.95509_kind_phys, 0.95698_kind_phys, 0.96105_kind_phys, 0.96583_kind_phys, 0.97048_kind_phys, 0.97460_kind_phys, 0.97796_kind_phys, 0.98050_kind_phys, &
      0.98218_kind_phys, 0.98304_kind_phys/
    data(bm3i(5, 8, ibeta), ibeta=1, 10)/ &
      0.98596_kind_phys, 0.98660_kind_phys, 0.98794_kind_phys, 0.98943_kind_phys, 0.99074_kind_phys, 0.99172_kind_phys, 0.99227_kind_phys, 0.99235_kind_phys, &
      0.99192_kind_phys, 0.99096_kind_phys/
    data(bm3i(5, 9, ibeta), ibeta=1, 10)/ &
      0.99581_kind_phys, 0.99600_kind_phys, 0.99637_kind_phys, 0.99672_kind_phys, 0.99691_kind_phys, 0.99687_kind_phys, 0.99653_kind_phys, 0.99585_kind_phys, &
      0.99478_kind_phys, 0.99329_kind_phys/
    data(bm3i(5, 10, ibeta), ibeta=1, 10)/ &
      0.99873_kind_phys, 0.99878_kind_phys, 0.99884_kind_phys, 0.99883_kind_phys, 0.99869_kind_phys, 0.99834_kind_phys, 0.99774_kind_phys, 0.99684_kind_phys, &
      0.99558_kind_phys, 0.99394_kind_phys/
    data(bm3i(6, 1, ibeta), ibeta=1, 10)/ &
      0.93335_kind_phys, 0.93777_kind_phys, 0.94711_kind_phys, 0.95764_kind_phys, 0.96741_kind_phys, 0.97562_kind_phys, 0.98210_kind_phys, 0.98701_kind_phys, &
      0.99064_kind_phys, 0.99327_kind_phys/
    data(bm3i(6, 2, ibeta), ibeta=1, 10)/ &
      0.92142_kind_phys, 0.92646_kind_phys, 0.93723_kind_phys, 0.94947_kind_phys, 0.96096_kind_phys, 0.97069_kind_phys, 0.97842_kind_phys, 0.98431_kind_phys, &
      0.98868_kind_phys, 0.99186_kind_phys/
    data(bm3i(6, 3, ibeta), ibeta=1, 10)/ &
      0.88678_kind_phys, 0.89351_kind_phys, 0.90810_kind_phys, 0.92508_kind_phys, 0.94138_kind_phys, 0.95549_kind_phys, 0.96693_kind_phys, 0.97578_kind_phys, &
      0.98243_kind_phys, 0.98731_kind_phys/
    data(bm3i(6, 4, ibeta), ibeta=1, 10)/ &
      0.83249_kind_phys, 0.84124_kind_phys, 0.86051_kind_phys, 0.88357_kind_phys, 0.90655_kind_phys, 0.92728_kind_phys, 0.94477_kind_phys, 0.95880_kind_phys, &
      0.96964_kind_phys, 0.97779_kind_phys/
    data(bm3i(6, 5, ibeta), ibeta=1, 10)/ &
      0.79593_kind_phys, 0.80444_kind_phys, 0.82355_kind_phys, 0.84725_kind_phys, 0.87211_kind_phys, 0.89593_kind_phys, 0.91735_kind_phys, 0.93566_kind_phys, &
      0.95066_kind_phys, 0.96255_kind_phys/
    data(bm3i(6, 6, ibeta), ibeta=1, 10)/ &
      0.84124_kind_phys, 0.84695_kind_phys, 0.85980_kind_phys, 0.87575_kind_phys, 0.89256_kind_phys, 0.90885_kind_phys, 0.92383_kind_phys, 0.93704_kind_phys, &
      0.94830_kind_phys, 0.95761_kind_phys/
    data(bm3i(6, 7, ibeta), ibeta=1, 10)/ &
      0.92721_kind_phys, 0.93011_kind_phys, 0.93647_kind_phys, 0.94406_kind_phys, 0.95166_kind_phys, 0.95862_kind_phys, 0.96460_kind_phys, 0.96949_kind_phys, &
      0.97326_kind_phys, 0.97595_kind_phys/
    data(bm3i(6, 8, ibeta), ibeta=1, 10)/ &
      0.97573_kind_phys, 0.97681_kind_phys, 0.97913_kind_phys, 0.98175_kind_phys, 0.98421_kind_phys, 0.98624_kind_phys, 0.98772_kind_phys, 0.98860_kind_phys, &
      0.98885_kind_phys, 0.98847_kind_phys/
    data(bm3i(6, 9, ibeta), ibeta=1, 10)/ &
      0.99271_kind_phys, 0.99304_kind_phys, 0.99373_kind_phys, 0.99444_kind_phys, 0.99499_kind_phys, 0.99528_kind_phys, 0.99522_kind_phys, 0.99477_kind_phys, &
      0.99390_kind_phys, 0.99258_kind_phys/
    data(bm3i(6, 10, ibeta), ibeta=1, 10)/ &
      0.99782_kind_phys, 0.99791_kind_phys, 0.99807_kind_phys, 0.99817_kind_phys, 0.99813_kind_phys, 0.99788_kind_phys, 0.99737_kind_phys, 0.99653_kind_phys, &
      0.99533_kind_phys, 0.99374_kind_phys/
    data(bm3i(7, 1, ibeta), ibeta=1, 10)/ &
      0.95858_kind_phys, 0.96158_kind_phys, 0.96780_kind_phys, 0.97460_kind_phys, 0.98073_kind_phys, 0.98575_kind_phys, 0.98963_kind_phys, 0.99252_kind_phys, &
      0.99463_kind_phys, 0.99615_kind_phys/
    data(bm3i(7, 2, ibeta), ibeta=1, 10)/ &
      0.95091_kind_phys, 0.95438_kind_phys, 0.96163_kind_phys, 0.96962_kind_phys, 0.97688_kind_phys, 0.98286_kind_phys, 0.98751_kind_phys, 0.99099_kind_phys, &
      0.99353_kind_phys, 0.99536_kind_phys/
    data(bm3i(7, 3, ibeta), ibeta=1, 10)/ &
      0.92751_kind_phys, 0.93233_kind_phys, 0.94255_kind_phys, 0.95406_kind_phys, 0.96473_kind_phys, 0.97366_kind_phys, 0.98070_kind_phys, 0.98602_kind_phys, &
      0.98994_kind_phys, 0.99278_kind_phys/
    data(bm3i(7, 4, ibeta), ibeta=1, 10)/ &
      0.88371_kind_phys, 0.89075_kind_phys, 0.90595_kind_phys, 0.92351_kind_phys, 0.94028_kind_phys, 0.95474_kind_phys, 0.96642_kind_phys, 0.97544_kind_phys, &
      0.98220_kind_phys, 0.98715_kind_phys/
    data(bm3i(7, 5, ibeta), ibeta=1, 10)/ &
      0.82880_kind_phys, 0.83750_kind_phys, 0.85671_kind_phys, 0.87980_kind_phys, 0.90297_kind_phys, 0.92404_kind_phys, 0.94195_kind_phys, 0.95644_kind_phys, &
      0.96772_kind_phys, 0.97625_kind_phys/
    data(bm3i(7, 6, ibeta), ibeta=1, 10)/ &
      0.81933_kind_phys, 0.82655_kind_phys, 0.84279_kind_phys, 0.86295_kind_phys, 0.88412_kind_phys, 0.90449_kind_phys, 0.92295_kind_phys, 0.93890_kind_phys, &
      0.95215_kind_phys, 0.96281_kind_phys/
    data(bm3i(7, 7, ibeta), ibeta=1, 10)/ &
      0.89099_kind_phys, 0.89519_kind_phys, 0.90448_kind_phys, 0.91577_kind_phys, 0.92732_kind_phys, 0.93820_kind_phys, 0.94789_kind_phys, 0.95616_kind_phys, &
      0.96297_kind_phys, 0.96838_kind_phys/
    data(bm3i(7, 8, ibeta), ibeta=1, 10)/ &
      0.95886_kind_phys, 0.96064_kind_phys, 0.96448_kind_phys, 0.96894_kind_phys, 0.97324_kind_phys, 0.97701_kind_phys, 0.98004_kind_phys, 0.98228_kind_phys, &
      0.98371_kind_phys, 0.98435_kind_phys/
    data(bm3i(7, 9, ibeta), ibeta=1, 10)/ &
      0.98727_kind_phys, 0.98786_kind_phys, 0.98908_kind_phys, 0.99043_kind_phys, 0.99160_kind_phys, 0.99245_kind_phys, 0.99288_kind_phys, 0.99285_kind_phys, &
      0.99234_kind_phys, 0.99131_kind_phys/
    data(bm3i(7, 10, ibeta), ibeta=1, 10)/ &
      0.99621_kind_phys, 0.99638_kind_phys, 0.99671_kind_phys, 0.99700_kind_phys, 0.99715_kind_phys, 0.99707_kind_phys, 0.99670_kind_phys, 0.99599_kind_phys, &
      0.99489_kind_phys, 0.99338_kind_phys/
    data(bm3i(8, 1, ibeta), ibeta=1, 10)/ &
      0.97470_kind_phys, 0.97666_kind_phys, 0.98064_kind_phys, 0.98491_kind_phys, 0.98867_kind_phys, 0.99169_kind_phys, 0.99399_kind_phys, 0.99569_kind_phys, &
      0.99691_kind_phys, 0.99779_kind_phys/
    data(bm3i(8, 2, ibeta), ibeta=1, 10)/ &
      0.96996_kind_phys, 0.97225_kind_phys, 0.97693_kind_phys, 0.98196_kind_phys, 0.98643_kind_phys, 0.99003_kind_phys, 0.99279_kind_phys, 0.99482_kind_phys, &
      0.99630_kind_phys, 0.99735_kind_phys/
    data(bm3i(8, 3, ibeta), ibeta=1, 10)/ &
      0.95523_kind_phys, 0.95848_kind_phys, 0.96522_kind_phys, 0.97260_kind_phys, 0.97925_kind_phys, 0.98468_kind_phys, 0.98888_kind_phys, 0.99200_kind_phys, &
      0.99427_kind_phys, 0.99590_kind_phys/
    data(bm3i(8, 4, ibeta), ibeta=1, 10)/ &
      0.92524_kind_phys, 0.93030_kind_phys, 0.94098_kind_phys, 0.95294_kind_phys, 0.96397_kind_phys, 0.97317_kind_phys, 0.98038_kind_phys, 0.98582_kind_phys, &
      0.98981_kind_phys, 0.99270_kind_phys/
    data(bm3i(8, 5, ibeta), ibeta=1, 10)/ &
      0.87576_kind_phys, 0.88323_kind_phys, 0.89935_kind_phys, 0.91799_kind_phys, 0.93583_kind_phys, 0.95126_kind_phys, 0.96377_kind_phys, 0.97345_kind_phys, &
      0.98072_kind_phys, 0.98606_kind_phys/
    data(bm3i(8, 6, ibeta), ibeta=1, 10)/ &
      0.83078_kind_phys, 0.83894_kind_phys, 0.85705_kind_phys, 0.87899_kind_phys, 0.90126_kind_phys, 0.92179_kind_phys, 0.93950_kind_phys, 0.95404_kind_phys, &
      0.96551_kind_phys, 0.97430_kind_phys/
    data(bm3i(8, 7, ibeta), ibeta=1, 10)/ &
      0.85727_kind_phys, 0.86294_kind_phys, 0.87558_kind_phys, 0.89111_kind_phys, 0.90723_kind_phys, 0.92260_kind_phys, 0.93645_kind_phys, 0.94841_kind_phys, &
      0.95838_kind_phys, 0.96643_kind_phys/
    data(bm3i(8, 8, ibeta), ibeta=1, 10)/ &
      0.93337_kind_phys, 0.93615_kind_phys, 0.94220_kind_phys, 0.94937_kind_phys, 0.95647_kind_phys, 0.96292_kind_phys, 0.96840_kind_phys, 0.97283_kind_phys, &
      0.97619_kind_phys, 0.97854_kind_phys/
    data(bm3i(8, 9, ibeta), ibeta=1, 10)/ &
      0.97790_kind_phys, 0.97891_kind_phys, 0.98105_kind_phys, 0.98346_kind_phys, 0.98569_kind_phys, 0.98751_kind_phys, 0.98879_kind_phys, 0.98950_kind_phys, &
      0.98961_kind_phys, 0.98912_kind_phys/
    data(bm3i(8, 10, ibeta), ibeta=1, 10)/ &
      0.99337_kind_phys, 0.99367_kind_phys, 0.99430_kind_phys, 0.99493_kind_phys, 0.99541_kind_phys, 0.99562_kind_phys, 0.99551_kind_phys, 0.99501_kind_phys, &
      0.99410_kind_phys, 0.99274_kind_phys/
    data(bm3i(9, 1, ibeta), ibeta=1, 10)/ &
      0.98470_kind_phys, 0.98594_kind_phys, 0.98844_kind_phys, 0.99106_kind_phys, 0.99334_kind_phys, 0.99514_kind_phys, 0.99650_kind_phys, 0.99749_kind_phys, &
      0.99821_kind_phys, 0.99872_kind_phys/
    data(bm3i(9, 2, ibeta), ibeta=1, 10)/ &
      0.98184_kind_phys, 0.98330_kind_phys, 0.98624_kind_phys, 0.98934_kind_phys, 0.99205_kind_phys, 0.99420_kind_phys, 0.99582_kind_phys, 0.99701_kind_phys, &
      0.99787_kind_phys, 0.99848_kind_phys/
    data(bm3i(9, 3, ibeta), ibeta=1, 10)/ &
      0.97288_kind_phys, 0.97498_kind_phys, 0.97927_kind_phys, 0.98385_kind_phys, 0.98789_kind_phys, 0.99113_kind_phys, 0.99360_kind_phys, 0.99541_kind_phys, &
      0.99673_kind_phys, 0.99766_kind_phys/
    data(bm3i(9, 4, ibeta), ibeta=1, 10)/ &
      0.95403_kind_phys, 0.95741_kind_phys, 0.96440_kind_phys, 0.97202_kind_phys, 0.97887_kind_phys, 0.98444_kind_phys, 0.98872_kind_phys, 0.99190_kind_phys, &
      0.99421_kind_phys, 0.99586_kind_phys/
    data(bm3i(9, 5, ibeta), ibeta=1, 10)/ &
      0.91845_kind_phys, 0.92399_kind_phys, 0.93567_kind_phys, 0.94873_kind_phys, 0.96076_kind_phys, 0.97079_kind_phys, 0.97865_kind_phys, 0.98457_kind_phys, &
      0.98892_kind_phys, 0.99206_kind_phys/
    data(bm3i(9, 6, ibeta), ibeta=1, 10)/ &
      0.86762_kind_phys, 0.87533_kind_phys, 0.89202_kind_phys, 0.91148_kind_phys, 0.93027_kind_phys, 0.94669_kind_phys, 0.96013_kind_phys, 0.97062_kind_phys, &
      0.97855_kind_phys, 0.98441_kind_phys/
    data(bm3i(9, 7, ibeta), ibeta=1, 10)/ &
      0.84550_kind_phys, 0.85253_kind_phys, 0.86816_kind_phys, 0.88721_kind_phys, 0.90671_kind_phys, 0.92490_kind_phys, 0.94083_kind_phys, 0.95413_kind_phys, &
      0.96481_kind_phys, 0.97314_kind_phys/
    data(bm3i(9, 8, ibeta), ibeta=1, 10)/ &
      0.90138_kind_phys, 0.90544_kind_phys, 0.91437_kind_phys, 0.92513_kind_phys, 0.93602_kind_phys, 0.94615_kind_phys, 0.95506_kind_phys, 0.96258_kind_phys, &
      0.96868_kind_phys, 0.97347_kind_phys/
    data(bm3i(9, 9, ibeta), ibeta=1, 10)/ &
      0.96248_kind_phys, 0.96415_kind_phys, 0.96773_kind_phys, 0.97187_kind_phys, 0.97583_kind_phys, 0.97925_kind_phys, 0.98198_kind_phys, 0.98394_kind_phys, &
      0.98514_kind_phys, 0.98559_kind_phys/
    data(bm3i(9, 10, ibeta), ibeta=1, 10)/ &
      0.98837_kind_phys, 0.98892_kind_phys, 0.99005_kind_phys, 0.99127_kind_phys, 0.99232_kind_phys, 0.99306_kind_phys, 0.99339_kind_phys, 0.99328_kind_phys, &
      0.99269_kind_phys, 0.99161_kind_phys/
    data(bm3i(10, 1, ibeta), ibeta=1, 10)/ &
      0.99080_kind_phys, 0.99158_kind_phys, 0.99311_kind_phys, 0.99471_kind_phys, 0.99607_kind_phys, 0.99715_kind_phys, 0.99795_kind_phys, 0.99853_kind_phys, &
      0.99895_kind_phys, 0.99925_kind_phys/
    data(bm3i(10, 2, ibeta), ibeta=1, 10)/ &
      0.98910_kind_phys, 0.99001_kind_phys, 0.99182_kind_phys, 0.99371_kind_phys, 0.99533_kind_phys, 0.99661_kind_phys, 0.99757_kind_phys, 0.99826_kind_phys, &
      0.99876_kind_phys, 0.99912_kind_phys/
    data(bm3i(10, 3, ibeta), ibeta=1, 10)/ &
      0.98374_kind_phys, 0.98506_kind_phys, 0.98772_kind_phys, 0.99051_kind_phys, 0.99294_kind_phys, 0.99486_kind_phys, 0.99630_kind_phys, 0.99736_kind_phys, &
      0.99812_kind_phys, 0.99866_kind_phys/
    data(bm3i(10, 4, ibeta), ibeta=1, 10)/ &
      0.97238_kind_phys, 0.97453_kind_phys, 0.97892_kind_phys, 0.98361_kind_phys, 0.98773_kind_phys, 0.99104_kind_phys, 0.99354_kind_phys, 0.99538_kind_phys, &
      0.99671_kind_phys, 0.99765_kind_phys/
    data(bm3i(10, 5, ibeta), ibeta=1, 10)/ &
      0.94961_kind_phys, 0.95333_kind_phys, 0.96103_kind_phys, 0.96941_kind_phys, 0.97693_kind_phys, 0.98303_kind_phys, 0.98772_kind_phys, 0.99119_kind_phys, &
      0.99371_kind_phys, 0.99551_kind_phys/
    data(bm3i(10, 6, ibeta), ibeta=1, 10)/ &
      0.90943_kind_phys, 0.91550_kind_phys, 0.92834_kind_phys, 0.94275_kind_phys, 0.95608_kind_phys, 0.96723_kind_phys, 0.97600_kind_phys, 0.98263_kind_phys, &
      0.98751_kind_phys, 0.99103_kind_phys/
    data(bm3i(10, 7, ibeta), ibeta=1, 10)/ &
      0.86454_kind_phys, 0.87200_kind_phys, 0.88829_kind_phys, 0.90749_kind_phys, 0.92630_kind_phys, 0.94300_kind_phys, 0.95687_kind_phys, 0.96785_kind_phys, &
      0.97626_kind_phys, 0.98254_kind_phys/
    data(bm3i(10, 8, ibeta), ibeta=1, 10)/ &
      0.87498_kind_phys, 0.88048_kind_phys, 0.89264_kind_phys, 0.90737_kind_phys, 0.92240_kind_phys, 0.93642_kind_phys, 0.94877_kind_phys, 0.95917_kind_phys, &
      0.96762_kind_phys, 0.97429_kind_phys/
    data(bm3i(10, 9, ibeta), ibeta=1, 10)/ &
      0.93946_kind_phys, 0.94209_kind_phys, 0.94781_kind_phys, 0.95452_kind_phys, 0.96111_kind_phys, 0.96704_kind_phys, 0.97203_kind_phys, 0.97602_kind_phys, &
      0.97900_kind_phys, 0.98106_kind_phys/
    data(bm3i(10, 10, ibeta), ibeta=1, 10)/ &
      0.97977_kind_phys, 0.98071_kind_phys, 0.98270_kind_phys, 0.98492_kind_phys, 0.98695_kind_phys, 0.98858_kind_phys, 0.98970_kind_phys, 0.99027_kind_phys, &
      0.99026_kind_phys, 0.98968_kind_phys/

! fsb fm correction for intramodal m2 coagulation
    data bm2ii/ &
      0.707107_kind_phys, 0.720583_kind_phys, 0.745310_kind_phys, 0.748056_kind_phys, 0.696935_kind_phys, &
      0.604164_kind_phys, 0.504622_kind_phys, 0.416559_kind_phys, 0.343394_kind_phys, 0.283641_kind_phys/

! *** total correction for intramodal m2 coagulation

    data bm2iitt/ &
      1.000000_kind_phys, 0.907452_kind_phys, 0.680931_kind_phys, 0.409815_kind_phys, 0.196425_kind_phys, &
      0.078814_kind_phys, 0.028473_kind_phys, 0.009800_kind_phys, 0.003322_kind_phys, 0.001129_kind_phys/

! fsb fm correction for m2 i to j coagulation

    data(bm2ij(1, 1, ibeta), ibeta=1, 10)/ &
      0.707107_kind_phys, 0.716828_kind_phys, 0.738240_kind_phys, 0.764827_kind_phys, 0.793610_kind_phys, &
      0.822843_kind_phys, 0.851217_kind_phys, 0.877670_kind_phys, 0.901404_kind_phys, 0.921944_kind_phys/
    data(bm2ij(1, 2, ibeta), ibeta=1, 10)/ &
      0.719180_kind_phys, 0.727975_kind_phys, 0.747638_kind_phys, 0.772334_kind_phys, 0.799234_kind_phys, &
      0.826666_kind_phys, 0.853406_kind_phys, 0.878482_kind_phys, 0.901162_kind_phys, 0.920987_kind_phys/
    data(bm2ij(1, 3, ibeta), ibeta=1, 10)/ &
      0.760947_kind_phys, 0.767874_kind_phys, 0.783692_kind_phys, 0.803890_kind_phys, 0.826015_kind_phys, &
      0.848562_kind_phys, 0.870498_kind_phys, 0.891088_kind_phys, 0.909823_kind_phys, 0.926400_kind_phys/
    data(bm2ij(1, 4, ibeta), ibeta=1, 10)/ &
      0.830926_kind_phys, 0.836034_kind_phys, 0.847708_kind_phys, 0.862528_kind_phys, 0.878521_kind_phys, &
      0.894467_kind_phys, 0.909615_kind_phys, 0.923520_kind_phys, 0.935959_kind_phys, 0.946858_kind_phys/
    data(bm2ij(1, 5, ibeta), ibeta=1, 10)/ &
      0.903643_kind_phys, 0.907035_kind_phys, 0.914641_kind_phys, 0.924017_kind_phys, 0.933795_kind_phys, &
      0.943194_kind_phys, 0.951806_kind_phys, 0.959449_kind_phys, 0.966087_kind_phys, 0.971761_kind_phys/
    data(bm2ij(1, 6, ibeta), ibeta=1, 10)/ &
      0.954216_kind_phys, 0.956094_kind_phys, 0.960211_kind_phys, 0.965123_kind_phys, 0.970068_kind_phys, &
      0.974666_kind_phys, 0.978750_kind_phys, 0.982277_kind_phys, 0.985268_kind_phys, 0.987775_kind_phys/
    data(bm2ij(1, 7, ibeta), ibeta=1, 10)/ &
      0.980546_kind_phys, 0.981433_kind_phys, 0.983343_kind_phys, 0.985568_kind_phys, 0.987751_kind_phys, &
      0.989735_kind_phys, 0.991461_kind_phys, 0.992926_kind_phys, 0.994150_kind_phys, 0.995164_kind_phys/
    data(bm2ij(1, 8, ibeta), ibeta=1, 10)/ &
      0.992142_kind_phys, 0.992524_kind_phys, 0.993338_kind_phys, 0.994272_kind_phys, 0.995174_kind_phys, &
      0.995981_kind_phys, 0.996675_kind_phys, 0.997257_kind_phys, 0.997740_kind_phys, 0.998137_kind_phys/
    data(bm2ij(1, 9, ibeta), ibeta=1, 10)/ &
      0.996868_kind_phys, 0.997026_kind_phys, 0.997361_kind_phys, 0.997742_kind_phys, 0.998106_kind_phys, &
      0.998430_kind_phys, 0.998705_kind_phys, 0.998935_kind_phys, 0.999125_kind_phys, 0.999280_kind_phys/
    data(bm2ij(1, 10, ibeta), ibeta=1, 10)/ &
      0.998737_kind_phys, 0.998802_kind_phys, 0.998939_kind_phys, 0.999094_kind_phys, 0.999241_kind_phys, &
      0.999371_kind_phys, 0.999481_kind_phys, 0.999573_kind_phys, 0.999648_kind_phys, 0.999709_kind_phys/
    data(bm2ij(2, 1, ibeta), ibeta=1, 10)/ &
      0.729600_kind_phys, 0.739948_kind_phys, 0.763059_kind_phys, 0.791817_kind_phys, 0.822510_kind_phys, &
      0.852795_kind_phys, 0.881000_kind_phys, 0.905999_kind_phys, 0.927206_kind_phys, 0.944532_kind_phys/
    data(bm2ij(2, 2, ibeta), ibeta=1, 10)/ &
      0.727025_kind_phys, 0.737116_kind_phys, 0.759615_kind_phys, 0.787657_kind_phys, 0.817740_kind_phys, &
      0.847656_kind_phys, 0.875801_kind_phys, 0.901038_kind_phys, 0.922715_kind_phys, 0.940643_kind_phys/
    data(bm2ij(2, 3, ibeta), ibeta=1, 10)/ &
      0.738035_kind_phys, 0.746779_kind_phys, 0.766484_kind_phys, 0.791340_kind_phys, 0.818324_kind_phys, &
      0.845546_kind_phys, 0.871629_kind_phys, 0.895554_kind_phys, 0.916649_kind_phys, 0.934597_kind_phys/
    data(bm2ij(2, 4, ibeta), ibeta=1, 10)/ &
      0.784185_kind_phys, 0.790883_kind_phys, 0.806132_kind_phys, 0.825501_kind_phys, 0.846545_kind_phys, &
      0.867745_kind_phys, 0.888085_kind_phys, 0.906881_kind_phys, 0.923705_kind_phys, 0.938349_kind_phys/
    data(bm2ij(2, 5, ibeta), ibeta=1, 10)/ &
      0.857879_kind_phys, 0.862591_kind_phys, 0.873238_kind_phys, 0.886539_kind_phys, 0.900645_kind_phys, &
      0.914463_kind_phys, 0.927360_kind_phys, 0.939004_kind_phys, 0.949261_kind_phys, 0.958125_kind_phys/
    data(bm2ij(2, 6, ibeta), ibeta=1, 10)/ &
      0.925441_kind_phys, 0.928304_kind_phys, 0.934645_kind_phys, 0.942324_kind_phys, 0.950181_kind_phys, &
      0.957600_kind_phys, 0.964285_kind_phys, 0.970133_kind_phys, 0.975147_kind_phys, 0.979388_kind_phys/
    data(bm2ij(2, 7, ibeta), ibeta=1, 10)/ &
      0.966728_kind_phys, 0.968176_kind_phys, 0.971323_kind_phys, 0.975027_kind_phys, 0.978705_kind_phys, &
      0.982080_kind_phys, 0.985044_kind_phys, 0.987578_kind_phys, 0.989710_kind_phys, 0.991485_kind_phys/
    data(bm2ij(2, 8, ibeta), ibeta=1, 10)/ &
      0.986335_kind_phys, 0.986980_kind_phys, 0.988362_kind_phys, 0.989958_kind_phys, 0.991511_kind_phys, &
      0.992912_kind_phys, 0.994122_kind_phys, 0.995143_kind_phys, 0.995992_kind_phys, 0.996693_kind_phys/
    data(bm2ij(2, 9, ibeta), ibeta=1, 10)/ &
      0.994547_kind_phys, 0.994817_kind_phys, 0.995391_kind_phys, 0.996046_kind_phys, 0.996677_kind_phys, &
      0.997238_kind_phys, 0.997719_kind_phys, 0.998122_kind_phys, 0.998454_kind_phys, 0.998727_kind_phys/
    data(bm2ij(2, 10, ibeta), ibeta=1, 10)/ &
      0.997817_kind_phys, 0.997928_kind_phys, 0.998163_kind_phys, 0.998429_kind_phys, 0.998683_kind_phys, &
      0.998908_kind_phys, 0.999099_kind_phys, 0.999258_kind_phys, 0.999389_kind_phys, 0.999497_kind_phys/
    data(bm2ij(3, 1, ibeta), ibeta=1, 10)/ &
      0.783612_kind_phys, 0.793055_kind_phys, 0.814468_kind_phys, 0.841073_kind_phys, 0.868769_kind_phys, &
      0.894963_kind_phys, 0.918118_kind_phys, 0.937527_kind_phys, 0.953121_kind_phys, 0.965244_kind_phys/
    data(bm2ij(3, 2, ibeta), ibeta=1, 10)/ &
      0.772083_kind_phys, 0.781870_kind_phys, 0.803911_kind_phys, 0.831238_kind_phys, 0.859802_kind_phys, &
      0.887036_kind_phys, 0.911349_kind_phys, 0.931941_kind_phys, 0.948649_kind_phys, 0.961751_kind_phys/
    data(bm2ij(3, 3, ibeta), ibeta=1, 10)/ &
      0.755766_kind_phys, 0.765509_kind_phys, 0.787380_kind_phys, 0.814630_kind_phys, 0.843526_kind_phys, &
      0.871670_kind_phys, 0.897443_kind_phys, 0.919870_kind_phys, 0.938557_kind_phys, 0.953576_kind_phys/
    data(bm2ij(3, 4, ibeta), ibeta=1, 10)/ &
      0.763816_kind_phys, 0.772145_kind_phys, 0.790997_kind_phys, 0.814784_kind_phys, 0.840434_kind_phys, &
      0.865978_kind_phys, 0.890034_kind_phys, 0.911671_kind_phys, 0.930366_kind_phys, 0.945963_kind_phys/
    data(bm2ij(3, 5, ibeta), ibeta=1, 10)/ &
      0.813597_kind_phys, 0.819809_kind_phys, 0.833889_kind_phys, 0.851618_kind_phys, 0.870640_kind_phys, &
      0.889514_kind_phys, 0.907326_kind_phys, 0.923510_kind_phys, 0.937768_kind_phys, 0.950003_kind_phys/
    data(bm2ij(3, 6, ibeta), ibeta=1, 10)/ &
      0.886317_kind_phys, 0.890437_kind_phys, 0.899643_kind_phys, 0.910955_kind_phys, 0.922730_kind_phys, &
      0.934048_kind_phys, 0.944422_kind_phys, 0.953632_kind_phys, 0.961624_kind_phys, 0.968444_kind_phys/
    data(bm2ij(3, 7, ibeta), ibeta=1, 10)/ &
      0.944565_kind_phys, 0.946855_kind_phys, 0.951872_kind_phys, 0.957854_kind_phys, 0.963873_kind_phys, &
      0.969468_kind_phys, 0.974438_kind_phys, 0.978731_kind_phys, 0.982372_kind_phys, 0.985424_kind_phys/
    data(bm2ij(3, 8, ibeta), ibeta=1, 10)/ &
      0.976358_kind_phys, 0.977435_kind_phys, 0.979759_kind_phys, 0.982467_kind_phys, 0.985125_kind_phys, &
      0.987540_kind_phys, 0.989642_kind_phys, 0.991425_kind_phys, 0.992916_kind_phys, 0.994150_kind_phys/
    data(bm2ij(3, 9, ibeta), ibeta=1, 10)/ &
      0.990471_kind_phys, 0.990932_kind_phys, 0.991917_kind_phys, 0.993048_kind_phys, 0.994142_kind_phys, &
      0.995121_kind_phys, 0.995964_kind_phys, 0.996671_kind_phys, 0.997258_kind_phys, 0.997740_kind_phys/
    data(bm2ij(3, 10, ibeta), ibeta=1, 10)/ &
      0.996199_kind_phys, 0.996389_kind_phys, 0.996794_kind_phys, 0.997254_kind_phys, 0.997694_kind_phys, &
      0.998086_kind_phys, 0.998420_kind_phys, 0.998699_kind_phys, 0.998929_kind_phys, 0.999117_kind_phys/
    data(bm2ij(4, 1, ibeta), ibeta=1, 10)/ &
      0.844355_kind_phys, 0.852251_kind_phys, 0.869914_kind_phys, 0.891330_kind_phys, 0.912823_kind_phys, &
      0.932259_kind_phys, 0.948642_kind_phys, 0.961767_kind_phys, 0.971897_kind_phys, 0.979510_kind_phys/
    data(bm2ij(4, 2, ibeta), ibeta=1, 10)/ &
      0.831550_kind_phys, 0.839954_kind_phys, 0.858754_kind_phys, 0.881583_kind_phys, 0.904592_kind_phys, &
      0.925533_kind_phys, 0.943309_kind_phys, 0.957647_kind_phys, 0.968779_kind_phys, 0.977185_kind_phys/
    data(bm2ij(4, 3, ibeta), ibeta=1, 10)/ &
      0.803981_kind_phys, 0.813288_kind_phys, 0.834060_kind_phys, 0.859400_kind_phys, 0.885285_kind_phys, &
      0.909286_kind_phys, 0.930084_kind_phys, 0.947193_kind_phys, 0.960714_kind_phys, 0.971078_kind_phys/
    data(bm2ij(4, 4, ibeta), ibeta=1, 10)/ &
      0.781787_kind_phys, 0.791080_kind_phys, 0.811931_kind_phys, 0.837749_kind_phys, 0.864768_kind_phys, &
      0.890603_kind_phys, 0.913761_kind_phys, 0.933477_kind_phys, 0.949567_kind_phys, 0.962261_kind_phys/
    data(bm2ij(4, 5, ibeta), ibeta=1, 10)/ &
      0.791591_kind_phys, 0.799355_kind_phys, 0.816916_kind_phys, 0.838961_kind_phys, 0.862492_kind_phys, &
      0.885595_kind_phys, 0.907003_kind_phys, 0.925942_kind_phys, 0.942052_kind_phys, 0.955310_kind_phys/
    data(bm2ij(4, 6, ibeta), ibeta=1, 10)/ &
      0.844933_kind_phys, 0.850499_kind_phys, 0.863022_kind_phys, 0.878593_kind_phys, 0.895038_kind_phys, &
      0.911072_kind_phys, 0.925939_kind_phys, 0.939227_kind_phys, 0.950765_kind_phys, 0.960550_kind_phys/
    data(bm2ij(4, 7, ibeta), ibeta=1, 10)/ &
      0.912591_kind_phys, 0.916022_kind_phys, 0.923607_kind_phys, 0.932777_kind_phys, 0.942151_kind_phys, &
      0.951001_kind_phys, 0.958976_kind_phys, 0.965950_kind_phys, 0.971924_kind_phys, 0.976965_kind_phys/
    data(bm2ij(4, 8, ibeta), ibeta=1, 10)/ &
      0.959859_kind_phys, 0.961617_kind_phys, 0.965433_kind_phys, 0.969924_kind_phys, 0.974382_kind_phys, &
      0.978472_kind_phys, 0.982063_kind_phys, 0.985134_kind_phys, 0.987716_kind_phys, 0.989865_kind_phys/
    data(bm2ij(4, 9, ibeta), ibeta=1, 10)/ &
      0.983377_kind_phys, 0.984162_kind_phys, 0.985844_kind_phys, 0.987788_kind_phys, 0.989681_kind_phys, &
      0.991386_kind_phys, 0.992860_kind_phys, 0.994104_kind_phys, 0.995139_kind_phys, 0.995991_kind_phys/
    data(bm2ij(4, 10, ibeta), ibeta=1, 10)/ &
      0.993343_kind_phys, 0.993672_kind_phys, 0.994370_kind_phys, 0.995169_kind_phys, 0.995937_kind_phys, &
      0.996622_kind_phys, 0.997209_kind_phys, 0.997700_kind_phys, 0.998106_kind_phys, 0.998439_kind_phys/
    data(bm2ij(5, 1, ibeta), ibeta=1, 10)/ &
      0.895806_kind_phys, 0.901918_kind_phys, 0.915233_kind_phys, 0.930783_kind_phys, 0.945768_kind_phys, &
      0.958781_kind_phys, 0.969347_kind_phys, 0.977540_kind_phys, 0.983697_kind_phys, 0.988225_kind_phys/
    data(bm2ij(5, 2, ibeta), ibeta=1, 10)/ &
      0.885634_kind_phys, 0.892221_kind_phys, 0.906629_kind_phys, 0.923540_kind_phys, 0.939918_kind_phys, &
      0.954213_kind_phys, 0.965873_kind_phys, 0.974951_kind_phys, 0.981794_kind_phys, 0.986840_kind_phys/
    data(bm2ij(5, 3, ibeta), ibeta=1, 10)/ &
      0.860120_kind_phys, 0.867858_kind_phys, 0.884865_kind_phys, 0.904996_kind_phys, 0.924724_kind_phys, &
      0.942177_kind_phys, 0.956602_kind_phys, 0.967966_kind_phys, 0.976616_kind_phys, 0.983043_kind_phys/
    data(bm2ij(5, 4, ibeta), ibeta=1, 10)/ &
      0.827462_kind_phys, 0.836317_kind_phys, 0.855885_kind_phys, 0.879377_kind_phys, 0.902897_kind_phys, &
      0.924232_kind_phys, 0.942318_kind_phys, 0.956900_kind_phys, 0.968222_kind_phys, 0.976774_kind_phys/
    data(bm2ij(5, 5, ibeta), ibeta=1, 10)/ &
      0.805527_kind_phys, 0.814279_kind_phys, 0.833853_kind_phys, 0.857892_kind_phys, 0.882726_kind_phys, &
      0.906095_kind_phys, 0.926690_kind_phys, 0.943938_kind_phys, 0.957808_kind_phys, 0.968615_kind_phys/
    data(bm2ij(5, 6, ibeta), ibeta=1, 10)/ &
      0.820143_kind_phys, 0.827223_kind_phys, 0.843166_kind_phys, 0.863002_kind_phys, 0.883905_kind_phys, &
      0.904128_kind_phys, 0.922585_kind_phys, 0.938687_kind_phys, 0.952222_kind_phys, 0.963255_kind_phys/
    data(bm2ij(5, 7, ibeta), ibeta=1, 10)/ &
      0.875399_kind_phys, 0.880208_kind_phys, 0.890929_kind_phys, 0.904065_kind_phys, 0.917699_kind_phys, &
      0.930756_kind_phys, 0.942656_kind_phys, 0.953131_kind_phys, 0.962113_kind_phys, 0.969657_kind_phys/
    data(bm2ij(5, 8, ibeta), ibeta=1, 10)/ &
      0.934782_kind_phys, 0.937520_kind_phys, 0.943515_kind_phys, 0.950656_kind_phys, 0.957840_kind_phys, &
      0.964516_kind_phys, 0.970446_kind_phys, 0.975566_kind_phys, 0.979905_kind_phys, 0.983534_kind_phys/
    data(bm2ij(5, 9, ibeta), ibeta=1, 10)/ &
      0.971369_kind_phys, 0.972679_kind_phys, 0.975505_kind_phys, 0.978797_kind_phys, 0.982029_kind_phys, &
      0.984964_kind_phys, 0.987518_kind_phys, 0.989685_kind_phys, 0.991496_kind_phys, 0.992994_kind_phys/
    data(bm2ij(5, 10, ibeta), ibeta=1, 10)/ &
      0.988329_kind_phys, 0.988893_kind_phys, 0.990099_kind_phys, 0.991485_kind_phys, 0.992825_kind_phys, &
      0.994025_kind_phys, 0.995058_kind_phys, 0.995925_kind_phys, 0.996643_kind_phys, 0.997234_kind_phys/
    data(bm2ij(6, 1, ibeta), ibeta=1, 10)/ &
      0.933384_kind_phys, 0.937784_kind_phys, 0.947130_kind_phys, 0.957655_kind_phys, 0.967430_kind_phys, &
      0.975639_kind_phys, 0.982119_kind_phys, 0.987031_kind_phys, 0.990657_kind_phys, 0.993288_kind_phys/
    data(bm2ij(6, 2, ibeta), ibeta=1, 10)/ &
      0.926445_kind_phys, 0.931227_kind_phys, 0.941426_kind_phys, 0.952975_kind_phys, 0.963754_kind_phys, &
      0.972845_kind_phys, 0.980044_kind_phys, 0.985514_kind_phys, 0.989558_kind_phys, 0.992498_kind_phys/
    data(bm2ij(6, 3, ibeta), ibeta=1, 10)/ &
      0.907835_kind_phys, 0.913621_kind_phys, 0.926064_kind_phys, 0.940308_kind_phys, 0.953745_kind_phys, &
      0.965189_kind_phys, 0.974327_kind_phys, 0.981316_kind_phys, 0.986510_kind_phys, 0.990297_kind_phys/
    data(bm2ij(6, 4, ibeta), ibeta=1, 10)/ &
      0.879088_kind_phys, 0.886306_kind_phys, 0.901945_kind_phys, 0.920079_kind_phys, 0.937460_kind_phys, &
      0.952509_kind_phys, 0.964711_kind_phys, 0.974166_kind_phys, 0.981265_kind_phys, 0.986484_kind_phys/
    data(bm2ij(6, 5, ibeta), ibeta=1, 10)/ &
      0.846500_kind_phys, 0.854862_kind_phys, 0.873189_kind_phys, 0.894891_kind_phys, 0.916264_kind_phys, &
      0.935315_kind_phys, 0.951197_kind_phys, 0.963812_kind_phys, 0.973484_kind_phys, 0.980715_kind_phys/
    data(bm2ij(6, 6, ibeta), ibeta=1, 10)/ &
      0.828137_kind_phys, 0.836250_kind_phys, 0.854310_kind_phys, 0.876287_kind_phys, 0.898710_kind_phys, &
      0.919518_kind_phys, 0.937603_kind_phys, 0.952560_kind_phys, 0.964461_kind_phys, 0.973656_kind_phys/
    data(bm2ij(6, 7, ibeta), ibeta=1, 10)/ &
      0.848595_kind_phys, 0.854886_kind_phys, 0.868957_kind_phys, 0.886262_kind_phys, 0.904241_kind_phys, &
      0.921376_kind_phys, 0.936799_kind_phys, 0.950096_kind_phys, 0.961172_kind_phys, 0.970145_kind_phys/
    data(bm2ij(6, 8, ibeta), ibeta=1, 10)/ &
      0.902919_kind_phys, 0.906922_kind_phys, 0.915760_kind_phys, 0.926427_kind_phys, 0.937312_kind_phys, &
      0.947561_kind_phys, 0.956758_kind_phys, 0.964747_kind_phys, 0.971525_kind_phys, 0.977175_kind_phys/
    data(bm2ij(6, 9, ibeta), ibeta=1, 10)/ &
      0.952320_kind_phys, 0.954434_kind_phys, 0.959021_kind_phys, 0.964418_kind_phys, 0.969774_kind_phys, &
      0.974688_kind_phys, 0.979003_kind_phys, 0.982690_kind_phys, 0.985789_kind_phys, 0.988364_kind_phys/
    data(bm2ij(6, 10, ibeta), ibeta=1, 10)/ &
      0.979689_kind_phys, 0.980650_kind_phys, 0.982712_kind_phys, 0.985093_kind_phys, 0.987413_kind_phys, &
      0.989502_kind_phys, 0.991308_kind_phys, 0.992831_kind_phys, 0.994098_kind_phys, 0.995142_kind_phys/
    data(bm2ij(7, 1, ibeta), ibeta=1, 10)/ &
      0.958611_kind_phys, 0.961598_kind_phys, 0.967817_kind_phys, 0.974620_kind_phys, 0.980752_kind_phys, &
      0.985771_kind_phys, 0.989650_kind_phys, 0.992543_kind_phys, 0.994653_kind_phys, 0.996171_kind_phys/
    data(bm2ij(7, 2, ibeta), ibeta=1, 10)/ &
      0.954225_kind_phys, 0.957488_kind_phys, 0.964305_kind_phys, 0.971795_kind_phys, 0.978576_kind_phys, &
      0.984144_kind_phys, 0.988458_kind_phys, 0.991681_kind_phys, 0.994034_kind_phys, 0.995728_kind_phys/
    data(bm2ij(7, 3, ibeta), ibeta=1, 10)/ &
      0.942147_kind_phys, 0.946158_kind_phys, 0.954599_kind_phys, 0.963967_kind_phys, 0.972529_kind_phys, &
      0.979612_kind_phys, 0.985131_kind_phys, 0.989271_kind_phys, 0.992301_kind_phys, 0.994487_kind_phys/
    data(bm2ij(7, 4, ibeta), ibeta=1, 10)/ &
      0.921821_kind_phys, 0.927048_kind_phys, 0.938140_kind_phys, 0.950598_kind_phys, 0.962118_kind_phys, &
      0.971752_kind_phys, 0.979326_kind_phys, 0.985046_kind_phys, 0.989254_kind_phys, 0.992299_kind_phys/
    data(bm2ij(7, 5, ibeta), ibeta=1, 10)/ &
      0.893419_kind_phys, 0.900158_kind_phys, 0.914598_kind_phys, 0.931070_kind_phys, 0.946584_kind_phys, &
      0.959795_kind_phys, 0.970350_kind_phys, 0.978427_kind_phys, 0.984432_kind_phys, 0.988811_kind_phys/
    data(bm2ij(7, 6, ibeta), ibeta=1, 10)/ &
      0.863302_kind_phys, 0.871111_kind_phys, 0.888103_kind_phys, 0.907990_kind_phys, 0.927305_kind_phys, &
      0.944279_kind_phys, 0.958245_kind_phys, 0.969211_kind_phys, 0.977540_kind_phys, 0.983720_kind_phys/
    data(bm2ij(7, 7, ibeta), ibeta=1, 10)/ &
      0.850182_kind_phys, 0.857560_kind_phys, 0.873890_kind_phys, 0.893568_kind_phys, 0.913408_kind_phys, &
      0.931591_kind_phys, 0.947216_kind_phys, 0.960014_kind_phys, 0.970121_kind_phys, 0.977886_kind_phys/
    data(bm2ij(7, 8, ibeta), ibeta=1, 10)/ &
      0.875837_kind_phys, 0.881265_kind_phys, 0.893310_kind_phys, 0.907936_kind_phys, 0.922910_kind_phys, &
      0.936977_kind_phys, 0.949480_kind_phys, 0.960154_kind_phys, 0.968985_kind_phys, 0.976111_kind_phys/
    data(bm2ij(7, 9, ibeta), ibeta=1, 10)/ &
      0.926228_kind_phys, 0.929445_kind_phys, 0.936486_kind_phys, 0.944868_kind_phys, 0.953293_kind_phys, &
      0.961108_kind_phys, 0.968028_kind_phys, 0.973973_kind_phys, 0.978974_kind_phys, 0.983118_kind_phys/
    data(bm2ij(7, 10, ibeta), ibeta=1, 10)/ &
      0.965533_kind_phys, 0.967125_kind_phys, 0.970558_kind_phys, 0.974557_kind_phys, 0.978484_kind_phys, &
      0.982050_kind_phys, 0.985153_kind_phys, 0.987785_kind_phys, 0.989982_kind_phys, 0.991798_kind_phys/
    data(bm2ij(8, 1, ibeta), ibeta=1, 10)/ &
      0.974731_kind_phys, 0.976674_kind_phys, 0.980660_kind_phys, 0.984926_kind_phys, 0.988689_kind_phys, &
      0.991710_kind_phys, 0.994009_kind_phys, 0.995703_kind_phys, 0.996929_kind_phys, 0.997805_kind_phys/
    data(bm2ij(8, 2, ibeta), ibeta=1, 10)/ &
      0.972062_kind_phys, 0.974192_kind_phys, 0.978571_kind_phys, 0.983273_kind_phys, 0.987432_kind_phys, &
      0.990780_kind_phys, 0.993333_kind_phys, 0.995218_kind_phys, 0.996581_kind_phys, 0.997557_kind_phys/
    data(bm2ij(8, 3, ibeta), ibeta=1, 10)/ &
      0.964662_kind_phys, 0.967300_kind_phys, 0.972755_kind_phys, 0.978659_kind_phys, 0.983921_kind_phys, &
      0.988181_kind_phys, 0.991444_kind_phys, 0.993859_kind_phys, 0.995610_kind_phys, 0.996863_kind_phys/
    data(bm2ij(8, 4, ibeta), ibeta=1, 10)/ &
      0.951782_kind_phys, 0.955284_kind_phys, 0.962581_kind_phys, 0.970559_kind_phys, 0.977737_kind_phys, &
      0.983593_kind_phys, 0.988103_kind_phys, 0.991454_kind_phys, 0.993889_kind_phys, 0.995635_kind_phys/
    data(bm2ij(8, 5, ibeta), ibeta=1, 10)/ &
      0.931947_kind_phys, 0.936723_kind_phys, 0.946751_kind_phys, 0.957843_kind_phys, 0.967942_kind_phys, &
      0.976267_kind_phys, 0.982734_kind_phys, 0.987571_kind_phys, 0.991102_kind_phys, 0.993642_kind_phys/
    data(bm2ij(8, 6, ibeta), ibeta=1, 10)/ &
      0.905410_kind_phys, 0.911665_kind_phys, 0.924950_kind_phys, 0.939908_kind_phys, 0.953798_kind_phys, &
      0.965469_kind_phys, 0.974684_kind_phys, 0.981669_kind_phys, 0.986821_kind_phys, 0.990556_kind_phys/
    data(bm2ij(8, 7, ibeta), ibeta=1, 10)/ &
      0.878941_kind_phys, 0.886132_kind_phys, 0.901679_kind_phys, 0.919688_kind_phys, 0.936970_kind_phys, &
      0.951980_kind_phys, 0.964199_kind_phys, 0.973709_kind_phys, 0.980881_kind_phys, 0.986174_kind_phys/
    data(bm2ij(8, 8, ibeta), ibeta=1, 10)/ &
      0.871653_kind_phys, 0.878218_kind_phys, 0.892652_kind_phys, 0.909871_kind_phys, 0.927034_kind_phys, &
      0.942592_kind_phys, 0.955836_kind_phys, 0.966604_kind_phys, 0.975065_kind_phys, 0.981545_kind_phys/
    data(bm2ij(8, 9, ibeta), ibeta=1, 10)/ &
      0.900693_kind_phys, 0.905239_kind_phys, 0.915242_kind_phys, 0.927232_kind_phys, 0.939335_kind_phys, &
      0.950555_kind_phys, 0.960420_kind_phys, 0.968774_kind_phys, 0.975651_kind_phys, 0.981188_kind_phys/
    data(bm2ij(8, 10, ibeta), ibeta=1, 10)/ &
      0.944922_kind_phys, 0.947435_kind_phys, 0.952894_kind_phys, 0.959317_kind_phys, 0.965689_kind_phys, &
      0.971529_kind_phys, 0.976645_kind_phys, 0.981001_kind_phys, 0.984641_kind_phys, 0.987642_kind_phys/
    data(bm2ij(9, 1, ibeta), ibeta=1, 10)/ &
      0.984736_kind_phys, 0.985963_kind_phys, 0.988453_kind_phys, 0.991078_kind_phys, 0.993357_kind_phys, &
      0.995161_kind_phys, 0.996519_kind_phys, 0.997512_kind_phys, 0.998226_kind_phys, 0.998734_kind_phys/
    data(bm2ij(9, 2, ibeta), ibeta=1, 10)/ &
      0.983141_kind_phys, 0.984488_kind_phys, 0.987227_kind_phys, 0.990119_kind_phys, 0.992636_kind_phys, &
      0.994632_kind_phys, 0.996137_kind_phys, 0.997238_kind_phys, 0.998030_kind_phys, 0.998595_kind_phys/
    data(bm2ij(9, 3, ibeta), ibeta=1, 10)/ &
      0.978726_kind_phys, 0.980401_kind_phys, 0.983819_kind_phys, 0.987450_kind_phys, 0.990626_kind_phys, &
      0.993157_kind_phys, 0.995071_kind_phys, 0.996475_kind_phys, 0.997486_kind_phys, 0.998206_kind_phys/
    data(bm2ij(9, 4, ibeta), ibeta=1, 10)/ &
      0.970986_kind_phys, 0.973224_kind_phys, 0.977818_kind_phys, 0.982737_kind_phys, 0.987072_kind_phys, &
      0.990546_kind_phys, 0.993184_kind_phys, 0.995124_kind_phys, 0.996523_kind_phys, 0.997521_kind_phys/
    data(bm2ij(9, 5, ibeta), ibeta=1, 10)/ &
      0.958579_kind_phys, 0.961700_kind_phys, 0.968149_kind_phys, 0.975116_kind_phys, 0.981307_kind_phys, &
      0.986301_kind_phys, 0.990112_kind_phys, 0.992923_kind_phys, 0.994954_kind_phys, 0.996404_kind_phys/
    data(bm2ij(9, 6, ibeta), ibeta=1, 10)/ &
      0.940111_kind_phys, 0.944479_kind_phys, 0.953572_kind_phys, 0.963506_kind_phys, 0.972436_kind_phys, &
      0.979714_kind_phys, 0.985313_kind_phys, 0.989468_kind_phys, 0.992483_kind_phys, 0.994641_kind_phys/
    data(bm2ij(9, 7, ibeta), ibeta=1, 10)/ &
      0.916127_kind_phys, 0.921878_kind_phys, 0.934003_kind_phys, 0.947506_kind_phys, 0.959899_kind_phys, &
      0.970199_kind_phys, 0.978255_kind_phys, 0.984314_kind_phys, 0.988755_kind_phys, 0.991960_kind_phys/
    data(bm2ij(9, 8, ibeta), ibeta=1, 10)/ &
      0.893848_kind_phys, 0.900364_kind_phys, 0.914368_kind_phys, 0.930438_kind_phys, 0.945700_kind_phys, &
      0.958824_kind_phys, 0.969416_kind_phys, 0.977603_kind_phys, 0.983746_kind_phys, 0.988262_kind_phys/
    data(bm2ij(9, 9, ibeta), ibeta=1, 10)/ &
      0.892161_kind_phys, 0.897863_kind_phys, 0.910315_kind_phys, 0.925021_kind_phys, 0.939523_kind_phys, &
      0.952544_kind_phys, 0.963544_kind_phys, 0.972442_kind_phys, 0.979411_kind_phys, 0.984742_kind_phys/
    data(bm2ij(9, 10, ibeta), ibeta=1, 10)/ &
      0.922260_kind_phys, 0.925966_kind_phys, 0.934047_kind_phys, 0.943616_kind_phys, 0.953152_kind_phys, &
      0.961893_kind_phys, 0.969506_kind_phys, 0.975912_kind_phys, 0.981167_kind_phys, 0.985394_kind_phys/
    data(bm2ij(10, 1, ibeta), ibeta=1, 10)/ &
      0.990838_kind_phys, 0.991598_kind_phys, 0.993128_kind_phys, 0.994723_kind_phys, 0.996092_kind_phys, &
      0.997167_kind_phys, 0.997969_kind_phys, 0.998552_kind_phys, 0.998969_kind_phys, 0.999265_kind_phys/
    data(bm2ij(10, 2, ibeta), ibeta=1, 10)/ &
      0.989892_kind_phys, 0.990727_kind_phys, 0.992411_kind_phys, 0.994167_kind_phys, 0.995678_kind_phys, &
      0.996864_kind_phys, 0.997751_kind_phys, 0.998396_kind_phys, 0.998858_kind_phys, 0.999186_kind_phys/
    data(bm2ij(10, 3, ibeta), ibeta=1, 10)/ &
      0.987287_kind_phys, 0.988327_kind_phys, 0.990428_kind_phys, 0.992629_kind_phys, 0.994529_kind_phys, &
      0.996026_kind_phys, 0.997148_kind_phys, 0.997965_kind_phys, 0.998551_kind_phys, 0.998967_kind_phys/
    data(bm2ij(10, 4, ibeta), ibeta=1, 10)/ &
      0.982740_kind_phys, 0.984130_kind_phys, 0.986952_kind_phys, 0.989926_kind_phys, 0.992508_kind_phys, &
      0.994551_kind_phys, 0.996087_kind_phys, 0.997208_kind_phys, 0.998012_kind_phys, 0.998584_kind_phys/
    data(bm2ij(10, 5, ibeta), ibeta=1, 10)/ &
      0.975380_kind_phys, 0.977330_kind_phys, 0.981307_kind_phys, 0.985529_kind_phys, 0.989216_kind_phys, &
      0.992147_kind_phys, 0.994358_kind_phys, 0.995975_kind_phys, 0.997136_kind_phys, 0.997961_kind_phys/
    data(bm2ij(10, 6, ibeta), ibeta=1, 10)/ &
      0.963911_kind_phys, 0.966714_kind_phys, 0.972465_kind_phys, 0.978614_kind_phys, 0.984022_kind_phys, &
      0.988346_kind_phys, 0.991620_kind_phys, 0.994020_kind_phys, 0.995747_kind_phys, 0.996974_kind_phys/
    data(bm2ij(10, 7, ibeta), ibeta=1, 10)/ &
      0.947187_kind_phys, 0.951161_kind_phys, 0.959375_kind_phys, 0.968258_kind_phys, 0.976160_kind_phys, &
      0.982540_kind_phys, 0.987409_kind_phys, 0.991000_kind_phys, 0.993592_kind_phys, 0.995441_kind_phys/
    data(bm2ij(10, 8, ibeta), ibeta=1, 10)/ &
      0.926045_kind_phys, 0.931270_kind_phys, 0.942218_kind_phys, 0.954297_kind_phys, 0.965273_kind_phys, &
      0.974311_kind_phys, 0.981326_kind_phys, 0.986569_kind_phys, 0.990394_kind_phys, 0.993143_kind_phys/
    data(bm2ij(10, 9, ibeta), ibeta=1, 10)/ &
      0.908092_kind_phys, 0.913891_kind_phys, 0.926288_kind_phys, 0.940393_kind_phys, 0.953667_kind_phys, &
      0.964987_kind_phys, 0.974061_kind_phys, 0.981038_kind_phys, 0.986253_kind_phys, 0.990078_kind_phys/
    data(bm2ij(10, 10, ibeta), ibeta=1, 10)/ &
      0.911143_kind_phys, 0.915972_kind_phys, 0.926455_kind_phys, 0.938721_kind_phys, 0.950701_kind_phys, &
      0.961370_kind_phys, 0.970329_kind_phys, 0.977549_kind_phys, 0.983197_kind_phys, 0.987518_kind_phys/

! fsb total correction factor for m2 coagulation j from i

    data(bm2ji(1, 1, ibeta), ibeta=1, 10)/ &
      0.753466_kind_phys, 0.756888_kind_phys, 0.761008_kind_phys, 0.759432_kind_phys, 0.748675_kind_phys, &
      0.726951_kind_phys, 0.693964_kind_phys, 0.650915_kind_phys, 0.600227_kind_phys, 0.545000_kind_phys/
    data(bm2ji(1, 2, ibeta), ibeta=1, 10)/ &
      0.824078_kind_phys, 0.828698_kind_phys, 0.835988_kind_phys, 0.838943_kind_phys, 0.833454_kind_phys, &
      0.817148_kind_phys, 0.789149_kind_phys, 0.750088_kind_phys, 0.701887_kind_phys, 0.647308_kind_phys/
    data(bm2ji(1, 3, ibeta), ibeta=1, 10)/ &
      1.007389_kind_phys, 1.014362_kind_phys, 1.028151_kind_phys, 1.041011_kind_phys, 1.047939_kind_phys, &
      1.045707_kind_phys, 1.032524_kind_phys, 1.007903_kind_phys, 0.972463_kind_phys, 0.927667_kind_phys/
    data(bm2ji(1, 4, ibeta), ibeta=1, 10)/ &
      1.246157_kind_phys, 1.255135_kind_phys, 1.274249_kind_phys, 1.295351_kind_phys, 1.313362_kind_phys, &
      1.325187_kind_phys, 1.329136_kind_phys, 1.324491_kind_phys, 1.311164_kind_phys, 1.289459_kind_phys/
    data(bm2ji(1, 5, ibeta), ibeta=1, 10)/ &
      1.450823_kind_phys, 1.459551_kind_phys, 1.478182_kind_phys, 1.499143_kind_phys, 1.518224_kind_phys, &
      1.533312_kind_phys, 1.543577_kind_phys, 1.548882_kind_phys, 1.549395_kind_phys, 1.545364_kind_phys/
    data(bm2ji(1, 6, ibeta), ibeta=1, 10)/ &
      1.575248_kind_phys, 1.581832_kind_phys, 1.595643_kind_phys, 1.610866_kind_phys, 1.624601_kind_phys, &
      1.635690_kind_phys, 1.643913_kind_phys, 1.649470_kind_phys, 1.652688_kind_phys, 1.653878_kind_phys/
    data(bm2ji(1, 7, ibeta), ibeta=1, 10)/ &
      1.638426_kind_phys, 1.642626_kind_phys, 1.651293_kind_phys, 1.660641_kind_phys, 1.668926_kind_phys, &
      1.675571_kind_phys, 1.680572_kind_phys, 1.684147_kind_phys, 1.686561_kind_phys, 1.688047_kind_phys/
    data(bm2ji(1, 8, ibeta), ibeta=1, 10)/ &
      1.669996_kind_phys, 1.672392_kind_phys, 1.677283_kind_phys, 1.682480_kind_phys, 1.687028_kind_phys, &
      1.690651_kind_phys, 1.693384_kind_phys, 1.695372_kind_phys, 1.696776_kind_phys, 1.697734_kind_phys/
    data(bm2ji(1, 9, ibeta), ibeta=1, 10)/ &
      1.686148_kind_phys, 1.687419_kind_phys, 1.689993_kind_phys, 1.692704_kind_phys, 1.695057_kind_phys, &
      1.696922_kind_phys, 1.698329_kind_phys, 1.699359_kind_phys, 1.700099_kind_phys, 1.700621_kind_phys/
    data(bm2ji(1, 10, ibeta), ibeta=1, 10)/ &
      1.694364_kind_phys, 1.695010_kind_phys, 1.696313_kind_phys, 1.697676_kind_phys, 1.698853_kind_phys, &
      1.699782_kind_phys, 1.700482_kind_phys, 1.700996_kind_phys, 1.701366_kind_phys, 1.701631_kind_phys/
    data(bm2ji(2, 1, ibeta), ibeta=1, 10)/ &
      0.783166_kind_phys, 0.779369_kind_phys, 0.768044_kind_phys, 0.747572_kind_phys, 0.716709_kind_phys, &
      0.675422_kind_phys, 0.624981_kind_phys, 0.567811_kind_phys, 0.507057_kind_phys, 0.445975_kind_phys/
    data(bm2ji(2, 2, ibeta), ibeta=1, 10)/ &
      0.848390_kind_phys, 0.847100_kind_phys, 0.840874_kind_phys, 0.826065_kind_phys, 0.800296_kind_phys, &
      0.762625_kind_phys, 0.713655_kind_phys, 0.655545_kind_phys, 0.591603_kind_phys, 0.525571_kind_phys/
    data(bm2ji(2, 3, ibeta), ibeta=1, 10)/ &
      1.039894_kind_phys, 1.043786_kind_phys, 1.049445_kind_phys, 1.049664_kind_phys, 1.039407_kind_phys, &
      1.015322_kind_phys, 0.975983_kind_phys, 0.922180_kind_phys, 0.856713_kind_phys, 0.783634_kind_phys/
    data(bm2ji(2, 4, ibeta), ibeta=1, 10)/ &
      1.345995_kind_phys, 1.356064_kind_phys, 1.376947_kind_phys, 1.398304_kind_phys, 1.412685_kind_phys, &
      1.414611_kind_phys, 1.400652_kind_phys, 1.369595_kind_phys, 1.322261_kind_phys, 1.260993_kind_phys/
    data(bm2ji(2, 5, ibeta), ibeta=1, 10)/ &
      1.675575_kind_phys, 1.689859_kind_phys, 1.720957_kind_phys, 1.756659_kind_phys, 1.788976_kind_phys, &
      1.812679_kind_phys, 1.824773_kind_phys, 1.824024_kind_phys, 1.810412_kind_phys, 1.784630_kind_phys/
    data(bm2ji(2, 6, ibeta), ibeta=1, 10)/ &
      1.919835_kind_phys, 1.933483_kind_phys, 1.962973_kind_phys, 1.996810_kind_phys, 2.028377_kind_phys, &
      2.054172_kind_phys, 2.072763_kind_phys, 2.083963_kind_phys, 2.088190_kind_phys, 2.086052_kind_phys/
    data(bm2ji(2, 7, ibeta), ibeta=1, 10)/ &
      2.064139_kind_phys, 2.074105_kind_phys, 2.095233_kind_phys, 2.118909_kind_phys, 2.140688_kind_phys, &
      2.158661_kind_phys, 2.172373_kind_phys, 2.182087_kind_phys, 2.188330_kind_phys, 2.191650_kind_phys/
    data(bm2ji(2, 8, ibeta), ibeta=1, 10)/ &
      2.144871_kind_phys, 2.150990_kind_phys, 2.163748_kind_phys, 2.177731_kind_phys, 2.190364_kind_phys, &
      2.200712_kind_phys, 2.208687_kind_phys, 2.214563_kind_phys, 2.218716_kind_phys, 2.221502_kind_phys/
    data(bm2ji(2, 9, ibeta), ibeta=1, 10)/ &
      2.189223_kind_phys, 2.192595_kind_phys, 2.199540_kind_phys, 2.207033_kind_phys, 2.213706_kind_phys, &
      2.219125_kind_phys, 2.223297_kind_phys, 2.226403_kind_phys, 2.228660_kind_phys, 2.230265_kind_phys/
    data(bm2ji(2, 10, ibeta), ibeta=1, 10)/ &
      2.212595_kind_phys, 2.214342_kind_phys, 2.217912_kind_phys, 2.221723_kind_phys, 2.225082_kind_phys, &
      2.227791_kind_phys, 2.229869_kind_phys, 2.231417_kind_phys, 2.232551_kind_phys, 2.233372_kind_phys/
    data(bm2ji(3, 1, ibeta), ibeta=1, 10)/ &
      0.837870_kind_phys, 0.824476_kind_phys, 0.793119_kind_phys, 0.750739_kind_phys, 0.700950_kind_phys, &
      0.646691_kind_phys, 0.590508_kind_phys, 0.534354_kind_phys, 0.479532_kind_phys, 0.426856_kind_phys/
    data(bm2ji(3, 2, ibeta), ibeta=1, 10)/ &
      0.896771_kind_phys, 0.885847_kind_phys, 0.859327_kind_phys, 0.821694_kind_phys, 0.775312_kind_phys, &
      0.722402_kind_phys, 0.665196_kind_phys, 0.605731_kind_phys, 0.545742_kind_phys, 0.486687_kind_phys/
    data(bm2ji(3, 3, ibeta), ibeta=1, 10)/ &
      1.076089_kind_phys, 1.071727_kind_phys, 1.058845_kind_phys, 1.036171_kind_phys, 1.002539_kind_phys, &
      0.957521_kind_phys, 0.901640_kind_phys, 0.836481_kind_phys, 0.764597_kind_phys, 0.689151_kind_phys/
    data(bm2ji(3, 4, ibeta), ibeta=1, 10)/ &
      1.409571_kind_phys, 1.415168_kind_phys, 1.425346_kind_phys, 1.432021_kind_phys, 1.428632_kind_phys, &
      1.409696_kind_phys, 1.371485_kind_phys, 1.312958_kind_phys, 1.236092_kind_phys, 1.145293_kind_phys/
    data(bm2ji(3, 5, ibeta), ibeta=1, 10)/ &
      1.862757_kind_phys, 1.880031_kind_phys, 1.918394_kind_phys, 1.963456_kind_phys, 2.004070_kind_phys, &
      2.030730_kind_phys, 2.036144_kind_phys, 2.016159_kind_phys, 1.970059_kind_phys, 1.900079_kind_phys/
    data(bm2ji(3, 6, ibeta), ibeta=1, 10)/ &
      2.289741_kind_phys, 2.313465_kind_phys, 2.366789_kind_phys, 2.431612_kind_phys, 2.495597_kind_phys, &
      2.549838_kind_phys, 2.588523_kind_phys, 2.608665_kind_phys, 2.609488_kind_phys, 2.591662_kind_phys/
    data(bm2ji(3, 7, ibeta), ibeta=1, 10)/ &
      2.597157_kind_phys, 2.618731_kind_phys, 2.666255_kind_phys, 2.722597_kind_phys, 2.777531_kind_phys, &
      2.825187_kind_phys, 2.862794_kind_phys, 2.889648_kind_phys, 2.906199_kind_phys, 2.913380_kind_phys/
    data(bm2ji(3, 8, ibeta), ibeta=1, 10)/ &
      2.797975_kind_phys, 2.813116_kind_phys, 2.845666_kind_phys, 2.882976_kind_phys, 2.918289_kind_phys, &
      2.948461_kind_phys, 2.972524_kind_phys, 2.990687_kind_phys, 3.003664_kind_phys, 3.012284_kind_phys/
    data(bm2ji(3, 9, ibeta), ibeta=1, 10)/ &
      2.920832_kind_phys, 2.929843_kind_phys, 2.948848_kind_phys, 2.970057_kind_phys, 2.989632_kind_phys, &
      3.006057_kind_phys, 3.019067_kind_phys, 3.028979_kind_phys, 3.036307_kind_phys, 3.041574_kind_phys/
    data(bm2ji(3, 10, ibeta), ibeta=1, 10)/ &
      2.989627_kind_phys, 2.994491_kind_phys, 3.004620_kind_phys, 3.015720_kind_phys, 3.025789_kind_phys, &
      3.034121_kind_phys, 3.040664_kind_phys, 3.045641_kind_phys, 3.049347_kind_phys, 3.052066_kind_phys/
    data(bm2ji(4, 1, ibeta), ibeta=1, 10)/ &
      0.893179_kind_phys, 0.870897_kind_phys, 0.820996_kind_phys, 0.759486_kind_phys, 0.695488_kind_phys, &
      0.634582_kind_phys, 0.579818_kind_phys, 0.532143_kind_phys, 0.490927_kind_phys, 0.454618_kind_phys/
    data(bm2ji(4, 2, ibeta), ibeta=1, 10)/ &
      0.948355_kind_phys, 0.927427_kind_phys, 0.880215_kind_phys, 0.821146_kind_phys, 0.758524_kind_phys, &
      0.697680_kind_phys, 0.641689_kind_phys, 0.591605_kind_phys, 0.546919_kind_phys, 0.506208_kind_phys/
    data(bm2ji(4, 3, ibeta), ibeta=1, 10)/ &
      1.109562_kind_phys, 1.093648_kind_phys, 1.056438_kind_phys, 1.007310_kind_phys, 0.951960_kind_phys, &
      0.894453_kind_phys, 0.837364_kind_phys, 0.781742_kind_phys, 0.727415_kind_phys, 0.673614_kind_phys/
    data(bm2ji(4, 4, ibeta), ibeta=1, 10)/ &
      1.423321_kind_phys, 1.417557_kind_phys, 1.402442_kind_phys, 1.379079_kind_phys, 1.347687_kind_phys, &
      1.308075_kind_phys, 1.259703_kind_phys, 1.201983_kind_phys, 1.134778_kind_phys, 1.058878_kind_phys/
    data(bm2ji(4, 5, ibeta), ibeta=1, 10)/ &
      1.933434_kind_phys, 1.944347_kind_phys, 1.968765_kind_phys, 1.997653_kind_phys, 2.023054_kind_phys, &
      2.036554_kind_phys, 2.029949_kind_phys, 1.996982_kind_phys, 1.934982_kind_phys, 1.845473_kind_phys/
    data(bm2ji(4, 6, ibeta), ibeta=1, 10)/ &
      2.547772_kind_phys, 2.577105_kind_phys, 2.645918_kind_phys, 2.735407_kind_phys, 2.830691_kind_phys, &
      2.917268_kind_phys, 2.981724_kind_phys, 3.013684_kind_phys, 3.007302_kind_phys, 2.961560_kind_phys/
    data(bm2ji(4, 7, ibeta), ibeta=1, 10)/ &
      3.101817_kind_phys, 3.139271_kind_phys, 3.225851_kind_phys, 3.336402_kind_phys, 3.453409_kind_phys, &
      3.563116_kind_phys, 3.655406_kind_phys, 3.724014_kind_phys, 3.766113_kind_phys, 3.781394_kind_phys/
    data(bm2ji(4, 8, ibeta), ibeta=1, 10)/ &
      3.540920_kind_phys, 3.573780_kind_phys, 3.647439_kind_phys, 3.737365_kind_phys, 3.828468_kind_phys, &
      3.911436_kind_phys, 3.981317_kind_phys, 4.036345_kind_phys, 4.076749_kind_phys, 4.103751_kind_phys/
    data(bm2ji(4, 9, ibeta), ibeta=1, 10)/ &
      3.856771_kind_phys, 3.879363_kind_phys, 3.928579_kind_phys, 3.986207_kind_phys, 4.042173_kind_phys, &
      4.091411_kind_phys, 4.132041_kind_phys, 4.164052_kind_phys, 4.188343_kind_phys, 4.206118_kind_phys/
    data(bm2ji(4, 10, ibeta), ibeta=1, 10)/ &
      4.053923_kind_phys, 4.067191_kind_phys, 4.095509_kind_phys, 4.127698_kind_phys, 4.158037_kind_phys, &
      4.184055_kind_phys, 4.205135_kind_phys, 4.221592_kind_phys, 4.234115_kind_phys, 4.243463_kind_phys/
    data(bm2ji(5, 1, ibeta), ibeta=1, 10)/ &
      0.935846_kind_phys, 0.906814_kind_phys, 0.843358_kind_phys, 0.768710_kind_phys, 0.695885_kind_phys, &
      0.631742_kind_phys, 0.579166_kind_phys, 0.538471_kind_phys, 0.508410_kind_phys, 0.486863_kind_phys/
    data(bm2ji(5, 2, ibeta), ibeta=1, 10)/ &
      0.988308_kind_phys, 0.959524_kind_phys, 0.896482_kind_phys, 0.821986_kind_phys, 0.748887_kind_phys, &
      0.684168_kind_phys, 0.630908_kind_phys, 0.589516_kind_phys, 0.558676_kind_phys, 0.536056_kind_phys/
    data(bm2ji(5, 3, ibeta), ibeta=1, 10)/ &
      1.133795_kind_phys, 1.107139_kind_phys, 1.048168_kind_phys, 0.977258_kind_phys, 0.906341_kind_phys, &
      0.842477_kind_phys, 0.789093_kind_phys, 0.746731_kind_phys, 0.713822_kind_phys, 0.687495_kind_phys/
    data(bm2ji(5, 4, ibeta), ibeta=1, 10)/ &
      1.405692_kind_phys, 1.385781_kind_phys, 1.340706_kind_phys, 1.284776_kind_phys, 1.227085_kind_phys, &
      1.173532_kind_phys, 1.127008_kind_phys, 1.087509_kind_phys, 1.052712_kind_phys, 1.018960_kind_phys/
    data(bm2ji(5, 5, ibeta), ibeta=1, 10)/ &
      1.884992_kind_phys, 1.879859_kind_phys, 1.868463_kind_phys, 1.854995_kind_phys, 1.841946_kind_phys, &
      1.829867_kind_phys, 1.816972_kind_phys, 1.799319_kind_phys, 1.771754_kind_phys, 1.729406_kind_phys/
    data(bm2ji(5, 6, ibeta), ibeta=1, 10)/ &
      2.592275_kind_phys, 2.612268_kind_phys, 2.661698_kind_phys, 2.731803_kind_phys, 2.815139_kind_phys, &
      2.901659_kind_phys, 2.978389_kind_phys, 3.031259_kind_phys, 3.048045_kind_phys, 3.021122_kind_phys/
    data(bm2ji(5, 7, ibeta), ibeta=1, 10)/ &
      3.390321_kind_phys, 3.435519_kind_phys, 3.545615_kind_phys, 3.698419_kind_phys, 3.876958_kind_phys, &
      4.062790_kind_phys, 4.236125_kind_phys, 4.378488_kind_phys, 4.475619_kind_phys, 4.519170_kind_phys/
    data(bm2ji(5, 8, ibeta), ibeta=1, 10)/ &
      4.161376_kind_phys, 4.216558_kind_phys, 4.346896_kind_phys, 4.519451_kind_phys, 4.711107_kind_phys, &
      4.902416_kind_phys, 5.077701_kind_phys, 5.226048_kind_phys, 5.341423_kind_phys, 5.421764_kind_phys/
    data(bm2ji(5, 9, ibeta), ibeta=1, 10)/ &
      4.843961_kind_phys, 4.892035_kind_phys, 5.001492_kind_phys, 5.138515_kind_phys, 5.281684_kind_phys, &
      5.416805_kind_phys, 5.535493_kind_phys, 5.634050_kind_phys, 5.712063_kind_phys, 5.770996_kind_phys/
    data(bm2ji(5, 10, ibeta), ibeta=1, 10)/ &
      5.352093_kind_phys, 5.385119_kind_phys, 5.458056_kind_phys, 5.545311_kind_phys, 5.632162_kind_phys, &
      5.710566_kind_phys, 5.777005_kind_phys, 5.830863_kind_phys, 5.873123_kind_phys, 5.905442_kind_phys/
    data(bm2ji(6, 1, ibeta), ibeta=1, 10)/ &
      0.964038_kind_phys, 0.930794_kind_phys, 0.859433_kind_phys, 0.777776_kind_phys, 0.700566_kind_phys, &
      0.634671_kind_phys, 0.582396_kind_phys, 0.543656_kind_phys, 0.517284_kind_phys, 0.501694_kind_phys/
    data(bm2ji(6, 2, ibeta), ibeta=1, 10)/ &
      1.013416_kind_phys, 0.979685_kind_phys, 0.907197_kind_phys, 0.824135_kind_phys, 0.745552_kind_phys, &
      0.678616_kind_phys, 0.625870_kind_phys, 0.587348_kind_phys, 0.561864_kind_phys, 0.547674_kind_phys/
    data(bm2ji(6, 3, ibeta), ibeta=1, 10)/ &
      1.145452_kind_phys, 1.111457_kind_phys, 1.038152_kind_phys, 0.953750_kind_phys, 0.873724_kind_phys, &
      0.805955_kind_phys, 0.753621_kind_phys, 0.717052_kind_phys, 0.694920_kind_phys, 0.684910_kind_phys/
    data(bm2ji(6, 4, ibeta), ibeta=1, 10)/ &
      1.376547_kind_phys, 1.345004_kind_phys, 1.276415_kind_phys, 1.196704_kind_phys, 1.121091_kind_phys, &
      1.058249_kind_phys, 1.012197_kind_phys, 0.983522_kind_phys, 0.970323_kind_phys, 0.968933_kind_phys/
    data(bm2ji(6, 5, ibeta), ibeta=1, 10)/ &
      1.778801_kind_phys, 1.755897_kind_phys, 1.706074_kind_phys, 1.649008_kind_phys, 1.597602_kind_phys, &
      1.560087_kind_phys, 1.540365_kind_phys, 1.538205_kind_phys, 1.549738_kind_phys, 1.568333_kind_phys/
    data(bm2ji(6, 6, ibeta), ibeta=1, 10)/ &
      2.447603_kind_phys, 2.445172_kind_phys, 2.443762_kind_phys, 2.451842_kind_phys, 2.475877_kind_phys, &
      2.519039_kind_phys, 2.580118_kind_phys, 2.653004_kind_phys, 2.727234_kind_phys, 2.789738_kind_phys/
    data(bm2ji(6, 7, ibeta), ibeta=1, 10)/ &
      3.368490_kind_phys, 3.399821_kind_phys, 3.481357_kind_phys, 3.606716_kind_phys, 3.772101_kind_phys, &
      3.969416_kind_phys, 4.184167_kind_phys, 4.396163_kind_phys, 4.582502_kind_phys, 4.721838_kind_phys/
    data(bm2ji(6, 8, ibeta), ibeta=1, 10)/ &
      4.426458_kind_phys, 4.489861_kind_phys, 4.648250_kind_phys, 4.877510_kind_phys, 5.160698_kind_phys, &
      5.477495_kind_phys, 5.803123_kind_phys, 6.111250_kind_phys, 6.378153_kind_phys, 6.586050_kind_phys/
    data(bm2ji(6, 9, ibeta), ibeta=1, 10)/ &
      5.568061_kind_phys, 5.644988_kind_phys, 5.829837_kind_phys, 6.081532_kind_phys, 6.371214_kind_phys, &
      6.672902_kind_phys, 6.963737_kind_phys, 7.226172_kind_phys, 7.449199_kind_phys, 7.627886_kind_phys/
    data(bm2ji(6, 10, ibeta), ibeta=1, 10)/ &
      6.639152_kind_phys, 6.707020_kind_phys, 6.863974_kind_phys, 7.065285_kind_phys, 7.281744_kind_phys, &
      7.492437_kind_phys, 7.683587_kind_phys, 7.847917_kind_phys, 7.983296_kind_phys, 8.090977_kind_phys/
    data(bm2ji(7, 1, ibeta), ibeta=1, 10)/ &
      0.980853_kind_phys, 0.945724_kind_phys, 0.871244_kind_phys, 0.787311_kind_phys, 0.708818_kind_phys, &
      0.641987_kind_phys, 0.588462_kind_phys, 0.547823_kind_phys, 0.518976_kind_phys, 0.500801_kind_phys/
    data(bm2ji(7, 2, ibeta), ibeta=1, 10)/ &
      1.026738_kind_phys, 0.990726_kind_phys, 0.914306_kind_phys, 0.828140_kind_phys, 0.747637_kind_phys, &
      0.679351_kind_phys, 0.625127_kind_phys, 0.584662_kind_phys, 0.556910_kind_phys, 0.540749_kind_phys/
    data(bm2ji(7, 3, ibeta), ibeta=1, 10)/ &
      1.146496_kind_phys, 1.108808_kind_phys, 1.028695_kind_phys, 0.938291_kind_phys, 0.854101_kind_phys, &
      0.783521_kind_phys, 0.728985_kind_phys, 0.690539_kind_phys, 0.667272_kind_phys, 0.657977_kind_phys/
    data(bm2ji(7, 4, ibeta), ibeta=1, 10)/ &
      1.344846_kind_phys, 1.306434_kind_phys, 1.224543_kind_phys, 1.132031_kind_phys, 1.046571_kind_phys, &
      0.976882_kind_phys, 0.926488_kind_phys, 0.896067_kind_phys, 0.884808_kind_phys, 0.891027_kind_phys/
    data(bm2ji(7, 5, ibeta), ibeta=1, 10)/ &
      1.670227_kind_phys, 1.634583_kind_phys, 1.558421_kind_phys, 1.472939_kind_phys, 1.396496_kind_phys, &
      1.339523_kind_phys, 1.307151_kind_phys, 1.300882_kind_phys, 1.319622_kind_phys, 1.360166_kind_phys/
    data(bm2ji(7, 6, ibeta), ibeta=1, 10)/ &
      2.224548_kind_phys, 2.199698_kind_phys, 2.148284_kind_phys, 2.095736_kind_phys, 2.059319_kind_phys, &
      2.050496_kind_phys, 2.075654_kind_phys, 2.136382_kind_phys, 2.229641_kind_phys, 2.347958_kind_phys/
    data(bm2ji(7, 7, ibeta), ibeta=1, 10)/ &
      3.104483_kind_phys, 3.105947_kind_phys, 3.118398_kind_phys, 3.155809_kind_phys, 3.230427_kind_phys, &
      3.350585_kind_phys, 3.519071_kind_phys, 3.731744_kind_phys, 3.976847_kind_phys, 4.235616_kind_phys/
    data(bm2ji(7, 8, ibeta), ibeta=1, 10)/ &
      4.288426_kind_phys, 4.331456_kind_phys, 4.447024_kind_phys, 4.633023_kind_phys, 4.891991_kind_phys, &
      5.221458_kind_phys, 5.610060_kind_phys, 6.036467_kind_phys, 6.471113_kind_phys, 6.880462_kind_phys/
    data(bm2ji(7, 9, ibeta), ibeta=1, 10)/ &
      5.753934_kind_phys, 5.837061_kind_phys, 6.048530_kind_phys, 6.363800_kind_phys, 6.768061_kind_phys, &
      7.241280_kind_phys, 7.755346_kind_phys, 8.276666_kind_phys, 8.771411_kind_phys, 9.210826_kind_phys/
    data(bm2ji(7, 10, ibeta), ibeta=1, 10)/ &
      7.466219_kind_phys, 7.568810_kind_phys, 7.819032_kind_phys, 8.168340_kind_phys, 8.582973_kind_phys, &
      9.030174_kind_phys, 9.478159_kind_phys, 9.899834_kind_phys, 10.275940_kind_phys, 10.595910_kind_phys/
    data(bm2ji(8, 1, ibeta), ibeta=1, 10)/ &
      0.990036_kind_phys, 0.954782_kind_phys, 0.880531_kind_phys, 0.797334_kind_phys, 0.719410_kind_phys, &
      0.652220_kind_phys, 0.596923_kind_phys, 0.552910_kind_phys, 0.519101_kind_phys, 0.494529_kind_phys/
    data(bm2ji(8, 2, ibeta), ibeta=1, 10)/ &
      1.032428_kind_phys, 0.996125_kind_phys, 0.919613_kind_phys, 0.833853_kind_phys, 0.753611_kind_phys, &
      0.684644_kind_phys, 0.628260_kind_phys, 0.583924_kind_phys, 0.550611_kind_phys, 0.527407_kind_phys/
    data(bm2ji(8, 3, ibeta), ibeta=1, 10)/ &
      1.141145_kind_phys, 1.102521_kind_phys, 1.021017_kind_phys, 0.929667_kind_phys, 0.844515_kind_phys, &
      0.772075_kind_phys, 0.714086_kind_phys, 0.670280_kind_phys, 0.639824_kind_phys, 0.621970_kind_phys/
    data(bm2ji(8, 4, ibeta), ibeta=1, 10)/ &
      1.314164_kind_phys, 1.273087_kind_phys, 1.186318_kind_phys, 1.089208_kind_phys, 0.999476_kind_phys, &
      0.924856_kind_phys, 0.867948_kind_phys, 0.829085_kind_phys, 0.807854_kind_phys, 0.803759_kind_phys/
    data(bm2ji(8, 5, ibeta), ibeta=1, 10)/ &
      1.580611_kind_phys, 1.538518_kind_phys, 1.449529_kind_phys, 1.350459_kind_phys, 1.260910_kind_phys, &
      1.190526_kind_phys, 1.143502_kind_phys, 1.121328_kind_phys, 1.124274_kind_phys, 1.151974_kind_phys/
    data(bm2ji(8, 6, ibeta), ibeta=1, 10)/ &
      2.016773_kind_phys, 1.977721_kind_phys, 1.895727_kind_phys, 1.806974_kind_phys, 1.732891_kind_phys, &
      1.685937_kind_phys, 1.673026_kind_phys, 1.697656_kind_phys, 1.761039_kind_phys, 1.862391_kind_phys/
    data(bm2ji(8, 7, ibeta), ibeta=1, 10)/ &
      2.750093_kind_phys, 2.723940_kind_phys, 2.672854_kind_phys, 2.628264_kind_phys, 2.612250_kind_phys, &
      2.640406_kind_phys, 2.723211_kind_phys, 2.866599_kind_phys, 3.071893_kind_phys, 3.335217_kind_phys/
    data(bm2ji(8, 8, ibeta), ibeta=1, 10)/ &
      3.881905_kind_phys, 3.887143_kind_phys, 3.913667_kind_phys, 3.981912_kind_phys, 4.111099_kind_phys, &
      4.316575_kind_phys, 4.608146_kind_phys, 4.988157_kind_phys, 5.449592_kind_phys, 5.974848_kind_phys/
    data(bm2ji(8, 9, ibeta), ibeta=1, 10)/ &
      5.438870_kind_phys, 5.492742_kind_phys, 5.640910_kind_phys, 5.886999_kind_phys, 6.241641_kind_phys, &
      6.710609_kind_phys, 7.289480_kind_phys, 7.960725_kind_phys, 8.693495_kind_phys, 9.446644_kind_phys/
    data(bm2ji(8, 10, ibeta), ibeta=1, 10)/ &
      7.521152_kind_phys, 7.624621_kind_phys, 7.892039_kind_phys, 8.300444_kind_phys, 8.839787_kind_phys, &
      9.493227_kind_phys, 10.231770_kind_phys, 11.015642_kind_phys, 11.799990_kind_phys, 12.542260_kind_phys/
    data(bm2ji(9, 1, ibeta), ibeta=1, 10)/ &
      0.994285_kind_phys, 0.960012_kind_phys, 0.887939_kind_phys, 0.807040_kind_phys, 0.730578_kind_phys, &
      0.663410_kind_phys, 0.606466_kind_phys, 0.559137_kind_phys, 0.520426_kind_phys, 0.489429_kind_phys/
    data(bm2ji(9, 2, ibeta), ibeta=1, 10)/ &
      1.033505_kind_phys, 0.998153_kind_phys, 0.923772_kind_phys, 0.840261_kind_phys, 0.761383_kind_phys, &
      0.692242_kind_phys, 0.633873_kind_phys, 0.585709_kind_phys, 0.546777_kind_phys, 0.516215_kind_phys/
    data(bm2ji(9, 3, ibeta), ibeta=1, 10)/ &
      1.132774_kind_phys, 1.094907_kind_phys, 1.015161_kind_phys, 0.925627_kind_phys, 0.841293_kind_phys, &
      0.767888_kind_phys, 0.706741_kind_phys, 0.657439_kind_phys, 0.619135_kind_phys, 0.591119_kind_phys/
    data(bm2ji(9, 4, ibeta), ibeta=1, 10)/ &
      1.286308_kind_phys, 1.245273_kind_phys, 1.158809_kind_phys, 1.061889_kind_phys, 0.971208_kind_phys, &
      0.893476_kind_phys, 0.830599_kind_phys, 0.782561_kind_phys, 0.748870_kind_phys, 0.729198_kind_phys/
    data(bm2ji(9, 5, ibeta), ibeta=1, 10)/ &
      1.511105_kind_phys, 1.467141_kind_phys, 1.374520_kind_phys, 1.271162_kind_phys, 1.175871_kind_phys, &
      1.096887_kind_phys, 1.037243_kind_phys, 0.997820_kind_phys, 0.978924_kind_phys, 0.980962_kind_phys/
    data(bm2ji(9, 6, ibeta), ibeta=1, 10)/ &
      1.857468_kind_phys, 1.812177_kind_phys, 1.717002_kind_phys, 1.612197_kind_phys, 1.519171_kind_phys, &
      1.448660_kind_phys, 1.405871_kind_phys, 1.393541_kind_phys, 1.413549_kind_phys, 1.467532_kind_phys/
    data(bm2ji(9, 7, ibeta), ibeta=1, 10)/ &
      2.430619_kind_phys, 2.388452_kind_phys, 2.301326_kind_phys, 2.210241_kind_phys, 2.139724_kind_phys, &
      2.104571_kind_phys, 2.114085_kind_phys, 2.174696_kind_phys, 2.291294_kind_phys, 2.467500_kind_phys/
    data(bm2ji(9, 8, ibeta), ibeta=1, 10)/ &
      3.385332_kind_phys, 3.357690_kind_phys, 3.306611_kind_phys, 3.269804_kind_phys, 3.274462_kind_phys, &
      3.340862_kind_phys, 3.484609_kind_phys, 3.717740_kind_phys, 4.048748_kind_phys, 4.481588_kind_phys/
    data(bm2ji(9, 9, ibeta), ibeta=1, 10)/ &
      4.850497_kind_phys, 4.858280_kind_phys, 4.896008_kind_phys, 4.991467_kind_phys, 5.171511_kind_phys, &
      5.459421_kind_phys, 5.873700_kind_phys, 6.426128_kind_phys, 7.119061_kind_phys, 7.942603_kind_phys/
    data(bm2ji(9, 10, ibeta), ibeta=1, 10)/ &
      6.957098_kind_phys, 7.020164_kind_phys, 7.197272_kind_phys, 7.499331_kind_phys, 7.946554_kind_phys, &
      8.555048_kind_phys, 9.330503_kind_phys, 10.263610_kind_phys, 11.327454_kind_phys, 12.478332_kind_phys/
    data(bm2ji(10, 1, ibeta), ibeta=1, 10)/ &
      0.994567_kind_phys, 0.961842_kind_phys, 0.892854_kind_phys, 0.814874_kind_phys, 0.740198_kind_phys, &
      0.673303_kind_phys, 0.615105_kind_phys, 0.565139_kind_phys, 0.522558_kind_phys, 0.486556_kind_phys/
    data(bm2ji(10, 2, ibeta), ibeta=1, 10)/ &
      1.031058_kind_phys, 0.997292_kind_phys, 0.926082_kind_phys, 0.845571_kind_phys, 0.768501_kind_phys, &
      0.699549_kind_phys, 0.639710_kind_phys, 0.588538_kind_phys, 0.545197_kind_phys, 0.508894_kind_phys/
    data(bm2ji(10, 3, ibeta), ibeta=1, 10)/ &
      1.122535_kind_phys, 1.086287_kind_phys, 1.009790_kind_phys, 0.923292_kind_phys, 0.840626_kind_phys, &
      0.766982_kind_phys, 0.703562_kind_phys, 0.650004_kind_phys, 0.605525_kind_phys, 0.569411_kind_phys/
    data(bm2ji(10, 4, ibeta), ibeta=1, 10)/ &
      1.261142_kind_phys, 1.221555_kind_phys, 1.137979_kind_phys, 1.043576_kind_phys, 0.953745_kind_phys, &
      0.874456_kind_phys, 0.807292_kind_phys, 0.752109_kind_phys, 0.708326_kind_phys, 0.675477_kind_phys/
    data(bm2ji(10, 5, ibeta), ibeta=1, 10)/ &
      1.456711_kind_phys, 1.413432_kind_phys, 1.322096_kind_phys, 1.219264_kind_phys, 1.122319_kind_phys, &
      1.038381_kind_phys, 0.969743_kind_phys, 0.916811_kind_phys, 0.879544_kind_phys, 0.858099_kind_phys/
    data(bm2ji(10, 6, ibeta), ibeta=1, 10)/ &
      1.741792_kind_phys, 1.695157_kind_phys, 1.596897_kind_phys, 1.487124_kind_phys, 1.385734_kind_phys, &
      1.301670_kind_phys, 1.238638_kind_phys, 1.198284_kind_phys, 1.181809_kind_phys, 1.190689_kind_phys/
    data(bm2ji(10, 7, ibeta), ibeta=1, 10)/ &
      2.190197_kind_phys, 2.141721_kind_phys, 2.040226_kind_phys, 1.929245_kind_phys, 1.832051_kind_phys, &
      1.760702_kind_phys, 1.721723_kind_phys, 1.719436_kind_phys, 1.757705_kind_phys, 1.840677_kind_phys/
    data(bm2ji(10, 8, ibeta), ibeta=1, 10)/ &
      2.940764_kind_phys, 2.895085_kind_phys, 2.801873_kind_phys, 2.707112_kind_phys, 2.638603_kind_phys, &
      2.613764_kind_phys, 2.644686_kind_phys, 2.741255_kind_phys, 2.912790_kind_phys, 3.168519_kind_phys/
    data(bm2ji(10, 9, ibeta), ibeta=1, 10)/ &
      4.186191_kind_phys, 4.155844_kind_phys, 4.101953_kind_phys, 4.069102_kind_phys, 4.089886_kind_phys, &
      4.189530_kind_phys, 4.389145_kind_phys, 4.707528_kind_phys, 5.161567_kind_phys, 5.765283_kind_phys/
    data(bm2ji(10, 10, ibeta), ibeta=1, 10)/ &
      6.119526_kind_phys, 6.127611_kind_phys, 6.171174_kind_phys, 6.286528_kind_phys, 6.508738_kind_phys, &
      6.869521_kind_phys, 7.396912_kind_phys, 8.113749_kind_phys, 9.034683_kind_phys, 10.162190_kind_phys/

! *** end of data statements.

! *** start calculations:

    constii = abs(half*(two)**two3rds - one)
    sqrttwo = sqrt(two)
    dlgsqt2 = one/log(sqrttwo)

    esat01 = exp(0.125_kind_phys*xxlsgat*xxlsgat)
    esac01 = exp(0.125_kind_phys*xxlsgac*xxlsgac)

    esat04 = esat01**4
    esac04 = esac01**4

    esat05 = esat04*esat01
    esac05 = esac04*esac01

    esat08 = esat04*esat04
    esac08 = esac04*esac04

    esat09 = esat08*esat01
    esac09 = esac08*esac01

    esat16 = esat08*esat08
    esac16 = esac08*esac08

    esat20 = esat16*esat04
    esac20 = esac16*esac04

    esat24 = esat20*esat04
    esac24 = esac20*esac04

    esat25 = esat20*esat05
    esac25 = esac20*esac05

    esat36 = esat20*esat16
    esac36 = esac20*esac16

    esat49 = esat24*esat25

    esat64 = esat20*esat20*esat24
    esac64 = esac20*esac20*esac24

    esat100 = esat64*esat36

    dgat2 = dgatk*dgatk
    dgat3 = dgatk*dgatk*dgatk
    dgac2 = dgacc*dgacc
    dgac3 = dgacc*dgacc*dgacc

    sqdgat = sqrt(dgatk)
    sqdgac = sqrt(dgacc)
    sqdgat5 = dgat2*sqdgat
    sqdgac5 = dgac2*sqdgac
    sqdgat7 = dgat3*sqdgat

    xm2at = dgat2*esat16
    xm3at = dgat3*esat36

    xm2ac = dgac2*esac16
    xm3ac = dgac3*esac36

! *** for the free molecular regime:  page h.3 of whitby et al. (1991)

    r = sqdgac/sqdgat
    r2 = r*r
    r3 = r2*r
    rx4 = r2*r2
    r5 = r3*r2
    r6 = r3*r3
    rx8 = rx4*rx4
    ri1 = one/r
    ri2 = one/r2
    ri3 = one/r3
    ri4 = ri2*ri2
    kngat = two*lamda/dgatk
    kngac = two*lamda/dgacc

! *** calculate ratio of geometric mean diameters
    rat = dgacc/dgatk
! *** trap subscripts for bm0 and bm0i, between 1 and 10
!     see page h.5 of whitby et al. (1991)

    n2n = max(1, min(10, &
                     nint(4.0_kind_phys*(sgatk - 0.75_kind_phys))))

    n2a = max(1, min(10, &
                     nint(4.0_kind_phys*(sgacc - 0.75_kind_phys))))

    n1 = max(1, min(10, &
                    1 + nint(dlgsqt2*log(rat))))

! *** intermodal coagulation

! *** set up for zeroeth moment

! *** near-continuum form:  equation h.10a of whitby et al. (1991)

    coagnc0 = knc*( &
              two + a*(kngat*(esat04 + r2*esat16*esac04) &
                       + kngac*(esac04 + ri2*esac16*esat04)) &
              + (r2 + ri2)*esat04*esac04)

! *** free-molecular form:  equation h.7a of whitby et al. (1991)

    coagfm0 = kfmatac*sqdgat*bm0ij(n1, n2n, n2a)*( &
              esat01 + r*esac01 + two*r2*esat01*esac04 &
              + rx4*esat09*esac16 + ri3*esat16*esac09 &
              + two*ri1*esat04 + esac01)

! *** loss to accumulation mode

! *** harmonic mean

    coagatac0 = coagnc0*coagfm0/(coagnc0 + coagfm0)

    qn12 = coagatac0

! *** set up for second moment
!      the second moment equations are new and begin with equations a1
!     through a4 of binkowski and shankar (1995). after some algebraic
!     rearrangement and application of the extended mean value theorem
!     of integral calculus, equations are obtained that can be solved
!     analytically with correction factors as has been done by
!     whitby et al. (1991)

! *** the term ( dp1 + dp2 ) ** (2/3) in equations a3 and a4 of
!     binkowski and shankar (1995) is approximated by
!     (dgat ** 3 + dgac **3 ) ** 2/3

! *** near-continuum form

    i1nc = knc*dgat2*( &
           two*esat16 &
           + r2*esat04*esac04 &
           + ri2*esat36*esac04 &
           + a*kngat*( &
           esat04 &
           + ri2*esat16*esac04 &
           + ri4*esat36*esac16 &
           + r2*esac04))

! *** free-molecular form

    i1fm = kfmatac*sqdgat5*bm2ij(n1, n2n, n2a)*( &
           esat25 &
           + two*r2*esat09*esac04 &
           + rx4*esat01*esac16 &
           + ri3*esat64*esac09 &
           + two*ri1*esat36*esac01 &
           + r*esat16*esac01)

! *** loss to accumulation mode

! *** harmonic mean

    i1 = (i1fm*i1nc)/(i1fm + i1nc)

    coagatac2 = i1

    qs12 = coagatac2

! *** gain by accumulation mode

    coagacat2 = ((one + r6)**two3rds - rx4)*i1

    qs21 = coagacat2*bm2ji(n1, n2n, n2a)

! *** set up for third moment

! *** near-continuum form: equation h.10b of whitby et al. (1991)

    coagnc3 = knc*dgat3*( &
              two*esat36 &
              + a*kngat*(esat16 + r2*esat04*esac04) &
              + a*kngac*(esat36*esac04 + ri2*esat64*esac16) &
              + r2*esat16*esac04 + ri2*esat64*esac04)

! *** free_molecular form: equation h.7b of whitby et al. (1991)

    coagfm3 = kfmatac*sqdgat7*bm3i(n1, n2n, n2a)*( &
              esat49 &
              + r*esat36*esac01 &
              + two*r2*esat25*esac04 &
              + rx4*esat09*esac16 &
              + ri3*esat100*esac09 &
              + two*ri1*esat64*esac01)

! *** gain by accumulation mode = loss from aitken mode

! *** harmonic mean

    coagatac3 = coagnc3*coagfm3/(coagnc3 + coagfm3)

    qv12 = coagatac3

! *** intramodal coagulation

! *** zeroeth moment

! *** aitken mode

! *** near-continuum form: equation h.12a of whitby et al. (1991)

    coagnc_at = knc*(one + esat08 + a*kngat*(esat20 + esat04))

! *** free-molecular form: equation h.11a of whitby et al. (1991)

    coagfm_at = kfmat*sqdgat*bm0(n2n)* &
                (esat01 + esat25 + two*esat05)

! *** harmonic mean

    coagatat0 = coagfm_at*coagnc_at/(coagfm_at + coagnc_at)

    qn11 = coagatat0

! *** accumulation mode

! *** near-continuum form: equation h.12a of whitby et al. (1991)

    coagnc_ac = knc*(one + esac08 + a*kngac*(esac20 + esac04))

! *** free-molecular form: equation h.11a of whitby et al. (1991)

    coagfm_ac = kfmac*sqdgac*bm0(n2a)* &
                (esac01 + esac25 + two*esac05)

! *** harmonic mean

    coagacac0 = coagfm_ac*coagnc_ac/(coagfm_ac + coagnc_ac)

    qn22 = coagacac0

! *** set up for second moment
!      the second moment equations are new and begin with 3.11a on page
!     45 of whitby et al. (1991). after some algebraic rearrangement and
!     application of the extended mean value theorem of integral calculus
!     equations are obtained that can be solved analytically with
!     correction factors as has been done by whitby et al. (1991)

    ! *** aitken mode
    ! *** near-continuum
    i1nc_at = knc*dgat2*( &
              two*esat16 &
              + esat04*esat04 &
              + esat36*esat04 &
              + a*kngat*( &
              two*esat04 &
              + esat16*esat04 &
              + esat36*esat16))

    ! *** free- molecular form
    i1fm_at = kfmat*sqdgat5*bm2ii(n2n)*( &
              esat25 &
              + two*esat09*esat04 &
              + esat01*esat16 &
              + esat64*esat09 &
              + two*esat36*esat01 &
              + esat16*esat01)

    i1_at = (i1nc_at*i1fm_at)/(i1nc_at + i1fm_at)

    coagatat2 = constii*i1_at

    qs11 = coagatat2*bm2iitt(n2n)

    ! *** accumulation mode
    ! *** near-continuum
    i1nc_ac = knc*dgac2*( &
              two*esac16 &
              + esac04*esac04 &
              + esac36*esac04 &
              + a*kngac*( &
              two*esac04 &
              + esac16*esac04 &
              + esac36*esac16))

    ! *** free- molecular form
    i1fm_ac = kfmac*sqdgac5*bm2ii(n2a)*( &
              esac25 &
              + two*esac09*esac04 &
              + esac01*esac16 &
              + esac64*esac09 &
              + two*esac36*esac01 &
              + esac16*esac01)

    i1_ac = (i1nc_ac*i1fm_ac)/(i1nc_ac + i1fm_ac)

    coagacac2 = constii*i1_ac

    qs22 = coagacac2*bm2iitt(n2a)

  end subroutine getcoags

end module modal_aero_coag

