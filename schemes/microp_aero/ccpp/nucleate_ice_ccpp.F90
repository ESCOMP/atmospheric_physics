! Unified CCPP scheme for ice nucleation via the nucleati kernel.
! Uses the abstract interface to support BAM, MAM, and CARMA.
!
! This scheme replaces nucleate_ice_cam_calc for CAM-SIMA.
module nucleate_ice_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: nucleate_ice_ccpp_init
  public :: nucleate_ice_ccpp_run

  ! Module state from init
  logical :: use_preexisting_ice_ = .false.
  logical :: nucleate_ice_use_troplev_ = .false.
  real(kind_phys) :: nucleate_ice_subgrid_
  real(kind_phys) :: nucleate_ice_subgrid_strat_
  real(kind_phys) :: nucleate_ice_strat_

  ! Hardcoded scheme constants (previously CAM-side parameters/flags).
  ! mincld_: floor for cloud fraction used in ice nucleation. CAM value was
  !          real(r8), parameter :: mincld = 0.0001_r8 in microp_aero.F90.
  ! use_nucleati_tendencies_: PUMAS v1.21+ (CAM7+) tendency path. Only supported
  !          path in the CCPP version; legacy path will abort as not implemented.
  real(kind_phys), parameter :: mincld_ = 0.0001_kind_phys
  logical,         parameter :: use_nucleati_tendencies_ = .true.

  ! Aerosol model flags (determined at init from aerosol_instances)
  logical :: clim_modal_aero_ = .false.
  logical :: clim_carma_aero_ = .false.

  ! index into aerosol_instances
  ! NOTE: this replicates the microp_aero flow where only one aerosol
  ! model/state is selected for use by microp_aero (ice nucleation, activation)
  ! even if multiple aerosol models are active at the same time.
  ! If modal or CARMA are active, they are used, otherwise bulk. (hplin, 4/7/26)
  !
  ! This module also can run without an aerosol model, in aquaplanet compsets
  ! with two-moment microphysics; nucleatei runs Meyers nucleation deposition
  ! even without aerosol, in which case iaermod_selected_ is -1.
  integer :: iaermod_selected_ = -1

  ! Modal/CARMA-specific: whether cloud-borne tendencies are needed
  logical :: clim_modal_carma_ = .false.
  logical :: prog_modal_aero_  = .false.

  ! Aerosol bins/species dimensions (set in init, used in run for transfer)
  integer :: nbins_   = 0
  integer :: nmaxspc_ = 0

  ! Constituent tendency indices for interstitial-to-cloud-borne aerosol transfer.
  ! Allocated in init if clim_modal_carma_ .and. use_preexisting_ice_.
  ! aer_cnst_idx_(ibin, 0)    = index for bin number/MMR (species 0)
  ! aer_cnst_idx_(ibin, ispc) = index for species ispc MMR
  ! Index is into CCPP constituent array; -1 means not transported.
  integer, allocatable :: aer_cnst_idx_(:,:)

contains

!> \section arg_table_nucleate_ice_ccpp_init Argument Table
!! \htmlinclude nucleate_ice_ccpp_init.html
  subroutine nucleate_ice_ccpp_init( &
    iulog, &
    const_props, &
    use_preexisting_ice, use_hetfrz_classnuc, &
    nucleate_ice_incloud, &
    pi, &
    nucleate_ice_subgrid, nucleate_ice_subgrid_strat, &
    nucleate_ice_strat, nucleate_ice_use_troplev, &
    prog_modal_aero, &
    errmsg, errflg)

    use nucleate_ice,           only: nucleati_init
    use aerosol_instances_mod,  only: aerosol_instances_is_active, &
                                      aerosol_instances_get_props, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,          only: ccpp_const_get_idx

    ! Input arguments
    integer,          intent(in)  :: iulog
    type(ccpp_constituent_prop_ptr_t), &
                      intent(in)  :: const_props(:)
    logical,          intent(in)  :: use_preexisting_ice
    logical,          intent(in)  :: use_hetfrz_classnuc
    logical,          intent(in)  :: nucleate_ice_incloud
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: nucleate_ice_subgrid
    real(kind_phys),  intent(in)  :: nucleate_ice_subgrid_strat
    real(kind_phys),  intent(in)  :: nucleate_ice_strat
    logical,          intent(in)  :: nucleate_ice_use_troplev
    logical,          intent(in)  :: prog_modal_aero

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: iaermod
    integer :: ibin, ispc, idxtmp
    class(aerosol_properties), pointer :: aprops
    character(len=32) :: tmpname

    errmsg = ''
    errflg = 0

    ! Store namelist parameters as module state
    use_preexisting_ice_        = use_preexisting_ice
    nucleate_ice_use_troplev_   = nucleate_ice_use_troplev
    nucleate_ice_subgrid_       = nucleate_ice_subgrid
    nucleate_ice_subgrid_strat_ = nucleate_ice_subgrid_strat
    nucleate_ice_strat_         = nucleate_ice_strat
    prog_modal_aero_            = prog_modal_aero

    ! Query aerosol model availability
    clim_modal_aero_ = aerosol_instances_is_active('modal')
    clim_carma_aero_ = aerosol_instances_is_active('carma')

    ! Select aerosol model (modal > CARMA > BAM > none)
    iaermod_selected_ = -1
    do iaermod = 1, aerosol_instances_get_num_models()
      aprops => aerosol_instances_get_props(iaermod, list_idx=0)
      if (.not. associated(aprops)) cycle

      if (aprops%model_is('modal') .or. aprops%model_is('CARMA')) then
        iaermod_selected_ = iaermod
        clim_modal_carma_ = .true.
        exit
      else if (aprops%model_is('BAM')) then
        iaermod_selected_ = iaermod
      end if
    end do

    ! Set up aerosol dimensions and constituent indices for transfer
    nbins_ = 0
    nmaxspc_ = 0
    if (iaermod_selected_ > 0) then
      aprops => aerosol_instances_get_props(iaermod_selected_, list_idx=0)
      if (associated(aprops)) then
        nbins_ = aprops%nbins()
        nmaxspc_ = maxval(aprops%nspecies())
      end if
    end if

    ! Set up constituent tendency indices for interstitial-to-cloud-borne
    ! aerosol transfer (MAM/CARMA + use_preexisting_ice only).
    ! Replaces nucleate_ice_cam_init L178-222.
    if (clim_modal_carma_ .and. use_preexisting_ice_ .and. associated(aprops)) then

      allocate(aer_cnst_idx_(nbins_, 0:nmaxspc_), stat=errflg, errmsg=errmsg)
      if(errflg /= 0) return
      aer_cnst_idx_ = -1

      do ibin = 1, nbins_
        if (aprops%icenuc_updates_num(ibin)) then

          ! Get constituent name for this bin (number or bin MMR)
          if (aprops%icenuc_updates_mmr(ibin, 0)) then
            call aprops%amb_mmr_name(ibin, 0, tmpname)
          else
            call aprops%amb_num_name(ibin, tmpname)
          end if

          ! Look up CCPP constituent index by name (standard_name = constituent name)
          call ccpp_const_get_idx(const_props, trim(tmpname), idxtmp, errmsg, errflg)
          if (errflg /= 0) return
          aer_cnst_idx_(ibin, 0) = idxtmp

          ! Iterate over species within the bin
          do ispc = 1, aprops%nspecies(ibin)
            if (aprops%icenuc_updates_mmr(ibin, ispc)) then
              call aprops%amb_mmr_name(ibin, ispc, tmpname)
              call ccpp_const_get_idx(const_props, trim(tmpname), idxtmp, errmsg, errflg)
              if (errflg /= 0) return
              aer_cnst_idx_(ibin, ispc) = idxtmp
            end if
          end do

        end if
      end do

    end if

    ! Initialize the nucleati kernel
    call nucleati_init(use_preexisting_ice, use_hetfrz_classnuc, &
                       nucleate_ice_incloud, iulog, pi, mincld_)

  end subroutine nucleate_ice_ccpp_init

!> \section arg_table_nucleate_ice_ccpp_run Argument Table
!! \htmlinclude nucleate_ice_ccpp_run.html
  subroutine nucleate_ice_ccpp_run( &
    ncol, pver, pcnst, top_lev, dtime,   &
    t, pmid, qv, qc, qi, ni,             &
    rair, tmelt, pi,                     &
    wsubi,                               &
    aist,                             &
    tropLev_chem,                        &
    qsatfac,                             &
    lat,                                 &
    naai, naai_hom,                      &
    nihf, niimm, nidep, nimey,           &
    regm, subgrid_diag, trop_pd,         &
    fhom, wice, weff,                    &
    INnso4, INnbc, INndust, INondust,    &
    INhet, INhom, INFrehom, INFreIN,     &
    ptend_q,                             &
    errmsg, errflg)

    ! portable core science code:
    use nucleate_ice,           only: nucleati

    ! to_be_ccppized dependency:
    use wv_saturation,          only: qsat_water

    ! abstract aerosol interface:
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_state, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod,      only: aerosol_state

    ! Input arguments
    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: pcnst
    integer,          intent(in)  :: top_lev        ! top vertical level for cloud physics [index]
    real(kind_phys),  intent(in)  :: dtime          ! timestep [s]
    real(kind_phys),  intent(in)  :: t(:, :)        ! temperature [K]
    real(kind_phys),  intent(in)  :: pmid(:, :)     ! pressure at layer midpoints [Pa]
    real(kind_phys),  intent(in)  :: qv(:, :)       ! water vapor mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)  :: qc(:, :)       ! cloud liquid mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)  :: qi(:, :)       ! cloud ice mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)  :: ni(:, :)       ! cloud ice number concentration [1 kg-1]

    real(kind_phys),  intent(in)  :: rair           ! gas constant for dry air [J kg-1 K-1]
    real(kind_phys),  intent(in)  :: tmelt          ! freezing point of water [K]
    real(kind_phys),  intent(in)  :: pi             ! pi

    real(kind_phys),  intent(in)  :: wsubi(:, :)    ! subgrid vertical velocity for ice nucleation [m s-1]
    real(kind_phys),  intent(in)  :: aist(:, :)     ! ice cloud fraction [fraction]
    integer,          intent(in)  :: tropLev_chem(:)! chemical tropopause level [index]
    real(kind_phys),  intent(in)  :: qsatfac(:, :)  ! subgrid saturation scaling factor [1]
    real(kind_phys),  intent(in)  :: lat(:)         ! latitude [radians]

    ! Note: below units for naai...nimey are tendencies (with s-1)
    ! as is convention for CAM7 PUMASv1.21+
    ! The legacy path is not supported in CCPP version of nucleate_ice; if it were
    ! to be supported the units will be concentrations (no s-1); they would not be
    ! divided by /dtime in this run phase.

    ! Output arguments
    real(kind_phys),  intent(out) :: naai(:, :)      ! nucleated ice number [kg-1 s-1]
    real(kind_phys),  intent(out) :: naai_hom(:, :)  ! nucleated ice number, hom. freezing only [kg-1 s-1]

    ! Diagnostic outputs
    real(kind_phys),  intent(out) :: nihf(:, :)     ! ice nuclei from hom. freezing [m-3 s-1]
    real(kind_phys),  intent(out) :: niimm(:, :)    ! ice nuclei from immersion freezing [m-3 s-1]
    real(kind_phys),  intent(out) :: nidep(:, :)    ! ice nuclei from deposition nucleation [m-3 s-1]
    real(kind_phys),  intent(out) :: nimey(:, :)    ! ice nuclei from Meyers deposition [m-3 s-1]
    real(kind_phys),  intent(out) :: regm(:, :)     ! nucleation regime temperature threshold [C]
    real(kind_phys),  intent(out) :: subgrid_diag(:, :) ! ice nucleation subgrid saturation factor [1]
    real(kind_phys),  intent(out) :: trop_pd(:, :)  ! chemical tropopause probability [1]

    ! Pre-existing ice diagnostics (output regardless, zeroed if not use_preexisting_ice_)
    real(kind_phys),  intent(out) :: fhom(:, :)     ! fraction of cirrus with hom. freezing [fraction]
    real(kind_phys),  intent(out) :: wice(:, :)     ! vertical velocity reduction from preexisting ice [m/s]
    real(kind_phys),  intent(out) :: weff(:, :)     ! effective vertical velocity for ice nucleation [m/s]

    ! Pre-existing ice aerosol diagnostics (zeroed if not use_preexisting_ice_)
    real(kind_phys),  intent(out) :: INnso4(:, :)   ! so4 number conc. tendency to ice nucleation [m-3 s-1]
    real(kind_phys),  intent(out) :: INnbc(:, :)    ! bc number conc. tendency to ice nucleation [m-3 s-1]
    real(kind_phys),  intent(out) :: INndust(:, :)  ! dust number conc. tendency to ice nucleation [m-3 s-1]
    real(kind_phys),  intent(out) :: INondust(:, :) ! dust number conc. tendency from ice nucleation [m-3 s-1]
    real(kind_phys),  intent(out) :: INhet(:, :)    ! heterogeneous IN number tendency [m-3 s-1]
    real(kind_phys),  intent(out) :: INhom(:, :)    ! homogeneous IN number tendency [m-3 s-1]
    real(kind_phys),  intent(out) :: INFrehom(:, :) ! homogeneous ice nucleation frequency [1]
    real(kind_phys),  intent(out) :: INFreIN(:, :)  ! ice nucleation frequency [1]

    ! Constituent tendencies: CCPP framework applies these via tendency updater.
    ! For MAM/CARMA + use_preexisting_ice: interstitial-to-cloud-borne transfer.
    ! For BAM and aquaplanet: zeroed (no constituent tendencies).
    real(kind_phys),  intent(out) :: ptend_q(:, :, :)  ! constituent tendencies [kg kg-1 s-1] (ncol, pver, pcnst)

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    class(aerosol_properties), pointer :: aeroprops
    class(aerosol_state),      pointer :: aerostate

    integer :: i, k, m, l
    integer :: idxtmp

    real(kind_phys) :: rho(ncol,pver)
    real(kind_phys) :: qs(ncol)
    real(kind_phys) :: es(ncol)
    real(kind_phys) :: gammas(ncol)
    real(kind_phys) :: relhum(ncol,pver)
    real(kind_phys) :: icldm(ncol,pver)
    real(kind_phys) :: subgrid(ncol,pver)

    real(kind_phys) :: dust_num_col(ncol,pver)
    real(kind_phys) :: sulf_num_col(ncol,pver)
    real(kind_phys) :: soot_num_col(ncol,pver)
    real(kind_phys) :: sulf_num_tot_col(ncol,pver)

    real(kind_phys) :: so4_num, dst_num, soot_num
    real(kind_phys) :: oso4_num, odst_num, osoot_num
    real(kind_phys) :: so4_num_st_cr_tot
    real(kind_phys) :: dso4_num, ramp

    ! Aerosol transfer working arrays (allocated only for MAM/CARMA)
    real(kind_phys), allocatable :: size_wght(:,:,:,:)  ! (ncol,pver,nbins,nmaxspc)
    real(kind_phys), allocatable :: amb_num_bins(:,:,:) ! (ncol,pver,nbins)
    real(kind_phys), pointer :: num_col(:,:)
    real(kind_phys), pointer :: amb_mmr(:,:)
    real(kind_phys), pointer :: cld_mmr(:,:)

    character(len=32) :: spectype
    real(kind_phys) :: wght
    real(kind_phys) :: delmmr, delmmr_sum
    real(kind_phys) :: delnum, delnum_sum

    real(kind_phys), parameter :: per_cm3 = 1.e-6_kind_phys  ! m-3 to cm-3

    !-------------------------------------------------------------------------------

    errmsg = ''
    errflg = 0

    nullify(aeroprops)
    nullify(aerostate)

    ! Initialize all outputs
    naai(:,:)     = 0._kind_phys
    naai_hom(:,:) = 0._kind_phys
    nihf(:,:)     = 0._kind_phys
    niimm(:,:)    = 0._kind_phys
    nidep(:,:)    = 0._kind_phys
    nimey(:,:)    = 0._kind_phys
    regm(:,:)     = 0._kind_phys
    trop_pd(:,:)  = 0._kind_phys
    fhom(:,:)     = 0._kind_phys
    wice(:,:)     = 0._kind_phys
    weff(:,:)     = 0._kind_phys
    ptend_q(:,:,:) = 0._kind_phys

    INnso4(:,:)   = 0.0_kind_phys
    INnbc(:,:)    = 0.0_kind_phys
    INndust(:,:)  = 0.0_kind_phys
    INondust(:,:) = 0.0_kind_phys
    INhet(:,:)    = 0.0_kind_phys
    INhom(:,:)    = 0.0_kind_phys
    INFrehom(:,:) = 0.0_kind_phys
    INFreIN(:,:)  = 0.0_kind_phys

    !---------------------------------------------------------------------------
    ! Prepare abstract aerosol properties and state object.
    !---------------------------------------------------------------------------
    if (iaermod_selected_ > 0) then
      aeroprops => aerosol_instances_get_props(iaermod_selected_, list_idx=0)
      aerostate => aerosol_instances_get_state(iaermod_selected_, list_idx=0)
    end if

    ! Allocate aerosol transfer arrays if needed
    if (associated(aeroprops) .and. nbins_ > 0) then
      allocate(size_wght(ncol, pver, nbins_, nmaxspc_))
      allocate(amb_num_bins(ncol, pver, nbins_))
    end if

    !---------------------------------------------------------------------------
    ! Prepare derived input quantities
    !---------------------------------------------------------------------------
    rho(:ncol,:) = pmid(:ncol,:) / (rair * t(:ncol,:))

    ! Compute subgrid saturation factor from tropopause level
    subgrid(:,:) = 0._kind_phys
    do k = top_lev, pver
      do i = 1, ncol
        trop_pd(i, tropLev_chem(i)) = 1._kind_phys

        if (k <= tropLev_chem(i)) then
          if (nucleate_ice_subgrid_strat_ == -1._kind_phys) then
            subgrid(i, k) = 1._kind_phys / qsatfac(i, k)
          else
            subgrid(i, k) = nucleate_ice_subgrid_strat_
          end if
        else
          if (nucleate_ice_subgrid_ == -1._kind_phys) then
            subgrid(i, k) = 1._kind_phys / qsatfac(i, k)
          else
            subgrid(i, k) = nucleate_ice_subgrid_
          end if
        end if
      end do
    end do
    subgrid_diag(:ncol,:pver) = subgrid(:ncol,:pver)

    ! Compute relative humidity and ice cloud fraction
    do k = top_lev, pver
      call qsat_water(t(1:ncol,k), pmid(1:ncol,k), &
                      es(1:ncol), qs(1:ncol), ncol, gam=gammas(1:ncol))
      do i = 1, ncol
        relhum(i,k) = qv(i,k) / qs(i)
        icldm(i,k) = max(aist(i,k), mincld_)
      end do
    end do

    ! Collect number densities [# cm-3] for dust, sulfate, and soot.
    ! Unified for all aerosol models via the abstract interface.
    dust_num_col = 0._kind_phys
    sulf_num_col = 0._kind_phys
    sulf_num_tot_col = 0._kind_phys
    soot_num_col = 0._kind_phys

    if (associated(aeroprops) .and. associated(aerostate)) then
      call aerostate%nuclice_get_numdens( aeroprops, use_preexisting_ice_, &
                                          ncol, pver, rho, &
                                          dust_num_col, sulf_num_col, &
                                          soot_num_col, sulf_num_tot_col )

      ! Prepare per-bin ambient numbers and size weights for aerosol transfer
      do m = 1, nbins_
        call aerostate%get_ambient_num(m, num_col)
        amb_num_bins(:ncol,:,m) = num_col(:ncol,:)
        do l = 1, aeroprops%nspecies(m)
          call aeroprops%species_type(m, l, spectype)
          call aerostate%icenuc_size_wght(m, ncol, pver, spectype, &
                                          use_preexisting_ice_, size_wght(:,:,m,l))
        end do
      end do
    end if

    !---------------------------------------------------------------------------
    ! Column loop: call nucleati kernel
    !---------------------------------------------------------------------------
    do k = top_lev, pver
      do i = 1, ncol
        so4_num_st_cr_tot = 0._kind_phys
        if (t(i,k) < tmelt - 5._kind_phys) then

          ! Set aerosol number for so4, soot, and dust [#/cm3]
          so4_num = sulf_num_col(i,k)
          dst_num = dust_num_col(i,k)
          so4_num_st_cr_tot = sulf_num_tot_col(i,k)

          ! *** Turn off soot nucleation ***
          ! this mod is reproduced verbatim from CAM as of cam6_4_162
          ! (hplin, 4/7/26)
          soot_num = 0.0_kind_phys

          if (use_nucleati_tendencies_) then
            ! PUMAS v1.21+ path (CAM7+)
            call nucleati( &
              wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k), &
              qc(i,k), qi(i,k), ni(i,k), rho(i,k), &
              so4_num, dst_num, soot_num, subgrid(i,k), &
              naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
              wice(i,k), weff(i,k), fhom(i,k), regm(i,k), &
              oso4_num, odst_num, osoot_num, &
              call_frm_zm_in=.false., add_preexisting_ice_in=.false.)
          else
            ! Legacy path
            ! (no optional args - defaults to: call_frm_zm = .false.
            !                       and   add_preexisting_ice = .true)
            ! this path is unsupported in CAM-SIMA (see below) (hplin, 4/7/26)
            call nucleati( &
              wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k), &
              qc(i,k), qi(i,k), ni(i,k), rho(i,k), &
              so4_num, dst_num, soot_num, subgrid(i,k), &
              naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
              wice(i,k), weff(i,k), fhom(i,k), regm(i,k), &
              oso4_num, odst_num, osoot_num)
          end if

          !---------------------------------------------------------------------
          ! MAM / CARMA only:
          ! Move aerosol used for nucleation from interstitial to cloud-borne,
          ! otherwise the same coarse mode aerosols will be available again in
          ! the next timestep and will suppress homogeneous freezing.
          ! Replaces nucleate_ice_cam L620-692.
          !---------------------------------------------------------------------
          if (clim_modal_carma_ .and. use_preexisting_ice_) then

            do m = 1, nbins_

              if (aeroprops%icenuc_updates_num(m)) then

                if (amb_num_bins(i,k,m) > 0._kind_phys) then
                  delmmr_sum = 0._kind_phys
                  delnum_sum = 0._kind_phys

                  ! iterate over the species within the bin
                  do l = 1, aeroprops%nspecies(m)
                    if (aeroprops%icenuc_updates_mmr(m, l)) then

                      call aeroprops%species_type(m, l, spectype)

                      wght = size_wght(i,k,m,l)

                      if (wght > 0._kind_phys) then

                        idxtmp = aer_cnst_idx_(m, l)

                        call aerostate%get_ambient_mmr(species_ndx=l, bin_ndx=m, mmr=amb_mmr)
                        call aerostate%get_cldbrne_mmr(species_ndx=l, bin_ndx=m, mmr=cld_mmr)

                        ! determine change in aerosol mass
                        delmmr = 0._kind_phys
                        delnum = 0._kind_phys
                        if (trim(spectype) == 'dust') then
                          if (dst_num > 0._kind_phys) then
                            delmmr = (odst_num / dst_num) * icldm(i,k) * amb_mmr(i,k) * wght
                            delnum = (odst_num * icldm(i,k)) / rho(i,k) / per_cm3
                          end if
                        else if (trim(spectype) == 'sulfate') then
                          if (so4_num > 0._kind_phys) then
                            delmmr = (oso4_num / so4_num) * icldm(i,k) * amb_mmr(i,k) * wght
                            delnum = (oso4_num * icldm(i,k)) / rho(i,k) / per_cm3
                          end if
                        end if

                        if (idxtmp > 0) then
                          ! constituent tendency for transported species
                          ptend_q(i,k,idxtmp) = -delmmr / dtime
                        else
                          ! apply change of mass to not-transported species directly
                          amb_mmr(i,k) = amb_mmr(i,k) - delmmr
                        end if
                        cld_mmr(i,k) = cld_mmr(i,k) + delmmr

                        delmmr_sum = delmmr_sum + delmmr
                        delnum_sum = delnum_sum + delnum
                      end if
                    end if
                  end do

                  idxtmp = aer_cnst_idx_(m, 0)

                  ! update aerosol state bin and tendency for grid box i,k
                  call aerostate%update_bin(m, i, k, delmmr_sum, delnum_sum, &
                                            idxtmp, dtime, ptend_q)

                end if

              end if
            end do

          end if ! clim_modal_carma_ .and. use_preexisting_ice_

          !---------------------------------------------------------------------
          ! Polar strat SO4 enhancement
          !---------------------------------------------------------------------
          ! Liu&Penner does not generate enough nucleation in the polar winter
          ! stratosphere, which affects surface area density, dehydration and
          ! ozone chemistry. Part of this is that there are a larger number of
          ! particles in the accumulation mode than in the Aitken mode. In volcanic
          ! periods, the coarse mode may also be important. As a short
          ! term work around, include the accumulation and coarse mode particles
          ! and assume a larger fraction of the sulfates nucleate in the polar
          ! stratosphere.
          !
          ! Do not include the tropopause level, as stratospheric aerosols
          ! only exist above the tropopause level.
          !
          ! NOTE: This may still not represent the proper particles that
          ! participate in nucleation, because it doesn't include STS and NAT
          ! particles. It may not represent the proper saturation threshold for
          ! nucleation, and wsubi from CLUBB is probably not representative of
          ! wave driven varaibility in the polar stratosphere.
          if (nucleate_ice_use_troplev_ .and. clim_modal_carma_) then
            if ((k < tropLev_chem(i)) .and. &
                (nucleate_ice_strat_ > 0._kind_phys) .and. &
                (oso4_num > 0._kind_phys)) then
              dso4_num = max(0._kind_phys, &
                (nucleate_ice_strat_ * so4_num_st_cr_tot - oso4_num) &
                * 1e6_kind_phys / rho(i,k))
              naai(i,k) = naai(i,k) + dso4_num
              nihf(i,k) = nihf(i,k) + dso4_num
            end if
          else
            ! Fallback: pressure/latitude heuristic (maintains backwards compatibility)
            if (pmid(i,k) <= 12500._kind_phys .and. &
                pmid(i,k) > 100._kind_phys .and. &
                abs(lat(i)) >= 60._kind_phys * pi / 180._kind_phys) then
              ramp = 1._kind_phys - min(1._kind_phys, &
                max(0._kind_phys, (pmid(i,k) - 10000._kind_phys) / 2500._kind_phys))

              if (oso4_num > 0._kind_phys) then
                dso4_num = (max(oso4_num, ramp * nucleate_ice_strat_ * so4_num) &
                  - oso4_num) * 1e6_kind_phys / rho(i,k)
                naai(i,k) = naai(i,k) + dso4_num
                nihf(i,k) = nihf(i,k) + dso4_num
              end if
            end if
          end if

          !---------------------------------------------------------------------
          ! Convert to tendencies for PUMAS (CAM7+)
          !---------------------------------------------------------------------
          if (use_nucleati_tendencies_) then
            ! ^^ this flag was previously CAM7+
            ! PUMAS v1.21+ path: convert to tendencies [1/kg/s] and [1/m3/s]

            naai_hom(i,k) = nihf(i,k) / dtime
            naai(i,k) = naai(i,k) / dtime

            ! output activated ice (convert from # kg-1 to # m-3 s-1)
            nihf(i,k)  = nihf(i,k)  * rho(i,k) / dtime
            niimm(i,k) = niimm(i,k) * rho(i,k) / dtime
            nidep(i,k) = nidep(i,k) * rho(i,k) / dtime
            nimey(i,k) = nimey(i,k) * rho(i,k) / dtime

            if (use_preexisting_ice_) then
              INnso4(i,k)  = so4_num  * 1e6_kind_phys / dtime  ! # cm-3 -> # m-3 s-1
              INnbc(i,k)   = soot_num * 1e6_kind_phys / dtime
              INndust(i,k) = dst_num  * 1e6_kind_phys / dtime
              INondust(i,k)= odst_num * 1e6_kind_phys / dtime
              INFreIN(i,k) = 1.0_kind_phys        ! 1, ice nucleation occurred
              INhet(i,k) = (niimm(i,k) + nidep(i,k))  ! # m-3 s-1, nimey not in cirrus
              INhom(i,k) = nihf(i,k)                   ! # m-3 s-1
              if (INhom(i,k) > 1e3_kind_phys) then     ! > 1/L
                INFrehom(i,k) = 1.0_kind_phys           ! 1, hom freezing occurred
              end if

              ! exclude no ice nucleation
              if ((INFrehom(i,k) < 0.5_kind_phys) .and. &
                  (INhet(i,k)    < 1.0_kind_phys)) then
                INnso4(i,k)   = 0.0_kind_phys
                INnbc(i,k)    = 0.0_kind_phys
                INndust(i,k)  = 0.0_kind_phys
                INondust(i,k) = 0.0_kind_phys
                INFreIN(i,k)  = 0.0_kind_phys
                INhet(i,k)    = 0.0_kind_phys
                INhom(i,k)    = 0.0_kind_phys
                INFrehom(i,k) = 0.0_kind_phys
                wice(i,k)     = 0.0_kind_phys
                weff(i,k)     = 0.0_kind_phys
                fhom(i,k)     = 0.0_kind_phys
              end if
            end if
          else
            ! Legacy path (CAM5/6): unsupported in CAM-SIMA CCPP.
            ! Kept commented for reference. If needed, uncomment and
            ! change output units from tendencies to concentrations.
            !
            ! naai_hom(i,k) = nihf(i,k)
            ! nihf(i,k)  = nihf(i,k)  * rho(i,k)
            ! niimm(i,k) = niimm(i,k) * rho(i,k)
            ! nidep(i,k) = nidep(i,k) * rho(i,k)
            ! nimey(i,k) = nimey(i,k) * rho(i,k)

            errflg = 1
            errmsg = 'non-CAM7 path for nucleate_ice is unsupported in CCPP'
            return
          end if ! use_nucleati_tendencies_
        end if ! T < tmelt - 5
      end do ! i
    end do ! k

    ! Clean up
    if (allocated(size_wght)) deallocate(size_wght)
    if (allocated(amb_num_bins)) deallocate(amb_num_bins)

  end subroutine nucleate_ice_ccpp_run

end module nucleate_ice_ccpp
