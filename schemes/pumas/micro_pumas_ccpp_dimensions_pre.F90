module micro_pumas_ccpp_dimensions_pre

  implicit none

  use pumas_kinds,       only: pumas_r8=>kind_r8

contains

  !> \section arg_table_micro_pumas_ccpp_dimensions_pre_init Argument Table
  !! \htmlinclude micro_pumas_ccpp_dimensions_pre_init.html
  subroutine micro_pumas_ccpp_dimensions_pre_init( ncol, nlev, nlevp1,              &
                             trop_cloud_top_lev, micro_ncol, micro_nlev,            &
                             micro_nlevp1, errmsg, errcode)

    !Host model dimensions/parameters:
    integer,         intent(in) :: ncol
    integer,         intent(in) :: nlev
    integer,         intent(in) :: nlevp1
    integer,         intent(in) :: trop_cloud_top_lev  !Index of the top model level for which
                                                       !cloud physics is applied (1 to nlev)

    !PUMAS dimensions/parameters:
    integer,         intent(out) :: micro_ncol         !Number of horizontal microphysics columns (count)
    integer,         intent(out) :: micro_nlev         !Number of microphysics vertical layers (count)
    integer,         intent(out) :: micro_nlevp1       !Number of microphysics vertical interfaces (count)

    !CCPP error handling:
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    !Initialize error message and error code:
    errmsg  = ''
    errcode = 0

    !Set PUMAS (cloud microphysics) dimensions:
    micro_ncol   = ncol
    micro_nlev   = nlev-trop_cloud_top_lev+1
    micro_nlevp1 = micro_nlev + 1

  end subroutine micro_pumas_ccpp_dimensions_pre_init

  !> \section arg_table_micro_pumas_ccpp_dimensions_pre_run Argument Table
  !! \htmlinclude micro_pumas_ccpp_dimensions_pre_run.html
  subroutine micro_pumas_ccpp_dimensions_pre_run(ncol, nlev, nlevp1,                &
                             micro_ncol, micro_nlev, micro_nlevp1, micro_dust_nbins,&
                             airT_in, pumas_airT, airq_in, pumas_airq,              &
                             cldliq_in, pumas_cldliq,                               &
                             cldice_in, pumas_cldice,                               &
                             micro_mg_num_steps, dtime, pumas_timestep,             &
                             numliq_in, pumas_numliq,                               &
                             numice_in, pumas_numice,                               &
                             rainliq_in, pumas_rainliq,                             &
                             snowice_in, pumas_snowice,                             &
                             numrain_in, pumas_numrain,                             &
                             numsnow_in, pumas_numsnow,                             &
                             graupice_in, pumas_graupice,                           &
                             numgraup_in, pumas_numgraup,                           &
                             relvar_in, pumas_relvar,                               &
                             accre_enhan_in, pumas_accre_enhan,           &
                             pmid_in, pumas_pmid,                                   &
                             pdel_in, pumas_pdel,                                   &
                             pint_in, pumas_pint,                                   &
                             strat_cldfrc_in, pumas_strat_cldfrc,                   &
                             strat_liq_cldfrc_in, pumas_strat_liq_cldfrc,           &
                             strat_ice_cldfrc_in, pumas_strat_ice_cldfrc,           &
                             qsatfac_in, pumas_qsatfac,                             &
                             naai_in, pumas_naai,                                   &
                             npccn_in, pumas_npccn,                                 &
                             rndst_in, pumas_rndst,                                 &
                             nacon_in, pumas_nacon,                                 &
                             snowice_tend_external_in, pumas_snowice_tend_external, &
                             numsnow_tend_external_in, pumas_numsnow_tend_external, &
                             effi_external_in, pumas_effi_external,                 &
                             frzimm_in, pumas_frzimm,                               &
                             frzcnt_in, pumas_frzcnt,                               &
                             frzdep_in, pumas_frzdep,                               &
                             errmsg, errcode)

    !External dependencies:
    use ccpp_kinds,        only: kind_phys

    !Host model dimensions/parameters:
    integer,         intent(in) :: ncol
    integer,         intent(in) :: nlev
    integer,         intent(in) :: nlevp1

    !PUMAS dimensions/parameters:
    integer,         intent(in) :: micro_ncol          !Number of horizontal microphysics columns (count)
    integer,         intent(in) :: micro_nlev          !Number of microphysics vertical layers (count)
    integer,         intent(in) :: micro_nlevp1        !Number of microphysics vertical interfaces (count)
    integer,         intent(in) :: micro_dust_nbins    !Number of dust bins

    ! Air temperature (K)
    real(kind_phys), intent(in)  :: airT_in(:, :)
    real(pumas_r8), intent(out) :: pumas_airT(:, :)
    ! Water vapor mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: airq_in(:, :)
    real(pumas_r8), intent(out) :: pumas_airq(:, :)
    ! Cloud liquid water mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: cldliq_in(:, :)
    real(pumas_r8), intent(out) :: pumas_cldliq(:, :)
    ! Cloud ice mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: cldice_in(:, :)
    real(pumas_r8), intent(out) :: pumas_cldice(:, :)
    ! Substepping in microphysics variables
    integer, intent(in)  :: micro_mg_num_steps
    real(kind_phys), intent(in)  :: dtime
    real(pumas_r8), intent(out) :: pumas_timestep
    ! Mass number concentration of cloud liquid water wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numliq_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numliq(:, :)
    ! Mass number concentration of cloud ice wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numice_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numice(:, :)
    ! Rain mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: rainliq_in(:, :)
    real(pumas_r8), intent(out) :: pumas_rainliq(:, :)
    ! Snow mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: snowice_in(:, :)
    real(pumas_r8), intent(out) :: pumas_snowice(:, :)
    ! Mass number concentration of rain wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numrain_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numrain(:, :)
    ! Mass number concentration of snow wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numsnow_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numsnow(:, :)
    ! Graupel mixing ratio wrt moist air and condensed water (kg kg-1)
    real(kind_phys), intent(in)  :: graupice_in(:, :)
    real(pumas_r8), intent(out) :: pumas_graupice(:, :)
    ! Mass number concentration of graupel wrt moist air and condensed water (kg-1)
    real(kind_phys), intent(in)  :: numgraup_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numgraup(:, :)
    ! Relative variance of cloud water (1)
    real(kind_phys), intent(in)  :: relvar_in(:, :)
    real(pumas_r8), intent(out) :: pumas_relvar(:, :)
    ! Accretion enhancement factor (1)
    real(kind_phys), intent(in)  :: accre_enhan_in(:, :)
    real(pumas_r8),  intent(out) :: pumas_accre_enhan(:, :)
    ! Air pressure (Pa)
    real(kind_phys), intent(in)  :: pmid_in(:, :)
    real(pumas_r8), intent(out) :: pumas_pmid(:, :)
    ! Air pressure thickness (Pa)
    real(kind_phys), intent(in)  :: pdel_in(:, :)
    real(pumas_r8), intent(out) :: pumas_pdel(:, :)
    ! Air pressure at interfaces (Pa)
    real(kind_phys), intent(in)  :: pint_in(:, :)
    real(pumas_r8), intent(out) :: pumas_pint(:, :)
    ! Stratiform cloud area fraction (fraction)
    real(kind_phys), intent(in)  :: strat_cldfrc_in(:, :)
    real(pumas_r8), intent(out) :: pumas_strat_cldfrc(:, :)
    ! Stratiform cloud liquid area fraction (fraction)
    real(kind_phys), intent(in)  :: strat_liq_cldfrc_in(:, :)
    real(pumas_r8), intent(out) :: pumas_strat_liq_cldfrc(:, :)
    ! Stratiform cloud ice area fraction (fraction)
    real(kind_phys), intent(in)  :: strat_ice_cldfrc_in(:, :)
    real(pumas_r8), intent(out) :: pumas_strat_ice_cldfrc(:, :)
    ! Subgrid cloud water saturation scaling factor (1)
    real(kind_phys), intent(in)  :: qsatfac_in(:, :)
    real(pumas_r8), intent(out) :: pumas_qsatfac(:, :)
    ! Tendency of activated ice nuclei mass number concentration (kg-1 s-1)
    real(kind_phys), intent(in)  :: naai_in(:, :)
    real(pumas_r8), intent(out) :: pumas_naai(:, :)
    ! Tendency of activated cloud condensation nuclei mass number concentration (kg-1 s-1)
    real(kind_phys), intent(in)  :: npccn_in(:, :)
    real(pumas_r8), intent(out) :: pumas_npccn(:, :)
    ! Dust radii by size bin  (m)
    real(kind_phys), intent(in)  :: rndst_in(:, :, :)
    real(pumas_r8), intent(out) :: pumas_rndst(:, micro_nlev, :)
    ! Dust number concentration by size bin (m-3)
    real(kind_phys), intent(in)  :: nacon_in(:, :, :)
    real(pumas_r8), intent(out) :: pumas_nacon(:, micro_nlev, :)
    ! Tendency of snow mixing ratio wrt moist air and condensed water from external microphysics (kg kg-1 s-1)
    real(kind_phys), intent(in)  :: snowice_tend_external_in(:, :)
    real(pumas_r8), intent(out) :: pumas_snowice_tend_external(:, :)
    ! Tendency of mass number concentration of snow wrt moist air and condensed water from external microphysics (kg-1 s-1)
    real(kind_phys), intent(in)  :: numsnow_tend_external_in(:, :)
    real(pumas_r8), intent(out) :: pumas_numsnow_tend_external(:, :)
    ! Effective radius of stratiform cloud ice particle from external microphysics (m)
    real(kind_phys), intent(in)  :: effi_external_in(:, :)
    real(pumas_r8), intent(out) :: pumas_effi_external(:, :)
    ! Tendency of cloud liquid droplet number concentration due to immersion freezing (cm-3)
    real(kind_phys), intent(in)  :: frzimm_in(:, :)
    real(pumas_r8), intent(out) :: pumas_frzimm(:, :)
    ! Tendency of cloud liquid droplet number concentration due to contact freezing (cm-3)
    real(kind_phys), intent(in)  :: frzcnt_in(:, :)
    real(pumas_r8), intent(out) :: pumas_frzcnt(:, :)
    ! Tendency of cloud ice number concentration due to deposition nucleation (cm-3)
    real(kind_phys), intent(in)  :: frzdep_in(:, :)
    real(pumas_r8), intent(out) :: pumas_frzdep(:, :)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    !Initialize error message and error code:
    errmsg  = ''
    errcode = 0

    pumas_timestep = dtime/micro_mg_num_steps

!+ IH
! For now we just use ncols = micro_ncol, but we need to constrain the vertical extent for the microphysical fields.
! Therefore micro_xxx(:ncol,:) = xxx(:,::)
!- IH
    pumas_airT(:ncol,:) = real(airT_in(:,::), pumas_r8)
    pumas_airq(:ncol,:) = real(airq_in(:,::), pumas_r8)
    pumas_cldliq(:ncol,:) = real(cldliq_in(:,::), pumas_r8)
    pumas_cldice(:ncol,:) = real(cldice_in(:,::), pumas_r8)
    pumas_numliq(:ncol,:) = real(numliq_in(:,::), pumas_r8)
    pumas_numice(:ncol,:) = real(numice_in(:,::), pumas_r8)
    pumas_rainliq(:ncol,:) = real(rainliq_in(:,::), pumas_r8)
    pumas_snowice(:ncol,:) = real(snowice_in(:,::), pumas_r8)
    pumas_numrain(:ncol,:) = real(numrain_in(:,::), pumas_r8)
    pumas_numsnow(:ncol,:) = real(numsnow_in(:,::), pumas_r8)
    pumas_graupice(:ncol,:) = real(graupice_in(:,::), pumas_r8)
    pumas_numgraup(:ncol,:) = real(numgraup_in(:,::), pumas_r8)
    pumas_relvar(:ncol,:) = real(relvar_in(:,::), pumas_r8)
    pumas_accre_enhan(:ncol,:) = real(accre_enhan_in(:,::), pumas_r8)
    pumas_pmid(:ncol,:) = real(pmid_in(:,::), pumas_r8)
    pumas_pdel(:ncol,:) = real(pdel_in(:,::), pumas_r8)
    pumas_pint(:ncol,:) = real(pint_in(:,:micro_nlevp1), pumas_r8)
    pumas_strat_cldfrc(:ncol,:) = real(strat_cldfrc_in(:,::), pumas_r8)
    pumas_strat_liq_cldfrc(:ncol,:) = real(strat_liq_cldfrc_in(:,::), pumas_r8)
    pumas_strat_ice_cldfrc(:ncol,:) = real(strat_ice_cldfrc_in(:,::), pumas_r8)
    pumas_qsatfac(:ncol,:) = real(qsatfac_in(:,::), pumas_r8)
    pumas_naai(:ncol,:) = real(naai_in(:,::), pumas_r8)
    pumas_npccn(:ncol,:) = real(npccn_in(:,::), pumas_r8)
    pumas_rndst(:ncol,:,:) = real(rndst_in(:,:micro_nlev,:), pumas_r8)
    pumas_nacon(:ncol,:,:) = real(nacon_in(:,:micro_nlev,:), pumas_r8)
    pumas_snowice_tend_external(:ncol,:) = real(snowice_tend_external_in(:,::), pumas_r8)
    pumas_numsnow_tend_external(:ncol,:) = real(numsnow_tend_external_in(:,::), pumas_r8)
    pumas_effi_external(:ncol,:) = real(effi_external_in(:,::), pumas_r8)
    pumas_frzcnt(:ncol,:) = real(frzcnt_in(:,::), pumas_r8)
    pumas_frzdep(:ncol,:) = real(frzdep_in(:,::), pumas_r8)


  end subroutine micro_pumas_ccpp_dimensions_pre_run

end module micro_pumas_ccpp_dimensions_pre
