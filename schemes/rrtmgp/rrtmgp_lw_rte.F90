!> \file rrtmgp_lw_rte.F90
!!

!> This module contains the call to the RRTMGP-LW radiation routine
module rrtmgp_lw_rte
  implicit none
  private

  public rrtmgp_lw_rte_run
contains

!> \section arg_table_rrtmgp_lw_rte_run Argument Table
!! \htmlinclude rrtmgp_lw_rte_run.html
!!
   subroutine rrtmgp_lw_rte_run(doLWrad, doLWclrsky, doGP_lwscat, use_LW_jacobian, use_LW_optimal_angles,   &
                                 nGauss_angles, lw_optical_props_clrsky, lw_optical_props_clouds,           &
                                 sources, sfc_emiss_byband, lw_gas_props, aerlw, fluxlwUP_jac, lw_Ds,       &
                                 flux_clrsky, flux_allsky, errmsg, errflg)
    use machine,                  only: kind_phys
    use mo_rte_lw,                only: rte_lw
    use ccpp_fluxes_byband,       only: ty_fluxes_byband_ccpp
    use ccpp_optical_props,       only: ty_optical_props_1scl_ccpp
    use ccpp_fluxes,              only: ty_fluxes_broadband_ccpp
    use ccpp_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp_ccpp
    use ccpp_source_functions,    only: ty_source_func_lw_ccpp
    use radiation_tools,          only: check_error_msg

    ! Inputs
    logical, intent(in) :: doLWrad               !< Flag to perform longwave calculation
    logical, intent(in) :: doLWclrsky            !< Flag to compute clear-sky fluxes
    logical, intent(in) :: doGP_lwscat           !< Flag to include scattering in clouds
    logical, intent(in) :: use_LW_jacobian       !< Flag to compute Jacobian
    logical, intent(in) :: use_LW_optimal_angles !< Flag to compute and use optimal angles

    integer, target, intent(in) :: nGauss_angles !< Number of gaussian quadrature angles used

    real(kind_phys), dimension(:,:),   intent(in) :: sfc_emiss_byband           !< Surface emissivity by band
    class(ty_source_func_lw_ccpp),     intent(in) :: sources                    !< Longwave sources object

    ! Outputs
    real(kind_phys), dimension(:,:), target, intent(inout) :: fluxlwUP_jac      !< Surface temperature flux Jacobian [W m-2 K-1]
    class(ty_fluxes_byband_ccpp),      intent(inout) :: flux_allsky             !< All-sky flux [W m-2]
    class(ty_fluxes_broadband_ccpp),   intent(inout) :: flux_clrsky             !< Clear-sky flux [W m-2]
    class(ty_optical_props_1scl_ccpp), intent(inout) :: aerlw                   !< Aerosol optical properties object
    class(ty_optical_props_1scl_ccpp), intent(inout) :: lw_optical_props_clrsky !< Clear-sky optical properties object
    class(ty_optical_props_1scl_ccpp), intent(inout) :: lw_optical_props_clouds !< Cloud optical properties object

    class(ty_gas_optics_rrtmgp_ccpp),  intent(inout) :: lw_gas_props            !< Gas optical properties object

    real(kind_phys), dimension(:,:), target, intent(out) :: lw_Ds               !< 1/cos of transport angle per column, g-point
    character(len=512),intent(out) :: errmsg                                    !< CCPP error message
    integer,           intent(out) :: errflg                                    !< CCPP error flag

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    !$acc data copyin(lw_optical_props_clrsky%optical_props,lw_optical_props_clrsky%optical_props%tau,   &
    !$acc             aerlw%optical_props,aerlw%optical_props%tau,          &
    !$acc             lw_optical_props_clouds%optical_props, lw_optical_props_clouds%optical_props%tau,        &
    !$acc             sources%sources,sources%sources%lay_source,     &
    !$acc             sources%sources%sfc_source,     &
    !$acc             sources%sources%lev_source,     &
    !$acc             sources%sources%sfc_source_jac, &
    !$acc             sfc_emiss_byband)                          &
    !$acc        copy(flux_clrsky%fluxes, flux_clrsky%fluxes%flux_net, flux_clrsky%fluxes%flux_up, &
    !$acc             flux_clrsky%fluxes%flux_dn, flux_allsky%fluxes, flux_allsky%fluxes%flux_net,  &
    !$acc             flux_allsky%fluxes%flux_up, flux_allsky%fluxes%flux_dn,    &
    !$acc             lw_Ds)

    ! ###################################################################################
    !
    ! Compute clear-sky fluxes (gaseous+aerosol) (optional)
    !
    ! ###################################################################################
    ! Increment
    errmsg = aerlw%optical_props%increment(lw_optical_props_clrsky%optical_props)
    call check_error_msg('rrtmgp_lw_rte_increment_aerosol_to_clrsky', errmsg)
    if (len_trim(errmsg) /= 0) then
        errflg = 1
        return
    end if

    ! Call RTE solver
    if (doLWclrsky) then
       if (use_lw_optimal_angles) then
          errmsg = lw_gas_props%gas_props%compute_optimal_angles(lw_optical_props_clrsky%optical_props,lw_Ds)
          call check_error_msg('rrtmgp_lw_rte_opt_angle', errmsg)
          if (len_trim(errmsg) /= 0) then
             errflg = 1
             return
          end if
          if (nGauss_angles > 1) then
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_clrsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles,        & ! IN  - Number of angles in Gaussian quadrature
                  lw_Ds = lw_Ds)                           ! IN  - 1/cos of transport angle per column and g-point
          else
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_clrsky%fluxes,                    & ! OUT - Fluxes
                  lw_Ds = lw_Ds)                           ! IN  - 1/cos of transport angle per column and g-point
          end if
       else
          if (nGauss_angles > 1) then
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_clrsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles)          ! IN  - Number of angles in Gaussian quadrature
          else
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_clrsky%fluxes)                      ! OUT - Fluxes
          end if
       end if
       call check_error_msg('rrtmgp_lw_rte_lw_rte_clrsky', errmsg)
       if (len_trim(errmsg) /= 0) then
          errflg = 1
          return
       end if
    end if

    ! ###################################################################################
    !
    ! All-sky fluxes (clear-sky + clouds + precipitation)
    ! *Note* CCPP does not allow for polymorphic types, they are ambiguous to the CCPP
    ! framework. rte-rrtmgp uses polymorphic types extensively, for example, querying the
    ! type to determine physics configuration/pathway/etc...
    !
    ! The logic in the code below is to satisfy the polymorphishm in the rte-rrtmgp code.
    ! The rte-rrtmgp "increment" procedures are utilized to provide the correct type to the
    ! rte solver (rte_lw). Rte_lw queries the type to determine if scattering is to be
    ! included in the calculation. The increment procedures are called so that the correct
    ! optical properties are inherited.
    ! 
    ! ###################################################################################

    ! Include LW cloud-scattering?
    if (doGP_lwscat) then 
       ! Increment
       errmsg = lw_optical_props_clrsky%optical_props%increment(lw_optical_props_clouds%optical_props)
       call check_error_msg('rrtmgp_lw_rte_increment_clrsky_to_clouds', errmsg)
       if (len_trim(errmsg) /= 0) then
           errflg = 1
           return
       end if
       if (use_LW_jacobian) then
          if (nGauss_angles > 1) then
             errmsg = rte_lw(           &
                  lw_optical_props_clouds%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles,        & ! IN  - Number of angles in Gaussian quadrature
                  flux_up_Jac    = fluxlwUP_jac)           ! OUT - surface temperature flux (upward) Jacobian (W m-2 K-1)
          else
             errmsg = rte_lw(           &
                  lw_optical_props_clouds%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  flux_up_Jac    = fluxlwUP_jac)           ! OUT - surface temperature flux (upward) Jacobian (W m-2 K-1)
          end if
       else
          if (nGauss_angles > 1) then
             errmsg = rte_lw(           &
                  lw_optical_props_clouds%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles)          ! IN  - Number of angles in Gaussian quadrature
          else
             errmsg = rte_lw(           &
                  lw_optical_props_clouds%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes)                      ! OUT - Fluxes
          end if
       end if
    ! No scattering in LW clouds.   
    else
       ! Increment
       errmsg = lw_optical_props_clouds%optical_props%increment(lw_optical_props_clrsky%optical_props)
       call check_error_msg('rrtmgp_lw_rte_increment_clouds_to_clrsky', errmsg)
       if (len_trim(errmsg) /= 0) then
           errflg = 1
           return
       end if
 
       if (use_LW_jacobian) then
          if (nGauss_angles > 1) then
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles,        & ! IN  - Number of angles in Gaussian quadrature
                  flux_up_Jac    = fluxlwUP_jac)           ! OUT - surface temperature flux (upward) Jacobian (W m-2 K-1)
          else
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  flux_up_Jac    = fluxlwUP_jac)           ! OUT - surface temperature flux (upward) Jacobian (W m-2 K-1)
          end if
       else
          if (nGauss_angles > 1) then
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes,                    & ! OUT - Fluxes
                  n_gauss_angles = nGauss_angles)          ! IN  - Number of angles in Gaussian quadrature
          else
             errmsg = rte_lw(           &
                  lw_optical_props_clrsky%optical_props, & ! IN  - optical-properties
                  sources%sources,                       & ! IN  - source function
                  sfc_emiss_byband,                      & ! IN  - surface emissivity in each LW band
                  flux_allsky%fluxes)                      ! OUT - Fluxes
          end if
       end if
    end if
    call check_error_msg('rrtmgp_lw_rte_lw_rte_allsky', errmsg)
    if (len_trim(errmsg) /= 0) then
       errflg = 1
    end if
    !$acc end data

  end subroutine rrtmgp_lw_rte_run
end module rrtmgp_lw_rte
