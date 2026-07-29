! Copyright (C) 2024-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

!> Top-level wrapper for MUSICA chemistry components
module musica_ccpp
  use musica_ccpp_micm,     only: micm_register, micm_init, micm_final
  use musica_ccpp_namelist, only: filename_of_tuvx_micm_mapping_configuration
  use musica_ccpp_tuvx,     only: tuvx_register, tuvx_init, tuvx_final
  use musica_ccpp_util,     only: do_tuvx
  use musica_util,          only: index_mappings_t

  implicit none
  private
  save

  public :: musica_ccpp_register, musica_ccpp_init, musica_ccpp_final

contains

  !> \section arg_table_musica_ccpp_register Argument Table
  !! \htmlinclude musica_ccpp_register.html
  subroutine musica_ccpp_register(constituent_props, number_of_micm_rate_parameters, &
                                  errmsg, errcode)
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_properties_t
    use musica_ccpp_micm,              only: rate_parameters_ordering
    use musica_ccpp_namelist,          only: micm_solver_type, filename_of_tuvx_configuration
    use musica_ccpp_species,           only: musica_species_t, register_musica_species
    use musica_ccpp_tuvx_load_species, only: check_tuvx_species_initialization
    use musica_ccpp_util,              only: set_do_tuvx

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    integer,                                          intent(out) :: number_of_micm_rate_parameters
    character(len=512),                               intent(out) :: errmsg
    integer,                                          intent(out) :: errcode

    ! local variables
    type(ccpp_constituent_properties_t), allocatable :: constituent_props_subset(:)
    type(musica_species_t),              allocatable :: micm_species(:)
    type(musica_species_t),              allocatable :: tuvx_species(:)

    call micm_register(micm_solver_type, constituent_props_subset, micm_species, &
                       errmsg, errcode)
    if (errcode /= 0) return
    constituent_props = constituent_props_subset
    deallocate(constituent_props_subset)

    ! Export the number of MICM rate parameters (photolysis, user-defined,
    ! emission, loss, and surface reaction parameters) so the framework can
    ! size the rate parameter array shared by the photolysis and chemistry
    ! run schemes
    number_of_micm_rate_parameters = rate_parameters_ordering%size()
    if (number_of_micm_rate_parameters < 0) then
      errmsg = "[MUSICA Error] Internal error: number_of_micm_rate_parameters < 0"
      errcode = 1
      return
    end if

    if (trim(filename_of_tuvx_configuration) /= "none") then
      call set_do_tuvx(.true.)
      call tuvx_register(micm_species, tuvx_species, constituent_props_subset, &
                         errmsg, errcode)
      if (errcode /= 0) return
      constituent_props = [ constituent_props, constituent_props_subset ]
    else
      call set_do_tuvx(.false.)
      allocate(tuvx_species(0))
    end if

    call register_musica_species(micm_species, tuvx_species)
    ! The TUV-x species indices are only assigned by tuvx_register, so they
    ! can only be validated when TUV-x is enabled
    if (do_tuvx) then
      call check_tuvx_species_initialization(errmsg, errcode)
      if (errcode /= 0) return
    end if

  end subroutine musica_ccpp_register

  !> \section arg_table_musica_ccpp_init Argument Table
  !! \htmlinclude musica_ccpp_init.html
  subroutine musica_ccpp_init(horizontal_dimension, vertical_layer_dimension, &
                              vertical_interface_dimension, &
                              photolysis_wavelength_grid_interfaces, &
                              constituent_props_ptr, molar_mass_dry_air__g_mol, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t, ccpp_constituent_properties_t
    use ccpp_kinds,                only: kind_phys
    use musica_ccpp_micm,          only: rate_parameters_ordering
    use musica_ccpp_namelist,      only: micm_solver_type
    use musica_ccpp_util,          only: has_error_occurred, m_to_nm, set_constants
    use musica_ccpp_species,       only: initialize_musica_species_indices, initialize_molar_mass_array, &
                                         check_initialization, musica_species_t

    integer,                           intent(in)  :: horizontal_dimension                     ! (count)
    integer,                           intent(in)  :: vertical_layer_dimension                 ! (count)
    integer,                           intent(in)  :: vertical_interface_dimension             ! (count)
    real(kind_phys),                   intent(in)  :: photolysis_wavelength_grid_interfaces(:) ! m
    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props_ptr(:)
    real(kind_phys),                   intent(in)  :: molar_mass_dry_air__g_mol                ! g mol-1
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errcode

    ! local variables
    type(ccpp_constituent_properties_t), allocatable :: constituent_props(:)
    type(musica_species_t),              allocatable :: micm_species(:)
    integer                                          :: number_of_grid_cells
    real(kind_phys), dimension(size(photolysis_wavelength_grid_interfaces)) &
                                                     :: photolysis_wavelength_grid_interfaces_nm ! nm

    number_of_grid_cells = horizontal_dimension * vertical_layer_dimension

    call set_constants(molar_mass_dry_air__g_mol * 1.0e-3_kind_phys) ! kg mol-1

    call micm_init(number_of_grid_cells, errmsg, errcode)
    if (errcode /= 0) return

    if (do_tuvx) then
      if (size(photolysis_wavelength_grid_interfaces) < 2) then
        errmsg = "MUSICA: Internal error: invalid photolysis_wavelength_grid_interfaces size."
        errcode = 1
        return
      end if
      photolysis_wavelength_grid_interfaces_nm(:) = &
          photolysis_wavelength_grid_interfaces(:) * m_to_nm
      call tuvx_init(vertical_layer_dimension, vertical_interface_dimension, &
                     photolysis_wavelength_grid_interfaces_nm,               &
                     rate_parameters_ordering, errmsg, errcode)
      if (errcode /= 0) return
    end if

    call initialize_musica_species_indices(constituent_props_ptr, errmsg, errcode)
    if (errcode /= 0) return

    call initialize_molar_mass_array(constituent_props_ptr, errmsg, errcode)
    if (errcode /= 0) return

    call check_initialization(errmsg, errcode)
    if (errcode /= 0) return

  end subroutine musica_ccpp_init

  !> \section arg_table_musica_ccpp_final Argument Table
  !! \htmlinclude musica_ccpp_final.html
  subroutine musica_ccpp_final(errmsg, errcode)
    use musica_ccpp_species, only: cleanup_musica_species
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    call cleanup_musica_species()

    call tuvx_final(errmsg, errcode)
    if (errcode /= 0) return

    call micm_final(errmsg, errcode)
    if (errcode /= 0) return

  end subroutine musica_ccpp_final

end module musica_ccpp