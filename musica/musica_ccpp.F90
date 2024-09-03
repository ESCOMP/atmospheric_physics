!> Top-level wrapper for MUSICA chemistry components
module musica_ccpp
  use musica_ccpp_micm, only : micm_register, micm_init, micm_run, micm_final
  use musica_ccpp_tuvx, only : tuvx_init, tuvx_run, tuvx_final

  implicit none
  private

  public :: musica_ccpp_register, musica_ccpp_init, musica_ccpp_run, musica_ccpp_final

contains

  subroutine musica_ccpp_register(constituents, solver_type, num_grid_cells, errmsg, errcode)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituents(:)
    integer,                                          intent(in)  :: solver_type
    integer,                                          intent(in)  :: num_grid_cells
    character(len=512),                               intent(out) :: errmsg
    integer,                                          intent(out) :: errcode

    call micm_register(constituents, solver_type, num_grid_cells, errmsg, errcode)

  end subroutine musica_ccpp_register

  !> \section arg_table_musica_ccpp_init Argument Table
  !! \htmlinclude musica_ccpp_init.html
  subroutine musica_ccpp_init(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    call tuvx_init(errmsg, errcode)
    call micm_init(errmsg, errcode)

  end subroutine musica_ccpp_init

  !> \section arg_table_musica_ccpp_run Argument Table
  !! \htmlinclude musica_ccpp_run.html
  subroutine musica_ccpp_run(time_step, temperature, pressure, dry_air_density, constituent_props, &
                      constituents, rate_params, height, errmsg, errcode)
    use micm_util,                 only: reshape_into_micm_arr, reshape_into_ccpp_arr
    use micm_util,                 only: convert_to_mol_per_cubic_meter, convert_to_mass_mixing_ratio
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_kinds,                only: kind_phys
    use iso_c_binding,             only: c_double
    real(kind_phys),                   intent(inout) :: time_step            ! s
    real(kind_phys), target,           intent(inout) :: temperature(:,:)     ! K
    real(kind_phys), target,           intent(inout) :: pressure(:,:)        ! Pa
    real(kind_phys), target,           intent(inout) :: dry_air_density(:,:) ! kg m-3
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
    real(kind_phys), target,           intent(inout) :: constituents(:,:,:)  ! kg kg-1
    real(kind_phys), target,           intent(inout) :: rate_params(:,:,:)
    real(kind_phys), target,           intent(in)    :: height(:,:)          ! km
    character(len=512),                intent(out)   :: errmsg
    integer,                           intent(out)   :: errcode

    ! local variables
    real(c_double), target, dimension(size(temperature, dim=1)     &
                                    * size(temperature, dim=2))      :: m_temperature
    real(c_double), target, dimension(size(pressure, dim=1)        &
                                    * size(pressure, dim=2))         :: m_pressure
    real(c_double), target, dimension(size(dry_air_density, dim=1) &
                                    * size(dry_air_density, dim=2))  :: m_dry_air_density
    real(c_double), target, dimension(size(constituents, dim=1)    &
                                    * size(constituents, dim=2)    & 
                                    * size(constituents, dim=3))     :: m_constituents ! mol m-3
    real(c_double), target, dimension(size(rate_params, dim=1)     &
                                    * size(rate_params, dim=2)     & 
                                    * size(rate_params, dim=3))      :: m_rate_params
    real(kind_phys), target, dimension(size(constituents, dim=3))    :: molar_mass_arr ! kg mol-1
    integer :: i_elem

    call tuvx_run(height, temperature, dry_air_density, errmsg, errcode)

    ! Get the molar_mass that is set in the call to instantiate()
    do i_elem = 1, size(molar_mass_arr)
      call constituent_props(i_elem)%molar_mass(molar_mass_arr(i_elem), errcode, errmsg)
      if (errcode /= 0) then
        errmsg = "[MUSICA Error] Unable to get molar mass."
        return
      end if
    end do

    ! TODO(jiwon) Check molar mass is non zero as it becomes a denominator for unit converison
    ! this code needs to go when ccpp framework does the check
    do i_elem = 1, size(molar_mass_arr)
      if (molar_mass_arr(i_elem) <= 0) then
        errcode = 1
        errmsg = "[MUSICA Error] Molar mass must be greater than zero."
        return
      end if
    end do

    ! Convert CAM-SIMA unit to MICM unit (kg kg-1  ->  mol m-3)
    call convert_to_mol_per_cubic_meter(dry_air_density, molar_mass_arr, constituents)

    ! Reshape array (3D -> 1D) and convert type (kind_phys -> c_double)
    call reshape_into_micm_arr(temperature, pressure, dry_air_density, constituents, rate_params, &
                      m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)

    call micm_run(time_step, m_temperature, m_pressure, m_dry_air_density, m_constituents, &
                      m_rate_params, errmsg, errcode)

    ! Reshape array (1D -> 3D) and convert type (c_double -> kind_phys)
    call reshape_into_ccpp_arr(temperature, pressure, dry_air_density, constituents, rate_params, &
                      m_temperature, m_pressure, m_dry_air_density, m_constituents, m_rate_params)

    ! Convert MICM unit back to CAM-SIMA unit (mol m-3  ->  kg kg-1)
    call convert_to_mass_mixing_ratio(dry_air_density, molar_mass_arr, constituents)

  end subroutine musica_ccpp_run

  !> \section arg_table_musica_ccpp_final Argument Table
  !! \htmlinclude musica_ccpp_final.html
  subroutine musica_ccpp_final(errmsg, errcode)
    integer,            intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    call tuvx_final(errmsg, errcode)
    call micm_final(errmsg, errcode)

  end subroutine musica_ccpp_final

end module musica_ccpp
