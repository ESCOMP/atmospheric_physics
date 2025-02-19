! This is a non-portable wrapper subroutine for cam_thermo_water_update
! in the cam_thermo module.
module thermo_water_update
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: thermo_water_update_run

contains

  ! Update water dependent properties
!> \section arg_table_thermo_water_update_run Argument Table
!! \htmlinclude arg_table_thermo_water_update_run.html
  subroutine thermo_water_update_run( &
    mmr,                   &
    ncol, pver,            &
    energy_formula_dycore, &
    pdel, pdeldry,         &
    cp_or_cv_dycore)

    ! This scheme is non-portable due to dependencies on cam_thermo
    ! for the actual logic of cam_thermo_water_update, which depends on air_composition
    ! and a series of other subroutines/module properties
    use cam_thermo, only: cam_thermo_water_update

    ! Input arguments
    real(kind_phys), intent(in)  :: mmr(:,:,:)            ! constituent mass mixing ratios [kg kg-1]
    integer,         intent(in)  :: ncol                  ! number of atmospheric columns
    integer,         intent(in)  :: pver                  ! number of vertical layers
    integer,         intent(in)  :: energy_formula_dycore ! total energy formulation used by dycore
    real(kind_phys), intent(in)  :: pdel(:,:)             ! layer thickness [Pa]
    real(kind_phys), intent(in)  :: pdeldry(:,:)          ! dry layer thickness [Pa]

    ! Output arguments
    real(kind_phys), intent(out) :: cp_or_cv_dycore(:,:)  ! enthalpy or heat capacity, dycore dependent [J K-1 kg-1]

    call cam_thermo_water_update( &
      mmr             = mmr,                            &  ! mmr*factor is a dry mixing ratio
      ncol            = ncol,                           &
      pver            = pver,                           &
      energy_formula  = energy_formula_dycore,          &
      cp_or_cv_dycore = cp_or_cv_dycore(:ncol,:),       &
      to_dry_factor   = pdel(:ncol,:)/pdeldry(:ncol,:)  &  ! factor to convert to dry
    )

  end subroutine thermo_water_update_run

end module thermo_water_update
