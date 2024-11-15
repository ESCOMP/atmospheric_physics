! MPAS and SE dynamical core specific
! 1) scaling of temperature for enforcing energy consistency
! 2) and to ensure correct computation of temperature dependent diagnostic tendencies, e.g., dtcore
module dycore_energy_consistency_adjust
  use ccpp_kinds, only: kind_phys
  implicit none
  private

  public :: dycore_energy_consistency_adjust_run

contains
!> \section arg_table_dycore_energy_consistency_adjust_run Argument Table
!! \htmlinclude arg_table_dycore_energy_consistency_adjust_run.html
  subroutine dycore_energy_consistency_adjust_run( &
      ncol, pver, &
      do_consistency_adjust, &
      scaling_dycore, &
      tend_dTdt, &
      T, &
      tend_dTdt_local)

    ! Input arguments
    integer,            intent(in)    :: ncol                  ! number of atmospheric columns
    integer,            intent(in)    :: pver                  ! number of vertical layers
    logical,            intent(in)    :: do_consistency_adjust ! do energy consistency adjustment?
    real(kind_phys),    intent(in)    :: scaling_dycore(:,:)   ! scaling for conversion of temperature increment [1]
    real(kind_phys),    intent(in)    :: tend_dTdt(:,:)        ! model physics temperature tendency [K s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: T(:,:)                ! air temperature [K]

    ! Output arguments
    real(kind_phys),    intent(out) :: tend_dTdt_local(:,:)  ! (scheme) temperature tendency [K s-1]

    if (do_consistency_adjust) then
      ! original formula for scaling of temperature:
      !   T(:ncol,:) = temp_ini(:ncol,:) + &
      !                scaling_dycore(:ncol,:) * (T(:ncol,:) - temp_ini(:ncol,:))
      ! and temperature tendency due to model physics:
      !   tend_dTdt(:ncol,:) = scaling_dycore(:ncol,:) * tend_dTdt(:ncol,:)
      !
      ! the terms can be arranged for this scaling to be applied through scheme tendencies
      ! at the cost of a round-off level difference
      tend_dTdt_local(:ncol,:) = (scaling_dycore(:ncol,:) - 1._kind_phys) * tend_dTdt(:ncol,:)
    endif
    ! do nothing for dynamical cores with energy consistent with CAM physics

  end subroutine dycore_energy_consistency_adjust_run

end module dycore_energy_consistency_adjust
