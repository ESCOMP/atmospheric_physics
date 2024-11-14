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
      temp_ini, &
      t, &
      tend_dtdt)

    ! Input arguments
    integer,            intent(in)    :: ncol                  ! number of atmospheric columns
    integer,            intent(in)    :: pver                  ! number of vertical layers
    logical,            intent(in)    :: do_consistency_adjust ! do energy consistency adjustment?
    real(kind_phys),    intent(in)    :: scaling_dycore(:,:)   ! scaling for conversion of temperature increment [1]
    real(kind_phys),    intent(in)    :: temp_ini(:,:)         ! initial temperature [K]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: T(:,:)                ! temperature [K]
    real(kind_phys),    intent(inout) :: tend_dtdt(:,:)        ! model phys temperature tendency [K s-1]

    if(do_consistency_adjust) then
      T(:ncol,:) = temp_ini(:ncol,:) + &
                   scaling_dycore(:ncol,:) * (T(:ncol,:) - temp_ini(:ncol,:))

      tend_dtdt(:ncol,:) = scaling_dycore(:ncol,:) * tend_dtdt(:ncol,:)
    endif
    ! do nothing for dynamical cores with energy consistent with CAM physics

  end subroutine dycore_energy_consistency_adjust_run

end module dycore_energy_consistency_adjust
