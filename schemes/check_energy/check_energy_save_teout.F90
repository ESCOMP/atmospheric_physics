! save total energy for global fixer in next timestep
! this must be called after the last parameterization and physics_update,
! and after a final check_energy_chng to compute te_cur.
module check_energy_save_teout
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: check_energy_save_teout_run

contains

!> \section arg_table_check_energy_save_teout_run Argument Table
!! \htmlinclude arg_table_check_energy_save_teout_run.html
  subroutine check_energy_save_teout_run(ncol, te_cur_dyn, teout)

    ! Input arguments
    integer,            intent(in)    :: ncol           ! number of atmospheric columns
    real(kind_phys),    intent(in)    :: te_cur_dyn (:) ! dycore  formula: current total energy [J m-2]

    ! Output arguments
    real(kind_phys),    intent(out)   :: teout(:)       ! total energy for global fixer in next timestep [J m-2]

    ! nb hplin: note that in physpkg.F90, the pbuf is updated to the previous dyn timestep
    ! through itim_old. Need to check if we need to replicate such pbuf functionality
    ! in the CAM-SIMA/CCPP infrastructure.
    teout(:ncol) = te_cur_dyn(:ncol)

  end subroutine check_energy_save_teout_run

end module check_energy_save_teout