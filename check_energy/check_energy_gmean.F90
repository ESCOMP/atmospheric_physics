module check_energy_gmean
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public  :: check_energy_gmean

contains

  ! Compute global mean total energy of physics input and output states
  ! computed consistently with dynamical core vertical coordinate
  ! (under hydrostatic assumption)
!> \section arg_table_check_energy_gmean_run Argument Table
!! \htmlinclude arg_table_check_energy_gmean_run.html
  subroutine check_energy_gmean_run( &
       ncol, pver, &
       pint, &
       te_ini_dyn, teout, &
       heat_glob)

    ! To-be-CCPP-ized
    use gmean_mod, only: gmean

    ! Input arguments
    integer,            intent(in)    :: ncol           ! number of atmospheric columns
    integer,            intent(in)    :: pver           ! number of vertical layers
    real(kind_phys),    intent(in)    :: pint(:,:)      ! interface pressure [Pa]
    real(kind_phys),    intent(in)    :: te_ini_dyn(:)  ! dycore  formula: initial total energy [J m-2]
    real(kind_phys),    intent(in)    :: teout(:)       ! total energy for global fixer in next timestep [J m-2]
    real(kind_phys),    intent(out)   :: heat_glob      ! global mean heating rate [J kg-1 s-1]

    ! Local variables
    real(kind_phys) :: te(ncol, 4)                      ! total energy of input/output states (copy)
    real(kind_phys) :: te_glob(4)                       ! global means of total energy

    ! These used to be module variables in check_energy and could be made into intent(out) if needed.
    real(kind_phys) :: teout_glob                       ! global mean energy of output state
    real(kind_phys) :: teinp_glob                       ! global mean energy of input state
    real(kind_phys) :: tedif_glob                       ! global mean energy difference
    real(kind_phys) :: psurf_glob                       ! global mean surface pressure
    real(kind_phys) :: ptopb_glob                       ! global mean top boundary pressure

    ! nb hplin: in CAM te(i, lchnk, 2) gets teout from pbuf and this is chunkized.
    ! so check_energy_gmean_run ccppized will not be used in CAM because of this chunk handling,
    ! and so will gmean be different between CAM and CAM-SIMA.

    ! Copy total energy out of input and output states.

    ! Compute global means of input and output energies and of surface pressure for heating rate.
    ! (assuming uniform ptop)

    ! Compute global mean total energy difference for check_energy_fix
  end subroutine check_energy_gmean_run

end module check_energy_gmean
