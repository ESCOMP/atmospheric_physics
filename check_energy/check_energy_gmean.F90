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
       ncol, pver, dtime, &
       gravit, &
       pint, &
       te_ini_dyn, teout, &
       tedif_glob, heat_glob)

    ! Dependency: Uses gmean from src/utils
    use gmean_mod, only: gmean

    ! Input arguments
    integer,            intent(in)    :: ncol           ! number of atmospheric columns
    integer,            intent(in)    :: pver           ! number of vertical layers
    real(kind_phys),    intent(in)    :: dtime          ! physics time step [s]
    real(kind_phys),    intent(in)    :: gravit         ! gravitational acceleration [m s-2]
    real(kind_phys),    intent(in)    :: pint(:,:)      ! interface pressure [Pa]
    real(kind_phys),    intent(in)    :: te_ini_dyn(:)  ! dycore  formula: initial total energy [J m-2]
    real(kind_phys),    intent(in)    :: teout(:)       ! total energy for global fixer in next timestep [J m-2]

    ! Output arguments
    real(kind_phys),    intent(out)   :: tedif_glob     ! global mean energy difference [J m-2]
    real(kind_phys),    intent(out)   :: heat_glob      ! global mean heating rate [J kg-1 s-1]

    ! Local variables
    real(kind_phys) :: te(ncol, 4)                      ! total energy of input/output states (copy)
    real(kind_phys) :: te_glob(4)                       ! global means of total energy

    real(kind_phys) :: teinp_glob                       ! global mean energy of input state [J m-2]
    real(kind_phys) :: teout_glob                       ! global mean energy of output state [J m-2]
    real(kind_phys) :: psurf_glob                       ! global mean surface pressure [Pa]
    real(kind_phys) :: ptopb_glob                       ! global mean top boundary pressure [Pa]

    ! DEVNOTE hplin: in CAM te(i, lchnk, 2) gets teout from pbuf and this is chunkized.
    ! so check_energy_gmean_run ccppized will not be used in CAM because of this chunk handling,
    ! and so will gmean be different between CAM and CAM-SIMA.

    ! Copy total energy out of input and output states.
    ! These four fields will have their global means calculated respectively
    te(:ncol, 1) = te_ini_dyn(:ncol)        ! Input energy using dycore energy formula [J m-2]
    te(:ncol, 2) = teout(:ncol)             ! Total energy from end of physics timestep [J m-2]
    te(:ncol, 3) = pint(:ncol, pver+1)      ! Surface pressure for heating rate [Pa]
    te(:ncol, 4) = pint(:ncol, 1)           ! Model top pressure for heating rate [Pa]
                                            ! not constant for z-based vertical coordinate

    ! Compute global means of input and output energies and of surface pressure for heating rate.
    ! (assuming uniform ptop)
    call gmean(te, te_glob, 4)
    teinp_glob = te_glob(1)
    teout_glob = te_glob(2)
    psurf_glob = te_glob(3)
    ptopb_glob = te_glob(4)

    ! Compute global mean total energy difference for check_energy_fix
    tedif_glob = teinp_glob - teout_glob
    heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)    ! [J kg-1 s-1]
  end subroutine check_energy_gmean_run

end module check_energy_gmean
