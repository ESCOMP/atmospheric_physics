! Clamp number concentration constituent tendencies
! so that after tendency application, values remain within [qmin, qmax].
!
! This scheme must be placed immediately after apply_constituent_tendencies
! in the SDF to replicate the same behavior as in physics_types.F90 in CAM.
module clamp_number_concentrations

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  public :: clamp_number_concentrations_init
  public :: clamp_number_concentrations_run

  ! Number of number-concentration species to clamp
  integer, parameter :: num_species = 4

  ! Constituent indices (-1 = not present; looked up in _init)
  integer :: ix_numliq  = -1
  integer :: ix_numrai  = -1
  integer :: ix_numice  = -1
  integer :: ix_numsno  = -1

  ! Clamp bounds
  real(kind_phys), parameter :: qmin = 1.0e-12_kind_phys ! minimum number concentration [kg-1]
  real(kind_phys), parameter :: qmax = 1.0e10_kind_phys  ! maximum number concentration [kg-1]

  ! Standard names for each species (order matches ix_ variables above)
  character(len=75), parameter :: std_names(num_species) = &
    (/'mass_number_concentration_of_cloud_liquid_wrt_moist_air_and_condensed_water', &
      'mass_number_concentration_of_rain_wrt_moist_air_and_condensed_water        ', &
      'mass_number_concentration_of_ice_wrt_moist_air_and_condensed_water         ', &
      'mass_number_concentration_of_snow_wrt_moist_air_and_condensed_water        '/)

contains

!> \section arg_table_clamp_number_concentrations_init Argument Table
!! \htmlinclude clamp_number_concentrations_init.html
  subroutine clamp_number_concentrations_init(const_props, errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,          only: ccpp_const_get_idx

    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)  :: const_props(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: ix_species(num_species)
    integer :: n
    character(len=512) :: local_errmsg
    integer            :: local_errflg

    errmsg = ''
    errflg = 0

    ix_species(:) = -1
    ! Look up constituent indices. A missing species is skipped
    do n = 1, num_species
      call ccpp_const_get_idx(const_props, trim(std_names(n)), &
                              ix_species(n), local_errmsg, local_errflg)
      if (local_errflg /= 0) then
        ! Constituent not found — mark as inactive, reset error
        ix_species(n) = -1
        local_errflg  = 0
        local_errmsg  = ''
      end if
    end do

    ! Store in module variables
    ix_numliq = ix_species(1)
    ix_numrai = ix_species(2)
    ix_numice = ix_species(3)
    ix_numsno = ix_species(4)

  end subroutine clamp_number_concentrations_init

!> \section arg_table_clamp_number_concentrations_run Argument Table
!! \htmlinclude clamp_number_concentrations_run.html
  subroutine clamp_number_concentrations_run( &
    ncol, pver, dt, &
    const_q, const_tend, &
    errmsg, errflg)

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: dt                  ! physics timestep [s]
    real(kind_phys),    intent(inout) :: const_q(:, :, :)    ! constituent mixing ratios
    real(kind_phys),    intent(inout) :: const_tend(:, :, :) ! constituent tendencies
    character(len=512), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables
    integer :: ix_species(num_species)
    integer :: n, i, k, ix
    real(kind_phys) :: q_projected  ! projected value after tendency application [kg-1]

    errmsg = ''
    errflg = 0

    ix_species = (/ix_numliq, ix_numrai, ix_numice, ix_numsno/)

    do n = 1, num_species
      ix = ix_species(n)
      if (ix <= 0) cycle

      do k = 1, pver
        do i = 1, ncol
          if (const_q(i, k, ix) < qmin) then
            const_q(i, k, ix)    = qmin
            const_tend(i, k, ix) = 0.0_kind_phys
          else if (const_q(i, k, ix) > qmax) then
            const_q(i, k, ix)    = qmax
            const_tend(i, k, ix) = 0.0_kind_phys
          end if
        end do
      end do
    end do

    ! The below alternative is a formulation based on the principle that we do not alter
    ! the constituent tendencies themselves.
    !
    ! However, the (qmin - q) / dt will cause floating point errors and introduce
    ! more harm than good (this is unfortunately because qmin and q may often have
    ! a large difference)
    !
    ! e.g.,
    ! (gdb) p const_q(356,23,1)
    ! $15 = 1826.1689801273931            <-- q >> qmin.
    ! (gdb) p const_tend(356,23,1)
    ! $16 = -1.0145383222929962
    ! (gdb) p dt
    ! $17 = 1800
    ! (gdb) p const_tend(356,23,1)*dt
    ! $18 = -1826.1689801273933
    ! (gdb) p const_q(356,23,1)+const_tend(356,23,1)*dt
    ! $19 = -2.2737367544323206e-13      <-- projected fall below qmin
    !                                        which triggers the clamp
    ! yet because q >> qmin, q + (qmin - q)/dt*dt will be far from qmin
    ! introducing numerical noise. (hplin, 3/18/26)
    !
    ! do n = 1, num_species
    !   ix = ix_species(n)
    !   if (ix <= 0) cycle

    !   do k = 1, pver
    !     do i = 1, ncol
    !       ! Project what the value would be after tendency application
    !       q_projected = const_q(i, k, ix) + const_tend(i, k, ix) * dt

    !       if (q_projected < qmin) then
    !         ! Rewrite tendency so the updater lands exactly on qmin
    !         const_tend(i, k, ix) = (qmin - const_q(i, k, ix)) / dt
    !       else if (q_projected > qmax) then
    !         ! Rewrite tendency so the updater lands exactly on qmax
    !         const_tend(i, k, ix) = (qmax - const_q(i, k, ix)) / dt
    !       end if
    !     end do
    !   end do
    ! end do

  end subroutine clamp_number_concentrations_run

end module clamp_number_concentrations
