module musica_ccpp_tuvx
  use iso_c_binding

  ! Note: "tuvx_t" is included in an external pre-built tuvx library that the host
  ! model is responsible for linking to during compilation
  use musica_tuvx, only: tuvx_t
  use musica_ccpp_util, only: has_error_occurred
  use ccpp_kinds, only: kind_phys
  use musica_ccpp_namelist, only: filename_of_tuvx_configuration

  implicit none
  private

  public :: tuvx_init, tuvx_run, tuvx_final

  type(tuvx_t), pointer  :: tuvx => null( )

contains

  !> Intitialize TUVX
  subroutine tuvx_init(n_vertical_levels, errmsg, errcode)
    use musica_tuvx, only: grid_map_t, grid_t, profile_map_t, profile_t
    use musica_util, only: error_t, mapping_t

    integer,            intent(in)  :: n_vertical_levels  ! Number of vertical levels per column
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    type(error_t) :: error

    errcode = 0
    errmsg = ''

    tuvx => tuvx_t( filename_of_tuvx_configuration, error )
    if (has_error_occurred( error, errmsg, errcode )) return

    grids => grid_map_t( error )
    height_grid => grid_t( "height", "km", n_vertical_levels )
    call grids%add( height, error )

    profiles => profile_map_t( error )
    temperature_profile => profile_t( "temperature", "K", height, error )
    call profiles%add( temperature_profile, error )

    ! ! unpack the core for each OMP thread on every MPI process
    ! allocate( tuvx_wrappers )
    ! do i_thread = 1, omp_threads
    ! associate( wrapper => tuvx_wrappers( i_thread ) )
    !   allocate( wrapper%core_ )
    !   pos = 0
    !   wrapper%height_      => wrapper%core_%get_updater( height      )
    !   wrapper%temperature_ => wrapper%core_%get_updater( temperature )
    !   wrapper%air_         => wrapper%core_%get_updater( air         )
    ! end associate
    ! end do

    deallocate( grids )
    deallocate( profiles )

  end subroutine tuvx_init

  subroutine tuvx_run( height, air_density, temperature, photolysis_rate_constants &
      photolysis_rate_constants, errmsg, errcode )
    ! Calculates photolysis rate constants for the current model conditions
    use musica_util, only: error_t

    real(kind=dk),      intent(in)  :: height(:,:) ! height above sea level [km] (layer, column)
    real(kind=dk),      intent(in)  :: air_density(:,:) ! number density of dry (?) air [molecule cm-3] (layer, column)
    real(kind=dk),      intent(in)  :: temperature(:,:) ! temperature [K] (layer, column)
    real(kind=dk),      intent(out) :: photolysis_rate_constants(:,:,:) ! photolysis rate constants [s-1] (reaction, layer, column)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode
    integer :: i_col

    errcode = 0
    errmsg = ''

    ! do i_col = 1, size( height, 2 )
    !   call wrapper%height_%update( height( :, i_col ) )
    !   call wrapper%air_%update( air_density( :, i_col ) )
    !   call wrapper%temperature_%update( temperature( :, i_col ) )
    !   call core%calculate( photolysis_rate_constants                        &
    !                       = photolysis_rate_constants( :, :, i_col ) )
    ! end do

  end subroutine tuvx_run

  !> Finalize tuvx
  subroutine tuvx_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errcode = 0
    errmsg = ''

  end subroutine tuvx_final

end module musica_ccpp_tuvx
