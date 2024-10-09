module musica_ccpp_tuvx
  use iso_c_binding

  ! Note: "tuvx_t" is included in an external pre-built tuvx library that the host
  ! model is responsible for linking to during compilation
  use musica_tuvx, only: tuvx_t, grid_t
  use musica_ccpp_util, only: has_error_occurred
  use ccpp_kinds, only: kind_phys
  use musica_ccpp_namelist, only: filename_of_tuvx_configuration

  implicit none
  private

  public :: tuvx_init, tuvx_run, tuvx_final

  type(tuvx_t), pointer  :: tuvx => null( )
  type(grid_t), pointer  :: height_grid => null( )

contains

  !> Intitialize TUVX
  subroutine tuvx_init(vertical_layer_dimension, &
      vertical_interface_dimension, errmsg, errcode)
    use musica_tuvx, only: grid_map_t, grid_t, profile_map_t, radiator_map_t
    use musica_util, only: error_t, mapping_t
    use musica_ccpp_tuvx_height_grid, only: create_height_grid, &
                                            height_grid_label, height_grid_units

    integer,            intent(in)  :: vertical_layer_dimension     ! (count)
    integer,            intent(in)  :: vertical_interface_dimension ! (count)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! local variables
    type(grid_map_t),     pointer :: grids
    type(profile_map_t),  pointer :: profiles
    type(radiator_map_t), pointer :: radiators
    type(error_t)                 :: error

    errcode = 0
    errmsg = ''

    grids => grid_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) return

    height_grid => create_height_grid( vertical_layer_dimension, &
        vertical_interface_dimension, errmsg, errcode )
    if (errcode /= 0) return
    call grids%add( height_grid, error )
    if (has_error_occurred( error, errmsg, errcode )) return

    profiles => profile_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( grids )
      return
    end if

    radiators => radiator_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( grids )
      deallocate( profiles )
      return
    end if

    tuvx => tuvx_t( filename_of_tuvx_configuration, grids, profiles, &
                    radiators, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( grids )
      deallocate( profiles )
      deallocate( radiators )
      return
    end if

    deallocate( height_grid )
    deallocate( grids )
    deallocate( profiles )
    deallocate( radiators )

    grids => tuvx%get_grids( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      tuvx => null()
      return
    end if

    height_grid => grids%get( height_grid_label, height_grid_units, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      deallocate( tuvx )
      tuvx => null()
      deallocate( grids )
      return
    end if

    deallocate( grids )

  end subroutine tuvx_init

  !> Calculates photolysis rate constants for the current model conditions
  subroutine tuvx_run( temperature, dry_air_density, height_midpoints, &
      height_interfaces, photolysis_rate_constants, errmsg, errcode )
    use musica_util, only: error_t
    use musica_ccpp_tuvx_height_grid, only: set_height_grid_values

    real(kind_phys),    intent(in)  :: temperature(:,:)       ! K (column, layer)
    real(kind_phys),    intent(in)  :: dry_air_density(:,:)   ! molecule cm-3 (column, layer)
    real(kind_phys),    intent(in)  :: height_midpoints(:,:)  ! km (column, layer)
    real(kind_phys),    intent(in)  :: height_interfaces(:,:) ! km (column, interface)
    ! temporarily set to Chapman mechanism and 1 dimension
    ! until mapping between MICM and TUV-x is implemented
    real(kind_phys),    intent(out) :: photolysis_rate_constants(:) ! s-1 (column, reaction)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! local variables
    type(error_t) :: error
    integer :: i_col

    errcode = 0
    errmsg = ''

    do i_col = 1, size(temperature, dim=1)
      call set_height_grid_values( height_grid, height_midpoints(i_col,:), &
          height_interfaces(i_col,:), errmsg, errcode )
    end do
    if (errcode /= 0) return

    ! stand-in until actual photolysis rate constants are calculated
    photolysis_rate_constants(:) = 1.0e-6_kind_phys

  end subroutine tuvx_run

  !> Finalize tuvx
  subroutine tuvx_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errcode = 0
    errmsg = ''
    deallocate( height_grid )

  end subroutine tuvx_final

end module musica_ccpp_tuvx