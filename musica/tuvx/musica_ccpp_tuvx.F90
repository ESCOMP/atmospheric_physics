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
  subroutine tuvx_init(errmsg, errcode)
    use musica_tuvx, only: grid_map_t, profile_map_t, radiator_map_t
    use musica_util, only: error_t, mapping_t

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
    if (has_error_occurred( error, errmsg, errcode )) then
      write(*,*) "[MUSICA Error] ", errmsg
      return
    end if

    profiles => profile_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      write(*,*) "[MUSICA Error] ", errmsg
      return
    end if

    radiators =>radiator_map_t( error )
    if (has_error_occurred( error, errmsg, errcode )) then
      write(*,*) "[MUSICA Error] ", errmsg
      return
    end if

    ! TODO(jiwon) - MUSICA TUVX constuctor needs update
    ! tuvx => tuvx_t( filename_of_tuvx_configuration, grids, profiles, radiators, error )
    tuvx => tuvx_t( filename_of_tuvx_configuration, error )
    if (has_error_occurred( error, errmsg, errcode )) then
      write(*,*) "[MUSICA Error] ", errmsg
      return
    end if

    deallocate( grids )
    deallocate( profiles )
    deallocate( radiators )
    deallocate( tuvx )

  end subroutine tuvx_init

  !> Calculates photolysis rate constants for the current model conditions
  subroutine tuvx_run( height, temperature, dry_air_density, errmsg, errcode )
    use musica_util, only: error_t

    real(kind_phys),    intent(in)  :: height(:,:)          ! km (layer, column)
    real(kind_phys),    intent(in)  :: temperature(:,:)     ! K (layer, column)
    real(kind_phys),    intent(in)  :: dry_air_density(:,:) ! molecule cm-3 (layer, column)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! local variables
    type(error_t) :: error

    errcode = 0
    errmsg = ''

  end subroutine tuvx_run

  !> Finalize tuvx
  subroutine tuvx_final(errmsg, errcode)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    errcode = 0
    errmsg = ''

  end subroutine tuvx_final

end module musica_ccpp_tuvx