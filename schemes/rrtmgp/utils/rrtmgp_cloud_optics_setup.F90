
!> This module contains two routines: The first initializes data and functions
!! needed to compute the longwave cloud radiative properties in RRTMGP. The second routine
!! is a ccpp scheme within the "radiation loop", where the shortwave optical properties
!! (optical-depth, single-scattering albedo, asymmetry parameter) are computed for ALL
!! cloud types visible to RRTMGP.
module rrtmgp_cloud_optics_setup
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  public :: rrtmgp_cloud_optics_setup_init

contains

  ! ######################################################################################
  ! SUBROUTINE rrtmgp_cloud_optics_setup_init()
  ! ######################################################################################
!> \section arg_table_rrtmgp_cloud_optics_setup_init Argument Table
!! \htmlinclude rrtmgp_cloud_optics_setup_init.html
!!
  subroutine rrtmgp_cloud_optics_setup_init(liq_filename, ice_filename, abs_lw_liq, abs_lw_ice,     &
                  ext_sw_liq, ext_sw_ice, ssa_sw_liq, ssa_sw_ice, asm_sw_liq, asm_sw_ice, g_lambda, &
                  g_mu, g_d_eff, errmsg, errflg)
    use ccpp_kinds,     only: kind_phys
    use ccpp_io_reader, only: abstract_netcdf_reader_t, create_netcdf_reader_t
    ! Inputs
    character(len=*),                   intent(in) :: liq_filename     ! Full file path for liquid optics file
    character(len=*),                   intent(in) :: ice_filename     ! Full file path for ice optics file
    ! Outputs
    real(kind_phys), dimension(:,:,:), allocatable, intent(out) :: abs_lw_liq    ! Longwave mass specific absorption for in-cloud liquid water path
    real(kind_phys), dimension(:,:,:), allocatable, intent(out) :: ext_sw_liq
    real(kind_phys), dimension(:,:,:), allocatable, intent(out) :: ssa_sw_liq
    real(kind_phys), dimension(:,:,:), allocatable, intent(out) :: asm_sw_liq
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: abs_lw_ice    ! Longwave mass specific absorption for in-cloud ice water path
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: ext_sw_ice
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: ssa_sw_ice
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: asm_sw_ice
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: g_lambda      ! lambda scale samples on grid
    real(kind_phys), dimension(:),     allocatable, intent(out) :: g_mu          ! Mu samples on grid
    real(kind_phys), dimension(:),     allocatable, intent(out) :: g_d_eff       ! Radiative effective diameter samples on grid

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    real(kind_phys), parameter :: liquid_water_density = 0.9970449e3_kind_phys
    class(abstract_netcdf_reader_t), allocatable :: file_reader
    character(len=256) :: alloc_errmsg
    character(len=*), parameter :: sub = 'rrtmgp_cloud_optics_setup_init'

    ! Set error variables
    errmsg = ''
    errflg = 0

    file_reader = create_netcdf_reader_t()

    ! Open liquid optics file
    call file_reader%open_file(liq_filename, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Read in variables
    call file_reader%get_var('mu', g_mu, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('lambda', g_lambda, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('k_ext_sw', ext_sw_liq, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('ssa_sw', ssa_sw_liq, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('asm_sw', asm_sw_liq, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('k_abs_lw', abs_lw_liq, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Close the liquid optics file
    call file_reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Convert kext from m^2/Volume to m^2/Kg
    ext_sw_liq = ext_sw_liq / liquid_water_density
    abs_lw_liq = abs_lw_liq / liquid_water_density

    ! Open the ice optics file
    call file_reader%open_file(ice_filename, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Read in variables
    call file_reader%get_var('d_eff', g_d_eff, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('sw_ext', ext_sw_ice, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('sw_ssa', ssa_sw_ice, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('sw_asm', asm_sw_ice, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call file_reader%get_var('lw_abs', abs_lw_ice, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Close the ice optics file
    call file_reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

  end subroutine rrtmgp_cloud_optics_setup_init

!==============================================================================

end module rrtmgp_cloud_optics_setup
