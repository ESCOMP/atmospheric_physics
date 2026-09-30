
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

  integer,         public              :: nmu, nlambda
  real(kind_phys), public, allocatable :: g_mu(:)
  real(kind_phys), public, allocatable :: g_lambda(:,:)
  real(kind_phys), public, allocatable :: abs_lw_liq(:,:,:)
  real(kind_phys), public, allocatable :: ext_sw_liq(:,:,:)
  real(kind_phys), public, allocatable :: asm_sw_liq(:,:,:)
  real(kind_phys), public, allocatable :: ssa_sw_liq(:,:,:)
  integer,         public              :: n_g_d
  real(kind_phys), public, allocatable :: g_d_eff(:)
  real(kind_phys), public, allocatable :: abs_lw_ice(:,:)
  real(kind_phys), public, allocatable :: ext_sw_ice(:,:)
  real(kind_phys), public, allocatable :: asm_sw_ice(:,:)
  real(kind_phys), public, allocatable :: ssa_sw_ice(:,:)

  ! Cloud optics parameterization selectors (from namelist), branched on by
  ! rrtmgp_sw_cloud_optics and rrtmgp_lw_cloud_optics.
  character(len=32), public :: liq_cld_optics = 'unset' ! gammadist or slingo
  character(len=32), public :: ice_cld_optics = 'unset' ! mitchell or ebertcurry

contains

  ! ######################################################################################
  ! SUBROUTINE rrtmgp_cloud_optics_setup_init()
  ! ######################################################################################
!> \section arg_table_rrtmgp_cloud_optics_setup_init Argument Table
!! \htmlinclude rrtmgp_cloud_optics_setup_init.html
!!
  subroutine rrtmgp_cloud_optics_setup_init(liq_filename, ice_filename, liq_cld_optics_nl, &
                  ice_cld_optics_nl, errmsg, errflg)
    use ccpp_io_reader, only: abstract_netcdf_reader_t, create_netcdf_reader_t
    ! Inputs
    character(len=*),                   intent(in) :: liq_filename      ! Full file path for liquid optics file
    character(len=*),                   intent(in) :: ice_filename      ! Full file path for ice optics file
    character(len=*),                   intent(in) :: liq_cld_optics_nl ! Liquid cloud optics type (namelist)
    character(len=*),                   intent(in) :: ice_cld_optics_nl ! Ice cloud optics type (namelist)
    ! Outputs
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    ! Local variables
    real(kind_phys), parameter :: liquid_water_density = 0.9970449e3_kind_phys
    class(abstract_netcdf_reader_t), pointer :: file_reader
    character(len=256) :: alloc_errmsg
    character(len=*), parameter :: sub = 'rrtmgp_cloud_optics_setup_init'

    ! Set error variables
    errmsg = ''
    errflg = 0

    ! Validate and store the cloud optics parameterization selectors.
    liq_cld_optics = liq_cld_optics_nl
    ice_cld_optics = ice_cld_optics_nl
    if (trim(liq_cld_optics) /= 'gammadist' .and. trim(liq_cld_optics) /= 'slingo') then
       write(errmsg,'(a,a,a,a)') sub, ': liq_cld_optics must be either slingo or gammadist, got "', trim(liq_cld_optics), '"'
       errflg = 1
       return
    end if
    if (trim(ice_cld_optics) /= 'mitchell' .and. trim(ice_cld_optics) /= 'ebertcurry') then
       write(errmsg,'(a,a,a,a)') sub, ': ice_cld_optics must be either ebertcurry or mitchell, got "', trim(ice_cld_optics), '"'
       errflg = 1
       return
    end if

    file_reader => create_netcdf_reader_t()

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
    deallocate(file_reader)
    nullify(file_reader)

    ! Set size module variables
    nmu = size(g_mu)
    nlambda = size(g_lambda, 2)
    n_g_d = size(g_d_eff)

  end subroutine rrtmgp_cloud_optics_setup_init

!==============================================================================

end module rrtmgp_cloud_optics_setup
