
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
!  subroutine rrtmgp_cloud_optics_setup_init(liq_filename, abs_lw_liq_out, &
!                  ext_sw_liq_out, ssa_sw_liq_out, asm_sw_liq_out, g_lambda_out, g_mu_out, errmsg, errflg)
  subroutine rrtmgp_cloud_optics_setup_init(liq_filename, ice_filename, abs_lw_liq_out, abs_lw_ice_out, &
                  ext_sw_liq_out, ext_sw_ice_out, ssa_sw_liq_out, ssa_sw_ice_out, asm_sw_liq_out,       &
                  asm_sw_ice_out, g_lambda_out, g_mu_out, g_d_eff_out, errmsg, errflg)
    use ccpp_kinds,     only: kind_phys
    use ccpp_io_reader, only: abstract_netcdf_reader_t, create_netcdf_reader_t
    ! Inputs
    character(len=*),                   intent(in) :: liq_filename     ! Full file path for liquid optics file
    character(len=*),                   intent(in) :: ice_filename     ! Full file path for ice optics file
    ! Outputs
    real(kind_phys), dimension(:,:,:), allocatable, intent(out) :: abs_lw_liq_out    ! Longwave mass specific absorption for in-cloud liquid water path
    real(kind_phys), dimension(:,:,:), allocatable, intent(out) :: ext_sw_liq_out
    real(kind_phys), dimension(:,:,:), allocatable, intent(out) :: ssa_sw_liq_out
    real(kind_phys), dimension(:,:,:), allocatable, intent(out) :: asm_sw_liq_out
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: abs_lw_ice_out    ! Longwave mass specific absorption for in-cloud ice water path
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: ext_sw_ice_out
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: ssa_sw_ice_out
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: asm_sw_ice_out
    real(kind_phys), dimension(:,:),   allocatable, intent(out) :: g_lambda_out      ! lambda scale samples on grid
    real(kind_phys), dimension(:),     allocatable, intent(out) :: g_mu_out          ! Mu samples on grid
    real(kind_phys), dimension(:),     allocatable, intent(out) :: g_d_eff_out       ! Radiative effective diameter samples on grid

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    class(abstract_netcdf_reader_t), allocatable :: pio_reader
    real(kind_phys), dimension(:),     pointer :: g_mu
    real(kind_phys), dimension(:),     pointer :: g_d_eff
    real(kind_phys), dimension(:,:),   pointer :: g_lambda
    real(kind_phys), dimension(:,:),   pointer :: ext_sw_ice
    real(kind_phys), dimension(:,:),   pointer :: ssa_sw_ice
    real(kind_phys), dimension(:,:),   pointer :: asm_sw_ice
    real(kind_phys), dimension(:,:),   pointer :: abs_lw_ice
    real(kind_phys), dimension(:,:,:), pointer :: ext_sw_liq
    real(kind_phys), dimension(:,:,:), pointer :: ssa_sw_liq
    real(kind_phys), dimension(:,:,:), pointer :: asm_sw_liq
    real(kind_phys), dimension(:,:,:), pointer :: abs_lw_liq
    character(len=256) :: alloc_errmsg
    character(len=*), parameter :: sub = 'rrtmgp_cloud_optics_setup_init'

    ! Set error variables
    errmsg = ''
    errflg = 0

    pio_reader = create_netcdf_reader_t()

    ! Open liquid optics file
    call pio_reader%open_file(liq_filename, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Read in variables
    call pio_reader%get_var('mu', g_mu, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('lambda', g_lambda, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('k_ext_sw', ext_sw_liq, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('ssa_sw', ssa_sw_liq, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('asm_sw', asm_sw_liq, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('k_abs_lw', abs_lw_liq, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Close the liquid optics file
    call pio_reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Convert kext from m^2/Volume to m^2/Kg
    ext_sw_liq = ext_sw_liq / 0.9970449e3_kind_phys
    abs_lw_liq = abs_lw_liq / 0.9970449e3_kind_phys

    ! Open the ice optics file
    call pio_reader%open_file(ice_filename, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Read in variables
    call pio_reader%get_var('d_eff', g_d_eff, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('sw_ext', ext_sw_ice, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('sw_ssa', ssa_sw_ice, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('sw_asm', asm_sw_ice, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if
    call pio_reader%get_var('lw_abs', abs_lw_ice, errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Close the ice optics file
    call pio_reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
       return
    end if

    ! Allocate output variables
    allocate(g_mu_out(size(g_mu)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating g_mu_out, message: ', alloc_errmsg
       return
    end if
    allocate(g_lambda_out(size(g_lambda,1), size(g_lambda,2)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating g_lambda_out, message: ', alloc_errmsg
       return
    end if
    allocate(g_d_eff_out(size(g_d_eff)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating g_d_eff_out, message: ', alloc_errmsg
       return
    end if
    allocate(ext_sw_liq_out(size(ext_sw_liq,1),size(ext_sw_liq,2),size(ext_sw_liq,3)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating ext_sw_liq_out, message: ', alloc_errmsg
       return
    end if
    allocate(ext_sw_ice_out(size(ext_sw_ice,1),size(ext_sw_ice,2)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating ext_sw_ice_out, message: ', alloc_errmsg
       return
    end if
    allocate(asm_sw_liq_out(size(asm_sw_liq,1),size(asm_sw_liq,2),size(asm_sw_liq,3)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating asm_sw_liq_out, message: ', alloc_errmsg
       return
    end if
    allocate(asm_sw_ice_out(size(asm_sw_ice,1),size(asm_sw_ice,2)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating asm_sw_ice_out, message: ', alloc_errmsg
       return
    end if
    allocate(ssa_sw_liq_out(size(ssa_sw_liq,1),size(ssa_sw_liq,2),size(ssa_sw_liq,3)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating ssa_sw_liq_out, message: ', alloc_errmsg
       return
    end if
    allocate(ssa_sw_ice_out(size(ssa_sw_ice,1),size(ssa_sw_ice,2)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating ssa_sw_ice_out, message: ', alloc_errmsg
       return
    end if
    allocate(abs_lw_liq_out(size(abs_lw_liq,1),size(abs_lw_liq,2),size(abs_lw_liq,3)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating abs_lw_liq_out, message: ', alloc_errmsg
       return
    end if
    allocate(abs_lw_ice_out(size(abs_lw_ice,1),size(abs_lw_ice,2)), stat=errflg, errmsg=alloc_errmsg)
    if (errflg /= 0) then
       write(errmsg, '(a,a,a)') sub, ': ERROR allocating abs_lw_ice_out, message: ', alloc_errmsg
       return
    end if

    ext_sw_liq_out = ext_sw_liq
    ext_sw_ice_out = ext_sw_ice
    ssa_sw_liq_out = ssa_sw_liq
    ssa_sw_ice_out = ssa_sw_ice
    asm_sw_liq_out = asm_sw_liq
    asm_sw_ice_out = asm_sw_ice
    abs_lw_liq_out = abs_lw_liq
    abs_lw_ice_out = abs_lw_ice
    g_mu_out       = g_mu
    g_lambda_out   = g_lambda
    g_d_eff_out    = g_d_eff

  end subroutine rrtmgp_cloud_optics_setup_init

!==============================================================================

end module rrtmgp_cloud_optics_setup
