!> Top-level wrapper for MUSICA chemistry components
module musica_ccpp
  use musica_ccpp_micm, only : micm_register, micm_init, micm_run, micm_final

  implicit none
  private

  public :: musica_register, musica_init, musica_run, musica_final

contains

  subroutine musica_register(constituents, errcode, errmsg)
    use ccpp_constituent_prop_mod, only : ccpp_constituent_properties_t
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituents(:)
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    call micm_register(constituents, errcode, errmsg)
  end subroutine musica_register

  !> \section arg_table_musica_init Argument Table
  !! \htmlinclude musica_init.html
  subroutine musica_init(errcode, errmsg)
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    call micm_init(errcode, errmsg)
  end subroutine musica_init

  !> \section arg_table_musica_run Argument Table
  !! \htmlinclude musica_run.html
  subroutine musica_run(time_step, temperature, pressure, dry_air_density, constituent_props, &
                      constituents, iulog, errcode, errmsg)
    use ccpp_kinds, only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    real(kind_phys),                   intent(in)    :: time_step            ! s
    real(kind_phys),                   intent(in)    :: temperature(:,:)     ! K
    real(kind_phys),                   intent(in)    :: pressure(:,:)        ! Pa
    real(kind_phys),                   intent(in)    :: dry_air_density(:,:) ! kg m-3
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: constituent_props(:)
    real(kind_phys),                   intent(inout) :: constituents(:,:,:)  ! kg kg-1
    integer,                           intent(in)    :: iulog
    integer,                           intent(out)   :: errcode
    character(len=512),                intent(out)   :: errmsg

    call micm_run(time_step, temperature, pressure, dry_air_density, constituent_props, &
                  constituents, iulog, errcode, errmsg)

  end subroutine musica_run

  !> \section arg_table_musica_final Argument Table
  !! \htmlinclude musica_final.html
  subroutine musica_final(errcode, errmsg)
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    call micm_final(errcode, errmsg)
  end subroutine musica_final

end module musica_ccpp