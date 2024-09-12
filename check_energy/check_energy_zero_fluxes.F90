! zeros input fluxes to check_energy
! before running other schemes
module check_energy_zero_fluxes

    use ccpp_kinds, only: kind_phys

    implicit none
    private

    public :: check_energy_zero_fluxes_run

contains

!> \section arg_table_check_energy_zero_fluxes_run Argument Table
!! \htmlinclude arg_table_check_energy_zero_fluxes_run.html
    subroutine check_energy_zero_fluxes_run(name, flx_vap, flx_cnd, flx_ice, flx_sen)
        character(len=*),   intent(out)    :: name           ! parameterization name for fluxes
        real(kind_phys),    intent(out)    :: flx_vap(:)     ! boundary flux of vapor [kg m-2 s-1]
        real(kind_phys),    intent(out)    :: flx_cnd(:)     ! boundary flux of liquid+ice (precip?) [m s-1]
        real(kind_phys),    intent(out)    :: flx_ice(:)     ! boundary flux of ice (snow?) [m s-1]
        real(kind_phys),    intent(out)    :: flx_sen(:)     ! boundary flux of sensible heat [W m-2]

        ! reset values to zero
        name = ''
        flx_vap(:) = 0._kind_phys
        flx_cnd(:) = 0._kind_phys
        flx_ice(:) = 0._kind_phys
        flx_sen(:) = 0._kind_phys
    end subroutine check_energy_zero_fluxes_run

end module check_energy_zero_fluxes