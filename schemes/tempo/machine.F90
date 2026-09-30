!> The TEMPO microphysics scheme depends on `machine` for real kind parameters
!> instead of `ccpp_kinds`.
module machine
    use ccpp_kinds, only: kind_phys
    use, intrinsic :: iso_fortran_env, only: kind_sngl_prec => real32, kind_dbl_prec => real64

    implicit none

    private
    public :: kind_phys
    public :: kind_sngl_prec, kind_dbl_prec
end module machine
