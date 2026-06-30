!> The mere existence of this module is to satisfy the misdirected dependency of the MMM physics,
!> which inexplicably depends on `ccpp_kind_types` instead of `ccpp_kinds`.
module ccpp_kind_types
    use ccpp_kinds, only: kind_phys
    use iso_fortran_env, only: kind_phys8 => real64

    implicit none

    private
    public :: kind_phys
    public :: kind_phys8
end module ccpp_kind_types
