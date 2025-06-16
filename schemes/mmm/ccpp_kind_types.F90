!> The mere existence of this module is to satisfy the misdirected dependency of MMM physics,
!> which inexplicably depends on `ccpp_kind_types` instead of `ccpp_kinds`.
module ccpp_kind_types
    use ccpp_kinds, only: kind_phys

    implicit none

    private
    public :: kind_phys
contains
end module ccpp_kind_types
