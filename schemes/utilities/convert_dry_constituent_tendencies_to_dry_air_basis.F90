! Convert dry constituent tendencies from moist to dry air basis
module convert_dry_constituent_tendencies_to_dry_air_basis
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: convert_dry_constituent_tendencies_to_dry_air_basis_run

contains

!> \section arg_table_convert_dry_constituent_tendencies_to_dry_air_basis_run Argument Table
!! \htmlinclude convert_dry_constituent_tendencies_to_dry_air_basis_run.html
  subroutine convert_dry_constituent_tendencies_to_dry_air_basis_run( &
    ncol, pver, pcnst, &
    pdel, pdeldry, &
    const_props, &
    tend_q, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! Input arguments
    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: pcnst
    real(kind_phys),  intent(in)    :: pdel(:, :)          ! Layer thickness (moist air) [Pa]
    real(kind_phys),  intent(in)    :: pdeldry(:, :)       ! Layer thickness (dry air) [Pa]
    
    ! Framework dependency for constituent properties
    type(ccpp_constituent_prop_ptr_t), &
                      intent(in)    :: const_props(:)      ! CCPP constituent properties pointer

    ! Input/Output arguments
    real(kind_phys),  intent(inout) :: tend_q(:, :, :)     ! Constituent tendencies [kg kg-1 s-1]

    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k, m
    logical :: const_is_dry

    errmsg = ''
    errflg = 0

    ! Convert the tendencies for the dry constituents to dry air basis
    do m = 1, pcnst
       ! Check if this constituent is dry type
       call const_props(m)%is_dry(const_is_dry, errflg, errmsg)
       if (errflg /= 0) return
       
       if (const_is_dry) then
          do k = 1, pver
             do i = 1, ncol
                tend_q(i, k, m) = tend_q(i, k, m) * pdel(i, k) / pdeldry(i, k)
             end do
          end do
       end if
    end do

  end subroutine convert_dry_constituent_tendencies_to_dry_air_basis_run

end module convert_dry_constituent_tendencies_to_dry_air_basis
