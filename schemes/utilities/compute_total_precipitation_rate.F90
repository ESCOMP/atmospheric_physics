! Compute the total surface precipitation rate (deep + shallow convective + stratiform)
module compute_total_precipitation_rate
  implicit none
  private

  public :: compute_total_precipitation_rate_run

contains

  !> \section arg_table_compute_total_precipitation_rate_run  Argument Table
  !! \htmlinclude compute_total_precipitation_rate_run.html
  subroutine compute_total_precipitation_rate_run(ncol, prec_dp, prec_sh, prec_str, prect, &
       errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,            intent(in)  :: ncol         ! number of columns
    real(kind_phys),    intent(in)  :: prec_dp(:)   ! LWE deep convective precipitation rate at surface [m s-1]
    real(kind_phys),    intent(in)  :: prec_sh(:)   ! LWE shallow convective precipitation rate at surface [m s-1]
    real(kind_phys),    intent(in)  :: prec_str(:)  ! LWE stratiform precipitation rate at surface [m s-1]
    real(kind_phys),    intent(out) :: prect(:)     ! LWE total precipitation rate at surface [m s-1]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    prect(:ncol) = prec_dp(:ncol) + prec_sh(:ncol) + prec_str(:ncol)

  end subroutine compute_total_precipitation_rate_run

end module compute_total_precipitation_rate
