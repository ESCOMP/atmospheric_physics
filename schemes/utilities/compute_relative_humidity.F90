! Compute relative humidity for the aerosol optics table lookup.
! In CAM this is computed inline in aerosol_optics_cam (radiation_tend);
! in CAM-SIMA it is a standalone scheme for suites where no upstream
! scheme provides relative humidity (compute_cloud_fraction does in
! cam4/cam5, but is not part of cam7).
module compute_relative_humidity

  implicit none
  private

  public :: compute_relative_humidity_run

contains

  !> \section arg_table_compute_relative_humidity_run  Argument Table
  !! \htmlinclude compute_relative_humidity_run.html
  subroutine compute_relative_humidity_run(ncol, pver, t, pmid, qv, relh, &
       errmsg, errflg)
    use ccpp_kinds,    only: kind_phys
    use wv_saturation, only: qsat

    integer,            intent(in)  :: ncol       ! number of columns
    integer,            intent(in)  :: pver       ! number of vertical layers
    real(kind_phys),    intent(in)  :: t(:,:)     ! air temperature [K]
    real(kind_phys),    intent(in)  :: pmid(:,:)  ! air pressure [Pa]
    real(kind_phys),    intent(in)  :: qv(:,:)    ! water vapor mixing ratio wrt moist air [kg kg-1]
    real(kind_phys),    intent(out) :: relh(:,:)  ! relative humidity [fraction]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(kind_phys) :: sate(ncol, pver)  ! saturation vapor pressure [Pa]
    real(kind_phys) :: satq(ncol, pver)  ! saturation specific humidity [kg kg-1]

    errmsg = ''
    errflg = 0

    ! calculate relative humidity for table lookup into rh grid
    call qsat(t, pmid, sate, satq, ncol, pver)
    relh = qv / satq
    relh = max(1.e-20_kind_phys, relh)

  end subroutine compute_relative_humidity_run

end module compute_relative_humidity
