! Seasalt emissions for modal aerosol model
! Seasalt section number fluxes accumulated into modal number and mass surface fluxes
! over the ocean fraction.
module modal_seasalt_emissions
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_seasalt_emissions_run

  ! Sea salt aerosol material density used for the number-to-mass flux
  ! conversion (value from CAM mo_constants).
  real(kind_phys), parameter :: seasalt_density = 2.2e+3_kind_phys  ! [kg m-3]

contains

  subroutine modal_seasalt_emissions_run(ncol, nslt, seasalt_indices, emis_scale, &
                                         u_bottom, v_bottom, zmid_bottom, &
                                         srf_temp, ocnfrc, pi, cflx)

    use sslt_sections, only: nsections, fluxes, Dg, rdry

    integer, intent(in) :: ncol
    integer, intent(in) :: nslt                        ! number of seasalt bins (mass species)
    integer, intent(in) :: seasalt_indices(:)          ! constituent indices: mass bins, then number bins
    real(kind_phys), intent(in) :: emis_scale          ! sea salt emission tuning factor
    real(kind_phys), intent(in) :: u_bottom(:)         ! bottom layer zonal wind (m/s)
    real(kind_phys), intent(in) :: v_bottom(:)         ! bottom layer meridional wind (m/s)
    real(kind_phys), intent(in) :: zmid_bottom(:)      ! bottom layer midpoint geopotential height above surface (m)
    real(kind_phys), intent(in) :: srf_temp(:)         ! sea surface temperature (K)
    real(kind_phys), intent(in) :: ocnfrc(:)           ! ocean fraction
    real(kind_phys), intent(in) :: pi
    real(kind_phys), intent(inout) :: cflx(:, :)        ! constituent surface fluxes (kg/m2/s or #/m2/s)

    ! local vars
    integer  :: mn, mm, ibin, i
    real(kind_phys) :: fi(ncol, nsections)
    real(kind_phys) :: u10cubed(ncol)
    real(kind_phys), parameter :: z0 = 0.0001_kind_phys  ! m roughness length over oceans--from ocean model

    real(kind_phys) :: sst_sz_range_lo(nslt)
    real(kind_phys) :: sst_sz_range_hi(nslt)

    u10cubed(:ncol) = sqrt(u_bottom(:ncol)**2 + v_bottom(:ncol)**2)

    ! move the winds to 10m high from the midpoint of the gridbox:
    ! follows Tie and Seinfeld and Pandis, p.859 with math.

    u10cubed(:ncol) = u10cubed(:ncol)*log(10._kind_phys/z0)/log(zmid_bottom(:ncol)/z0)

    ! we need them to the 3.41 power, according to Gong et al., 1997:
    u10cubed(:ncol) = u10cubed(:ncol)**3.41_kind_phys

    if (nslt == 4) then
      sst_sz_range_lo(:) = [0.08e-6_kind_phys, 0.02e-6_kind_phys, 0.3e-6_kind_phys, 1.0e-6_kind_phys] ! accu, aitken, fine, coarse
      sst_sz_range_hi(:) = [0.3e-6_kind_phys, 0.08e-6_kind_phys, 1.0e-6_kind_phys, 10.0e-6_kind_phys]
    else if (nslt == 3) then
      sst_sz_range_lo(:) = [0.08e-6_kind_phys, 0.02e-6_kind_phys, 1.0e-6_kind_phys]  ! accu, aitken, coarse
      sst_sz_range_hi(:) = [1.0e-6_kind_phys, 0.08e-6_kind_phys, 10.0e-6_kind_phys]
    end if

    fi(:ncol, :nsections) = fluxes(srf_temp, u10cubed, ncol)

    do ibin = 1, nslt
      mm = seasalt_indices(ibin)
      mn = seasalt_indices(nslt + ibin)

      if (mn > 0) then
        do i = 1, nsections
          if (Dg(i) >= sst_sz_range_lo(ibin) .and. Dg(i) < sst_sz_range_hi(ibin)) then
            cflx(:ncol, mn) = cflx(:ncol, mn) + fi(:ncol, i)*ocnfrc(:ncol)*emis_scale  !++ ag: scale sea-salt
          end if
        end do
      end if

      cflx(:ncol, mm) = 0.0_kind_phys
      do i = 1, nsections
        if (Dg(i) >= sst_sz_range_lo(ibin) .and. Dg(i) < sst_sz_range_hi(ibin)) then
          cflx(:ncol, mm) = cflx(:ncol, mm) + fi(:ncol, i)*ocnfrc(:ncol)*emis_scale &   !++ ag: scale sea-salt
                            *4._kind_phys/3._kind_phys*pi*rdry(i)**3*seasalt_density  ! should use dry size, convert from number to mass flux (kg/m2/s)
        end if
      end do

    end do

  end subroutine modal_seasalt_emissions_run

end module modal_seasalt_emissions
