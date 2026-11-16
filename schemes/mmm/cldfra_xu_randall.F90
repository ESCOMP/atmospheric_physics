module cldfra_xu_randall
    implicit none

    private
    public :: cldfra_xu_randall_run
contains
    !> This subroutine implements the Xu-Randall cloud fraction scheme from Xu and Randall (1996).
    !> See doi:10.1175/1520-0469(1996)053<3084:ASCPFU>2.0.CO;2 for details.
    !>
    !> This CCPP-conformant version is based on the implementations in MPAS (`calc_cldfraction`) and WRF (`cal_cldfra1`).
    !>
    !> \section arg_table_cldfra_xu_randall_run Argument Table
    !! \htmlinclude cldfra_xu_randall_run.html
    pure subroutine cldfra_xu_randall_run( &
            ncol, pver, epsilon, svpt0, &
            pres, temp, qv, qc, qi, qs, &
            cldfra, relhum, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys

        integer, intent(in) :: ncol, pver
        real(kind_phys), intent(in) :: epsilon, svpt0, &
                                       pres(:, :), temp(:, :), qv(:, :), qc(:, :), qi(:, :), qs(:, :)
        real(kind_phys), intent(out) :: cldfra(:, :), relhum(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        ! Constants from Xu and Randall (1996).
        real(kind_phys), parameter :: alpha0 = 100.0_kind_phys ! Correspond to the $\alpha_0$ constant.
        real(kind_phys), parameter :: gamma  = 0.49_kind_phys  ! Correspond to the $\gamma$ constant.
        real(kind_phys), parameter :: p      = 0.25_kind_phys  ! Correspond to the $p$ constant.

        ! Constants from Equation 6 in Murray (1967).
        ! See doi:10.1175/1520-0450(1967)006<0203:OTCOSV>2.0.CO;2 for details.
        real(kind_phys), parameter :: svp1  = 0.61078_kind_phys
        real(kind_phys), parameter :: svp2  = 17.2693882_kind_phys
        real(kind_phys), parameter :: svp3  = 35.86_kind_phys
        real(kind_phys), parameter :: svpi2 = 21.8745584_kind_phys
        real(kind_phys), parameter :: svpi3 = 7.66_kind_phys

        ! Implementation-defined thresholds.
        real(kind_phys), parameter :: qcldmin    = 1.0E-12_kind_phys
        real(kind_phys), parameter :: relhumcrit = 1.0_kind_phys ! Critical RH to be considered grid-scale saturation.

        integer :: i, k
        real(kind_phys) :: esw, esi, qvsw, qvsi, qcld, weight, qvs
        real(kind_phys) :: exponent, satdef

        cldfra(:, :) = 0.0_kind_phys
        relhum(:, :) = 0.0_kind_phys

        do k = 1, pver
            do i = 1, ncol
                ! Saturation vapor pressures wrt water and ice.
                esw = 1000.0_kind_phys * svp1 * exp(svp2  * (temp(i, k) - svpt0) / (temp(i, k) - svp3 ))
                esi = 1000.0_kind_phys * svp1 * exp(svpi2 * (temp(i, k) - svpt0) / (temp(i, k) - svpi3))

                ! Saturation vapor mixing ratios wrt water and ice.
                qvsw = epsilon * esw / (pres(i, k) - esw)
                qvsi = epsilon * esi / (pres(i, k) - esi)

                ! "Cloudy" mixing ratio, which corresponds to the $\bar{q_l}$ term in Xu and Randall (1996).
                qcld = qc(i, k) + qi(i, k) + qs(i, k)

                ! "Icy" weight of the cloudy mixing ratio.
                if (qcld < qcldmin) then
                    weight = 0.0_kind_phys
                else
                    weight = (qi(i, k) + qs(i, k)) / qcld
                end if

                ! Weighted saturation vapor mixing ratio and RH.
                qvs = (1.0_kind_phys - weight) * qvsw + weight * qvsi
                relhum(i, k) = qv(i, k) / qvs

                if (qcld < qcldmin) then
                    ! Assume that the cloud fraction is 0 when the cloudy mixing ratio is below the minimum threshold.
                    cldfra(i, k) = 0.0_kind_phys
                else if (relhum(i, k) >= relhumcrit) then
                    ! Assume that the cloud fraction is 1 when the RH exceeds the critical RH.
                    ! This is the RH >= 1 branch of Equation 4 in Xu and Randall (1996).
                    cldfra(i, k) = 1.0_kind_phys
                else
                    ! Saturation deficiency scaled by the critical RH.
                    ! 1.0E-10 is an implementation-defined limit.
                    satdef = max(1.0E-10_kind_phys, relhumcrit * qvs - qv(i, k))
                    relhum(i, k) = max(1.0E-10_kind_phys, relhum(i, k))

                    ! This is the RH < 1 branch of Equation 4 in Xu and Randall (1996).
                    ! 1.0E-3 is an implementation-defined limit.
                    exponent = -alpha0 * qcld / (satdef**gamma)
                    cldfra(i, k) = ((relhum(i, k) / relhumcrit)**p) * (1.0_kind_phys - max(1.0E-3_kind_phys, exp(exponent)))

                    if (cldfra(i, k) < 0.01_kind_phys) then
                        cldfra(i, k) = 0.0_kind_phys
                    end if
                end if
            end do
        end do

        errmsg = ''
        errflg = 0
    end subroutine cldfra_xu_randall_run
end module cldfra_xu_randall
