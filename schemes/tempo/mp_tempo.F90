!> This module contains the interstitial schemes for the TEMPO microphysics scheme.
module mp_tempo
    implicit none

    private
    public :: mp_tempo_init
    public :: mp_tempo_run
contains
    !> \section arg_table_mp_tempo_init Argument Table
    !! \htmlinclude mp_tempo_init.html
    subroutine mp_tempo_init( &
            num_subcyc, &
            aerosol_aware, hail_aware, refl_incl_melt, sat_adj, semi_lagr_sedi, &
            con_pi, con_hvap, con_hfus, con_rv, con_rd, con_cp, con_t0c, con_rgas, con_rhow, &
            tempo_cfgs, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys
        use module_mp_tempo_cfgs, only: ty_tempo_cfgs
        use module_mp_tempo_driver, only: tempo_init
        use module_mp_tempo_params, only: pi, lvap0, lfus, lsub, rv, rdry, cp, t0, r_uni, rho_w, &
                                          initialize_parameters

        integer, intent(in) :: num_subcyc
        logical, intent(in) :: aerosol_aware, hail_aware, refl_incl_melt, sat_adj, semi_lagr_sedi
        real(kind_phys), intent(in) :: con_pi, con_hvap, con_hfus, con_rv, con_rd, con_cp, con_t0c, con_rgas, con_rhow
        type(ty_tempo_cfgs), intent(out) :: tempo_cfgs
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        errmsg = ''
        errflg = 0

        ! Check for logic errors in namelist options.

        if (num_subcyc < 1) then
            errmsg = '"num_subcyc" must be positive for TEMPO microphysics.'
            errflg = 1

            return
        end if

        if (sat_adj) then
            if (aerosol_aware .or. hail_aware) then
                errmsg = '"sat_adj=.true." should only be run with "aerosol_aware=.false." and "hail_aware=.false."' // &
                    ' for TEMPO microphysics.'
                errflg = 1

                return
            end if
        end if

        call tempo_init( &
            aerosolaware_flag=aerosol_aware, &
            hailaware_flag=hail_aware, &
            refl10cm_from_melting_flag=refl_incl_melt, &
            cloud_condensation_flag=(.not. sat_adj), &
            semi_sedi_flag=semi_lagr_sedi, &
            tempo_cfgs=tempo_cfgs)

        pi = con_pi
        lvap0 = con_hvap
        lfus = con_hfus
        lsub = lvap0 + lfus
        rv = con_rv
        rdry = con_rd
        cp = con_cp
        t0 = con_t0c
        ! Convert from J "kmol-1" K-1 to J "mol-1" K-1.
        r_uni = con_rgas / 1000.0_kind_phys
        rho_w = con_rhow

        call initialize_parameters()
    end subroutine mp_tempo_init

    !> \section arg_table_mp_tempo_run Argument Table
    !! \htmlinclude mp_tempo_run.html
    subroutine mp_tempo_run( &
            ncol, pver, curr_step, num_subcyc, &
            initial_run, &
            dtime_phys, &
            dz, w, p, exner, th, &
            qv, qc, qr, qi, qs, qg, nr, ni, &
            nc, nwfa, nifa, nwfa2d, nifa2d, &
            ng, volg, &
            tempo_cfgs, &
            rain, ice, snow, graupel, prcp, &
            tend_th, tend_qv, tend_qc, tend_qr, tend_qi, tend_qs, tend_qg, tend_nr, tend_ni, &
            frozen_fraction, refl10cm, re_cloud, re_ice, re_snow, &
            tend_nc, tend_nwfa, tend_nifa, &
            tend_ng, tend_volg, &
            errmsg, errflg)
        use ccpp_kinds, only: kind_phys
        use module_mp_tempo_aerosols, only: init_water_friendly_aerosols, init_ice_friendly_aerosols
        use module_mp_tempo_cfgs, only: ty_tempo_cfgs
        use module_mp_tempo_driver, only: ty_tempo_driver_diags, &
                                          tempo_aerosol_surface_emissions, tempo_run
        use module_mp_tempo_params, only: eps, nain1

        integer, intent(in) :: ncol, pver, curr_step, num_subcyc
        logical, intent(in) :: initial_run
        real(kind_phys), intent(in) :: dtime_phys
        real(kind_phys), intent(in) :: dz(:, :), w(:, :), p(:, :), exner(:, :), th(:, :)
        real(kind_phys), intent(in) :: qv(:, :), qc(:, :), qr(:, :), qi(:, :), qs(:, :), qg(:, :), nr(:, :), ni(:, :)
        ! Only used when aerosol-aware.
        real(kind_phys), intent(in) :: nc(:, :)
        real(kind_phys), intent(inout) :: nwfa(:, :), nifa(:, :), nwfa2d(:), nifa2d(:)
        ! Only used when hail-aware.
        real(kind_phys), intent(in) :: ng(:, :), volg(:, :)
        type(ty_tempo_cfgs), intent(in) :: tempo_cfgs

        real(kind_phys), intent(out) :: rain(:), ice(:), snow(:), graupel(:), prcp(:)
        real(kind_phys), intent(out) :: tend_th(:, :)
        real(kind_phys), intent(out) :: tend_qv(:, :), tend_qc(:, :), tend_qr(:, :)
        real(kind_phys), intent(out) :: tend_qi(:, :), tend_qs(:, :), tend_qg(:, :)
        real(kind_phys), intent(out) :: tend_nr(:, :), tend_ni(:, :)
        real(kind_phys), intent(out) :: frozen_fraction(:), refl10cm(:, :), re_cloud(:, :), re_ice(:, :), re_snow(:, :)
        ! Only used when aerosol-aware.
        real(kind_phys), intent(out) :: tend_nc(:, :), tend_nwfa(:, :), tend_nifa(:, :)
        ! Only used when hail-aware.
        real(kind_phys), intent(out) :: tend_ng(:, :), tend_volg(:, :)
        character(*), intent(out) :: errmsg
        integer, intent(out) :: errflg

        integer :: i
        integer :: ids, ide, jds, jde, kds, kde, &
                   ims, ime, jms, jme, kms, kme, &
                   its, ite, jts, jte, kts, kte
        real(kind_phys) :: dtime_subcyc
        ! TEMPO expects vertical indexes to be in ascending order from bottom to top of atmosphere,
        ! which is the exact opposite to CAM-SIMA. Variables with the "_r" suffix have a reversed vertical dimension.
        ! Variables with the "_1" suffix have an extra dimension of size 1. In TEMPO, 2d arrays are dimensioned (i, j) while
        ! 3d arrays are dimensioned (i, k, j). Even if the j dimension is not used, it still needs to be 1.
        real(kind_phys) :: dz_r1(ncol, pver, 1), w_r1(ncol, pver, 1), p_r1(ncol, pver, 1), exner_r1(ncol, pver, 1)
        real(kind_phys) :: new_th_r1(ncol, pver, 1)
        real(kind_phys) :: new_qv_r1(ncol, pver, 1), new_qc_r1(ncol, pver, 1), new_qr_r1(ncol, pver, 1)
        real(kind_phys) :: new_qi_r1(ncol, pver, 1), new_qs_r1(ncol, pver, 1), new_qg_r1(ncol, pver, 1)
        real(kind_phys) :: new_nr_r1(ncol, pver, 1), new_ni_r1(ncol, pver, 1)
        ! Only used when aerosol-aware.
        real(kind_phys), allocatable :: nwfa_r1(:, :, :)
        real(kind_phys), allocatable :: nwfa2d_1(:, :)
        real(kind_phys), allocatable :: new_nc_r1(:, :, :)
        real(kind_phys), allocatable :: new_nwfa_r1(:, :, :)
        real(kind_phys), allocatable :: new_nifa_r1(:, :, :)
        ! Only used when hail-aware.
        real(kind_phys), allocatable :: new_ng_r1(:, :, :)
        real(kind_phys), allocatable :: new_volg_r1(:, :, :)
        type(ty_tempo_driver_diags) :: tempo_driver_diags

        errmsg = ''
        errflg = 0

        rain(:) = 0.0_kind_phys
        ice(:) = 0.0_kind_phys
        snow(:) = 0.0_kind_phys
        graupel(:) = 0.0_kind_phys
        prcp(:) = 0.0_kind_phys

        tend_th(:, :) = 0.0_kind_phys
        tend_qv(:, :) = 0.0_kind_phys
        tend_qc(:, :) = 0.0_kind_phys
        tend_qr(:, :) = 0.0_kind_phys
        tend_qi(:, :) = 0.0_kind_phys
        tend_qs(:, :) = 0.0_kind_phys
        tend_qg(:, :) = 0.0_kind_phys
        tend_nr(:, :) = 0.0_kind_phys
        tend_ni(:, :) = 0.0_kind_phys

        tend_nc(:, :) = 0.0_kind_phys
        tend_nwfa(:, :) = 0.0_kind_phys
        tend_nifa(:, :) = 0.0_kind_phys

        tend_ng(:, :) = 0.0_kind_phys
        tend_volg(:, :) = 0.0_kind_phys

        dz_r1(:, :, 1) = dz(:, pver:1:-1)
        w_r1(:, :, 1) = w(:, pver:1:-1)
        p_r1(:, :, 1) = p(:, pver:1:-1)
        exner_r1(:, :, 1) = exner(:, pver:1:-1)

        ! Check for existing CCN and IN aerosol data. If missing,
        ! fill in just a basic vertical profile that is somewhat boundary layer following.
        if (initial_run .and. tempo_cfgs % aerosolaware_flag) then
            ! Potential cloud condensation nuclei (CCN).
            if (maxval(nwfa) < eps) then
                do i = 1, ncol
                    call init_water_friendly_aerosols(dz_r1(i, :, 1), nwfa(i, :))

                    nwfa2d(i) = nwfa(i, 1) * 0.000196_kind_phys * (50.0_kind_phys / dz_r1(i, 1, 1))
                end do

                ! Flip upside down to be consistent with CAM-SIMA. Index 1 is at the top of atmosphere.
                nwfa(:, :) = nwfa(:, pver:1:-1)
            else
                if (maxval(nwfa2d) < eps) then
                    do i = 1, ncol
                        nwfa2d(i) = nwfa(i, pver) * 0.000196_kind_phys * (5.0_kind_phys / dz_r1(i, 1, 1))
                    end do
                end if
            end if

            ! Potential ice nuclei (IN).
            if (maxval(nifa) < eps) then
                do i = 1, ncol
                    call init_ice_friendly_aerosols(dz_r1(i, :, 1), nifa(i, :))

                    nifa2d(i) = 0.0_kind_phys
                end do

                ! Flip upside down to be consistent with CAM-SIMA. Index 1 is at the top of atmosphere.
                nifa(:, :) = nifa(:, pver:1:-1)
            else
                if (maxval(nifa2d) < eps) then
                    nifa2d(:) = 0.0_kind_phys
                end if
            end if

            ! Ensure non-negative aerosol number concentrations.
            where (nwfa <= 0.0_kind_phys)
                nwfa = 1.1E6_kind_phys
            end where

            where (nifa <= 0.0_kind_phys)
                nifa = nain1 * 0.01_kind_phys
            end where
        end if

        new_th_r1(:, :, 1) = th(:, pver:1:-1)
        new_qv_r1(:, :, 1) = qv(:, pver:1:-1)
        new_qc_r1(:, :, 1) = qc(:, pver:1:-1)
        new_qr_r1(:, :, 1) = qr(:, pver:1:-1)
        new_qi_r1(:, :, 1) = qi(:, pver:1:-1)
        new_qs_r1(:, :, 1) = qs(:, pver:1:-1)
        new_qg_r1(:, :, 1) = qg(:, pver:1:-1)
        new_nr_r1(:, :, 1) = nr(:, pver:1:-1)
        new_ni_r1(:, :, 1) = ni(:, pver:1:-1)

        if (tempo_cfgs % aerosolaware_flag) then
            allocate(nwfa_r1(ncol, pver, 1), nwfa2d_1(ncol, 1), &
                errmsg=errmsg, stat=errflg)

            if (errflg /= 0) then
                errmsg = 'mp_tempo_run: Failed to allocate "nwfa_r1", "nwfa2d_1"' // new_line('') // &
                    trim(adjustl(errmsg))

                return
            end if

            allocate(new_nc_r1(ncol, pver, 1), new_nwfa_r1(ncol, pver, 1), new_nifa_r1(ncol, pver, 1), &
                errmsg=errmsg, stat=errflg)

            if (errflg /= 0) then
                errmsg = 'mp_tempo_run: Failed to allocate "new_nc_r1", "new_nwfa_r1", "new_nifa_r1"' // new_line('') // &
                    trim(adjustl(errmsg))

                return
            end if

            new_nc_r1(:, :, 1) = nc(:, pver:1:-1)
            new_nwfa_r1(:, :, 1) = nwfa(:, pver:1:-1)
            new_nifa_r1(:, :, 1) = nifa(:, pver:1:-1)
        end if

        if (tempo_cfgs % hailaware_flag) then
            allocate(new_ng_r1(ncol, pver, 1), new_volg_r1(ncol, pver, 1), &
                errmsg=errmsg, stat=errflg)

            if (errflg /= 0) then
                errmsg = 'mp_tempo_run: Failed to allocate "new_ng_r1", "new_volg_r1"' // new_line('') // &
                    trim(adjustl(errmsg))

                return
            end if

            new_ng_r1(:, :, 1) = ng(:, pver:1:-1)
            new_volg_r1(:, :, 1) = volg(:, pver:1:-1)
        end if

        ! Set subcycle timestep.
        dtime_subcyc = dtime_phys / real(num_subcyc, kind_phys)

        ! Set internal dimensions.
        ids = 1
        ims = 1
        its = 1
        ide = ncol
        ime = ncol
        ite = ncol
        jds = 1
        jms = 1
        jts = 1
        jde = 1
        jme = 1
        jte = 1
        kds = 1
        kms = 1
        kts = 1
        kde = pver
        kme = pver
        kte = pver

        do i = 1, num_subcyc
            if (tempo_cfgs % aerosolaware_flag) then
                nwfa_r1(:, :, 1) = nwfa(:, pver:1:-1)
                nwfa2d_1(:, 1) = nwfa2d(:)

                call tempo_aerosol_surface_emissions( &
                    dt=dtime_subcyc, nwfa=nwfa_r1, nwfa2d=nwfa2d_1, &
                    ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, kts=kts)

                new_nwfa_r1(:, :, 1) = nwfa_r1(:, :, 1)
            end if

            call tempo_run( &
                itimestep=curr_step, dt=dtime_subcyc, &
                dz=dz_r1, w=w_r1, p=p_r1, pii=exner_r1, th=new_th_r1, &
                qv=new_qv_r1, qc=new_qc_r1, qr=new_qr_r1, qi=new_qi_r1, qs=new_qs_r1, qg=new_qg_r1, nr=new_nr_r1, ni=new_ni_r1, &
                nc=new_nc_r1, nwfa=new_nwfa_r1, nifa=new_nifa_r1, &
                ng=new_ng_r1, qb=new_volg_r1, &
                ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
                ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, &
                tempo_cfgs=tempo_cfgs, &
                tempo_diags=tempo_driver_diags)

            ! TEMPO precipitation output is in mm. Convert to m and accumulate.
            rain(:) = rain(:) + &
                max(0.0_kind_phys, tempo_driver_diags % rain_precip(:, 1)) / 1000.0_kind_phys
            ice(:) = ice(:) + &
                max(0.0_kind_phys, tempo_driver_diags % ice_liquid_equiv_precip(:, 1)) / 1000.0_kind_phys
            snow(:) = snow(:) + ( &
                max(0.0_kind_phys, tempo_driver_diags % ice_liquid_equiv_precip(:, 1)) + &
                max(0.0_kind_phys, tempo_driver_diags % snow_liquid_equiv_precip(:, 1)) &
                ) / 1000.0_kind_phys
            graupel(:) = graupel(:) + &
                max(0.0_kind_phys, tempo_driver_diags % graupel_liquid_equiv_precip(:, 1)) / 1000.0_kind_phys
            prcp(:) = prcp(:) + ( &
                max(0.0_kind_phys, tempo_driver_diags % rain_precip(:, 1)) + &
                max(0.0_kind_phys, tempo_driver_diags % ice_liquid_equiv_precip(:, 1)) + &
                max(0.0_kind_phys, tempo_driver_diags % snow_liquid_equiv_precip(:, 1)) + &
                max(0.0_kind_phys, tempo_driver_diags % graupel_liquid_equiv_precip(:, 1)) &
                ) / 1000.0_kind_phys
        end do

        tend_th(:, :) = (new_th_r1(:, pver:1:-1, 1) - th(:, :)) / dtime_phys
        tend_qv(:, :) = (new_qv_r1(:, pver:1:-1, 1) - qv(:, :)) / dtime_phys
        tend_qc(:, :) = (new_qc_r1(:, pver:1:-1, 1) - qc(:, :)) / dtime_phys
        tend_qr(:, :) = (new_qr_r1(:, pver:1:-1, 1) - qr(:, :)) / dtime_phys
        tend_qi(:, :) = (new_qi_r1(:, pver:1:-1, 1) - qi(:, :)) / dtime_phys
        tend_qs(:, :) = (new_qs_r1(:, pver:1:-1, 1) - qs(:, :)) / dtime_phys
        tend_qg(:, :) = (new_qg_r1(:, pver:1:-1, 1) - qg(:, :)) / dtime_phys
        tend_nr(:, :) = (new_nr_r1(:, pver:1:-1, 1) - nr(:, :)) / dtime_phys
        tend_ni(:, :) = (new_ni_r1(:, pver:1:-1, 1) - ni(:, :)) / dtime_phys

        if (tempo_cfgs % aerosolaware_flag) then
            tend_nc(:, :) = (new_nc_r1(:, pver:1:-1, 1) - nc(:, :)) / dtime_phys
            tend_nwfa(:, :) = (new_nwfa_r1(:, pver:1:-1, 1) - nwfa(:, :)) / dtime_phys
            tend_nifa(:, :) = (new_nifa_r1(:, pver:1:-1, 1) - nifa(:, :)) / dtime_phys

            deallocate(nwfa_r1)
            deallocate(nwfa2d_1)
            deallocate(new_nc_r1)
            deallocate(new_nwfa_r1)
            deallocate(new_nifa_r1)
        end if

        if (tempo_cfgs % hailaware_flag) then
            tend_ng(:, :) = (new_ng_r1(:, pver:1:-1, 1) - ng(:, :)) / dtime_phys
            tend_volg(:, :) = (new_volg_r1(:, pver:1:-1, 1) - volg(:, :)) / dtime_phys

            deallocate(new_ng_r1)
            deallocate(new_volg_r1)
        end if

        frozen_fraction(:) = tempo_driver_diags % frozen_fraction(:, 1)
        refl10cm(:, :) = tempo_driver_diags % refl10cm(:, pver:1:-1, 1)
        ! These are in um, not m.
        re_cloud(:, :) = tempo_driver_diags % re_cloud(:, pver:1:-1, 1)
        re_ice(:, :) = tempo_driver_diags % re_ice(:, pver:1:-1, 1)
        re_snow(:, :) = tempo_driver_diags % re_snow(:, pver:1:-1, 1)
    end subroutine mp_tempo_run
end module mp_tempo
