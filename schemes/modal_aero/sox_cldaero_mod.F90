module sox_cldaero_mod
  use ccpp_kinds, only: kind_phys
  use cldaero_mod, only: cldaero_conc_t, cldaero_allocate, cldaero_deallocate
  use cldaero_mod, only: cldaero_uptakerate
  use aerosol_properties_mod, only: aerosol_properties
  use aerosol_state_mod, only: aerosol_state

  implicit none
  private

  public :: sox_cldaero_init
  public :: sox_cldaero_create_obj
  public :: sox_cldaero_update
  public :: sox_cldaero_destroy_obj

  integer :: id_msa = -1, id_h2so4 = -1, id_so2 = -1, id_h2o2 = -1, id_nh3 = -1

  real(kind_phys), parameter :: small_value = 1.e-20_kind_phys

  integer :: ncnst_tot = -huge(1) ! total number of mode number conc + mode species
  integer, public, protected :: nbins = 0

  class(aerosol_properties), pointer :: aero_props => null()

  logical :: has_msa = .false.

  ! host value of pi (passed at init for bit-for-bit consistency)
  real(kind_phys) :: pi = -huge(1._kind_phys)

  ! apply the aqueous sulfur oxidation update to the aerosol/gas state here;
  ! .false. when the chemistry package does its own in-cloud sulfur oxidation
  ! (e.g. GEOS-Chem), to avoid double counting
  logical :: do_aqueous_sulfur_chemistry_aerosol_update = .true.

contains

  subroutine sox_cldaero_init(aero_props_in, id_msa_in, id_h2so4_in, id_so2_in, &
                              id_h2o2_in, id_nh3_in, pi_in, do_aqueous_sulfur_chemistry_aerosol_update_in, &
                              errmsg, errflg)
    class(aerosol_properties), target, intent(in) :: aero_props_in
    ! species indices in the chemistry solution array, resolved by the host
    integer, intent(in) :: id_msa_in, id_h2so4_in, id_so2_in, id_h2o2_in, id_nh3_in
    real(kind_phys), intent(in) :: pi_in                ! host value of pi
    logical, intent(in) :: do_aqueous_sulfur_chemistry_aerosol_update_in
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    errmsg = ''
    errflg = 0

    id_msa = id_msa_in
    id_h2so4 = id_h2so4_in
    id_so2 = id_so2_in
    id_h2o2 = id_h2o2_in
    id_nh3 = id_nh3_in
    has_msa = id_msa > 0

    pi = pi_in
    do_aqueous_sulfur_chemistry_aerosol_update = do_aqueous_sulfur_chemistry_aerosol_update_in

    if (id_so2 < 1) then
      errflg = 1
      errmsg = 'sox_cldaero_init: SO2 is not included in chemistry -- should not invoke sox_cldaero_mod...'
      return
    end if

    aero_props => aero_props_in

    ncnst_tot = aero_props%ncnst_tot()
    nbins = aero_props%nbins()

  end subroutine sox_cldaero_init

  function sox_cldaero_create_obj(cldfrc, qcw, lwc, cfact, ncol, pver) result(conc_obj)

    real(kind_phys), intent(in) :: cldfrc(:, :)
    real(kind_phys), intent(in) :: qcw(:, :, :)
    real(kind_phys), intent(in) :: lwc(:, :)
    real(kind_phys), intent(in) :: cfact(:, :)
    integer, intent(in) :: ncol
    integer, intent(in) :: pver

    type(cldaero_conc_t), pointer :: conc_obj

    character(len=32) :: spectype
    integer :: l, m
    integer :: i, k, mm, ntot_amode
    logical :: mode7

    conc_obj => cldaero_allocate(ncol, pver)

    if (aero_props%model_is('BAM')) then
      ! no cloud-borne aerosols
      conc_obj%xlwc(:ncol, :) = lwc(:ncol, :)*cfact(:ncol, :) ! cloud water L(water)/L(air)
      return
    end if

    if (aero_props%model_is('MAM')) then
      ntot_amode = aero_props%nbins()
      if (ntot_amode /= 7) then
        conc_obj%so4_fact = 1._kind_phys
      end if
    end if

    do k = 1, pver
      do i = 1, ncol
        if (cldfrc(i, k) > 0._kind_phys) then
          conc_obj%xlwc(i, k) = lwc(i, k)*cfact(i, k) ! cloud water L(water)/L(air)
          conc_obj%xlwc(i, k) = conc_obj%xlwc(i, k)/cldfrc(i, k) ! liquid water in the cloudy fraction of cell
        else
          conc_obj%xlwc(i, k) = 0._kind_phys
        end if
      end do
    end do

    conc_obj%no3c(:, :) = 0._kind_phys
    conc_obj%nh4c(:, :) = 0._kind_phys
    conc_obj%so4c(:, :) = 0._kind_phys

    do k = 1, pver
      do i = 1, ncol
        do m = 1, aero_props%nbins()
          do l = 1, aero_props%nspecies(m)
            mm = aero_props%indexer(m, l)
            call aero_props%get(m, l, spectype=spectype)
            if (trim(spectype) == 'sulfate') then
              conc_obj%so4c(i, k) = conc_obj%so4c(i, k) + qcw(i, k, mm)
            end if
            if (trim(spectype) == 'ammonium') then
              conc_obj%nh4c(i, k) = conc_obj%nh4c(i, k) + qcw(i, k, mm)
            end if
          end do
        end do
      end do
    end do

  end function sox_cldaero_create_obj

  ! Update the mixing ratios
  subroutine sox_cldaero_update(aero_state, &
                                ncol, pver, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, xlwc, &
                                gravit, &
                                delso4_hprxn, xh2so4, xso4, xso4_init, nh3g, xnh3, xnh4c, xmsa, xso2, xh2o2, qcw, qin, &
                                aqso4, aqh2so4, aqso4_h2o2, aqso4_o3, aqso4_h2o2_3d, aqso4_o3_3d)

    class(aerosol_state), intent(in) :: aero_state

    integer, intent(in) :: ncol
    integer, intent(in) :: pver

    real(kind_phys), intent(in) :: dtime ! time step (sec)

    real(kind_phys), intent(in) :: mbar(:, :) ! mean wet atmospheric mass ( amu )
    real(kind_phys), intent(in) :: pdel(:, :)
    real(kind_phys), intent(in) :: press(:, :)
    real(kind_phys), intent(in) :: tfld(:, :)

    real(kind_phys), intent(in) :: cldnum(:, :)
    real(kind_phys), intent(in) :: cldfrc(:, :)
    real(kind_phys), intent(in) :: cfact(:, :)
    real(kind_phys), intent(in) :: xlwc(:, :)

    real(kind_phys), intent(in) :: gravit ! gravitational acceleration (m/s2)

    real(kind_phys), intent(in) :: delso4_hprxn(:, :)
    real(kind_phys), intent(in) :: xh2so4(:, :)
    real(kind_phys), intent(in) :: xso4(:, :)
    real(kind_phys), intent(in) :: xso4_init(:, :)
    real(kind_phys), intent(in) :: nh3g(:, :)
    real(kind_phys), intent(in) :: xnh3(:, :)
    real(kind_phys), intent(in) :: xnh4c(:, :)
    real(kind_phys), intent(in) :: xmsa(:, :)
    real(kind_phys), intent(in) :: xso2(:, :)
    real(kind_phys), intent(in) :: xh2o2(:, :)

    real(kind_phys), intent(inout) :: qcw(:, :, :) ! cloud-borne aerosol (vmr)
    real(kind_phys), intent(inout) :: qin(:, :, :) ! xported species ( vmr )

    real(kind_phys), intent(out) :: aqso4(:, :)                   ! aqueous phase chemistry
    real(kind_phys), intent(out) :: aqh2so4(:, :)                 ! aqueous phase chemistry
    real(kind_phys), intent(out) :: aqso4_h2o2(:)                ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(kind_phys), intent(out) :: aqso4_o3(:)                  ! SO4 aqueous phase chemistry due to O3 (kg/m2)
    real(kind_phys), intent(out), optional :: aqso4_h2o2_3d(:, :)                ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(kind_phys), intent(out), optional :: aqso4_o3_3d(:, :)                  ! SO4 aqueous phase chemistry due to O3 (kg/m2)

    ! local vars ...

    real(kind_phys) :: dqdt_aqso4(ncol, pver, ncnst_tot), &
                       dqdt_aqh2so4(ncol, pver, ncnst_tot), &
                       dqdt_aqhprxn(ncol, pver), dqdt_aqo3rxn(ncol, pver)

    real(kind_phys) :: faqgain_msa(nbins, ncol, pver), faqgain_so4(nbins, ncol, pver)
    real(kind_phys) :: delso4_ox(ncol, pver)

    real(kind_phys) :: delnh3, delnh4
    real(kind_phys) :: dso4dt_aqrxn, dso4dt_hprxn, &
                       dso4dt_gasuptk, dmsadt_gasuptk, &
                       dmsadt_gasuptk_tomsa, dmsadt_gasuptk_toso4, &
                       dqdt_aq, dqdt_wr, dqdt

    real(kind_phys) :: fwetrem, uptkrate

    integer :: l, m, n, mm
    integer :: i, k
    real(kind_phys) :: xl
    real(kind_phys) :: mw_so4
    character(len=32) :: spectype
    character(len=32) :: specname

    ! make sure dqdt is zero initially, for budgets
    dqdt_aqso4(:, :, :) = 0.0_kind_phys
    dqdt_aqh2so4(:, :, :) = 0.0_kind_phys
    dqdt_aqhprxn(:, :) = 0.0_kind_phys
    dqdt_aqo3rxn(:, :) = 0.0_kind_phys

    aqso4 = 0.0_kind_phys
    aqh2so4 = 0.0_kind_phys
    aqso4_h2o2 = 0.0_kind_phys
    aqso4_o3 = 0.0_kind_phys
    delso4_ox = 0.0_kind_phys

    ! Avoid double counting in-cloud sulfur oxidation when running with
    ! GEOS-Chem. If running with GEOS-Chem then sulfur oxidation
    ! is performed internally to GEOS-Chem. Here, we just return to the
    ! parent routine and thus we do not apply tendencies calculated by MAM.
    ! (The host sets this flag; CAM passes .not. cam_chempkg_is('geoschem_mam4').)
    if (.not. do_aqueous_sulfur_chemistry_aerosol_update) return

    where (cldfrc(:ncol, :) >= 1.0e-5_kind_phys)
      delso4_ox(:ncol, :) = xso4(:ncol, :) - xso4_init(:ncol, :)
    end where

    !-------------------------------------------------------------------------
    ! Compute factors for partitioning aerosol mass gains among bins / modes.
    ! The factors are proportional to the activated particle MR for each
    ! bin, which is the MR of cloud drops "associated with" the mode
    ! thus we are assuming the cloud drop size is independent of the
    ! associated aerosol mode properties
    call aero_state%aqu_gain_binfraction(aero_props, 'sulfate', qcw, delso4_ox, faqgain_so4)
    if (has_msa) call aero_state%aqu_gain_binfraction(aero_props, 'msa', qcw, delso4_ox, faqgain_msa)

    lev_loop: do k = 1, pver
      col_loop: do i = 1, ncol
        cloud: if (cldfrc(i, k) >= 1.0e-5_kind_phys) then
          xl = xlwc(i, k)

          if (xl >= 1.e-8_kind_phys) then !! when cloud is present

            if (id_nh3 > 0) then
              delnh3 = nh3g(i, k) - xnh3(i, k)
              delnh4 = -delnh3
            end if

            ! faqgain_msa(n) = fraction of total msa_c gain going to mode n

            uptkrate = cldaero_uptakerate(xl, cldnum(i, k), cfact(i, k), cldfrc(i, k), tfld(i, k), press(i, k), pi)
            ! average uptake rate over dtime
            uptkrate = (1.0_kind_phys - exp(-min(100._kind_phys, dtime*uptkrate)))/dtime

            ! dso4dt_gasuptk = so4_c tendency from h2so4 gas uptake (mol/mol/s)
            ! dmsadt_gasuptk = msa_c tendency from msa gas uptake (mol/mol/s)
            dso4dt_gasuptk = xh2so4(i, k)*uptkrate
            if (has_msa) then
              dmsadt_gasuptk = xmsa(i, k)*uptkrate
            else
              dmsadt_gasuptk = 0.0_kind_phys
            end if

            ! if no modes have msa aerosol, then "rename" scavenged msa gas to so4
            if (has_msa) then
              dmsadt_gasuptk_toso4 = 0.0_kind_phys
              dmsadt_gasuptk_tomsa = dmsadt_gasuptk
            else
              ! no MSA
              dmsadt_gasuptk_tomsa = 0.0_kind_phys
              dmsadt_gasuptk_toso4 = dmsadt_gasuptk
            end if

            !-----------------------------------------------------------------------
            ! now compute TMR tendencies
            ! this includes the above aqueous so2 chemistry AND
            ! the uptake of highly soluble aerosol precursor gases (h2so4, msa, ...)
            ! AND the wetremoval of dissolved, unreacted so2 and h2o2

            dso4dt_aqrxn = (delso4_ox(i, k) + delso4_hprxn(i, k))/dtime
            dso4dt_hprxn = delso4_hprxn(i, k)/dtime

            ! fwetrem = fraction of in-cloud-water material that is wet removed
            ! fwetrem = max( 0.0_kind_phys, (1.0_kind_phys-exp(-min(100._kind_phys,dtime*clwlrat(i,k)))) )
            fwetrem = 0.0_kind_phys ! don't have so4 & msa wet removal here

            ! compute TMR tendencies for so4 and msa aerosol-in-cloud-water
            do m = 1, aero_props%nbins()
              do l = 1, aero_props%nspecies(m)
                mm = aero_props%indexer(m, l)
                call aero_props%get(m, l, spectype=spectype)
                if (trim(spectype) == 'sulfate') then

                  dqdt_aqso4(i, k, mm) = faqgain_so4(m, i, k)*dso4dt_aqrxn*cldfrc(i, k)

                  dqdt_aqh2so4(i, k, mm) = faqgain_so4(m, i, k)* &
                                           (dso4dt_gasuptk + dmsadt_gasuptk_toso4)*cldfrc(i, k)
                  dqdt_aq = dqdt_aqso4(i, k, mm) + dqdt_aqh2so4(i, k, mm)
                  dqdt_wr = -fwetrem*dqdt_aq
                  dqdt = dqdt_aq + dqdt_wr
                  qcw(i, k, mm) = qcw(i, k, mm) + dqdt*dtime

                end if
                if (trim(spectype) == 'msa') then
                  dqdt_aq = faqgain_msa(m, i, k)*dmsadt_gasuptk_tomsa*cldfrc(i, k)
                  dqdt_wr = -fwetrem*dqdt_aq
                  dqdt = dqdt_aq + dqdt_wr
                  qcw(i, k, mm) = qcw(i, k, mm) + dqdt*dtime
                end if
                if (trim(spectype) == 'ammonium') then
                  if (delnh4 > 0.0_kind_phys) then
                    dqdt_aq = faqgain_so4(m, i, k)*delnh4/dtime*cldfrc(i, k)
                    dqdt = dqdt_aq
                    qcw(i, k, mm) = qcw(i, k, mm) + dqdt*dtime
                  else
                    dqdt = (qcw(i, k, mm)/max(xnh4c(i, k), 1.0e-35_kind_phys)) &
                           *delnh4/dtime*cldfrc(i, k)
                    qcw(i, k, mm) = qcw(i, k, mm) + dqdt*dtime
                  end if
                end if
              end do
            end do

            ! For gas species, tendency includes
            ! reactive uptake to cloud water that essentially transforms the gas to
            ! a different species. Wet removal associated with this is applied
            ! to the "new" species (e.g., so4_c) rather than to the gas.
            ! wet removal of the unreacted gas that is dissolved in cloud water.
            ! Need to multiply both these parts by cldfrc

            ! h2so4 (g) & msa (g)
            qin(i, k, id_h2so4) = qin(i, k, id_h2so4) - dso4dt_gasuptk*dtime*cldfrc(i, k)
            if (has_msa) qin(i, k, id_msa) = qin(i, k, id_msa) - dmsadt_gasuptk*dtime*cldfrc(i, k)

            ! so2 -- the first order loss rate for so2 is frso2_c*clwlrat(i,k)
            ! fwetrem = max( 0.0_kind_phys, (1.0_kind_phys-exp(-min(100._kind_phys,dtime*frso2_c*clwlrat(i,k)))) )
            fwetrem = 0.0_kind_phys ! don't include so2 wet removal here

            dqdt_wr = -fwetrem*xso2(i, k)/dtime*cldfrc(i, k)
            dqdt_aq = -dso4dt_aqrxn*cldfrc(i, k)
            dqdt = dqdt_aq + dqdt_wr
            qin(i, k, id_so2) = qin(i, k, id_so2) + dqdt*dtime

            ! h2o2 -- the first order loss rate for h2o2 is frh2o2_c*clwlrat(i,k)
            ! fwetrem = max( 0.0_kind_phys, (1.0_kind_phys-exp(-min(100._kind_phys,dtime*frh2o2_c*clwlrat(i,k)))) )
            fwetrem = 0.0_kind_phys ! don't include h2o2 wet removal here

            dqdt_wr = -fwetrem*xh2o2(i, k)/dtime*cldfrc(i, k)
            dqdt_aq = -dso4dt_hprxn*cldfrc(i, k)
            dqdt = dqdt_aq + dqdt_wr
            qin(i, k, id_h2o2) = qin(i, k, id_h2o2) + dqdt*dtime

            ! NH3
            if (id_nh3 > 0) then
              dqdt_aq = delnh3/dtime*cldfrc(i, k)
              dqdt = dqdt_aq
              qin(i, k, id_nh3) = qin(i, k, id_nh3) + dqdt*dtime
            end if

            ! for SO4 from H2O2/O3 budgets
            dqdt_aqhprxn(i, k) = dso4dt_hprxn*cldfrc(i, k)
            dqdt_aqo3rxn(i, k) = (dso4dt_aqrxn - dso4dt_hprxn)*cldfrc(i, k)

          end if !! when cloud is present
        end if cloud
      end do col_loop
    end do lev_loop

    !==============================================================
    ! ... Update the mixing ratios
    !==============================================================
    do k = 1, pver

      do n = 1, aero_props%nbins()
        do l = 1, aero_props%nspecies(n)
          mm = aero_props%indexer(n, l)
          call aero_props%get(n, l, spectype=spectype)
          if (trim(spectype) == 'sulfate') then
            qcw(:ncol, k, mm) = MAX(qcw(:ncol, k, mm), small_value)
          end if
          if (trim(spectype) == 'msa') then
            qcw(:ncol, k, mm) = MAX(qcw(:ncol, k, mm), small_value)
          end if
          if (trim(spectype) == 'ammonium') then
            qcw(:ncol, k, mm) = MAX(qcw(:ncol, k, mm), small_value)
          end if
        end do
      end do

      qin(:ncol, k, id_so2) = MAX(qin(:ncol, k, id_so2), small_value)
      qin(:ncol, k, id_h2o2) = MAX(qin(:ncol, k, id_h2o2), small_value)
      qin(:ncol, k, id_h2so4) = MAX(qin(:ncol, k, id_h2so4), small_value)
      if (id_msa > 0) qin(:ncol, k, id_msa) = MAX(qin(:ncol, k, id_msa), small_value)
      if (id_nh3 > 0) qin(:ncol, k, id_nh3) = MAX(qin(:ncol, k, id_nh3), small_value)

    end do

    ! diagnostics
    mw_so4 = -huge(1._kind_phys)

    do n = 1, aero_props%nbins()
      ! while looking through all species, only dqdt_aqso4 from sulfates  is gt zero
      do l = 1, aero_props%nspecies(n)
        mm = aero_props%indexer(n, l)
        call aero_props%get(n, l, spectype=spectype, specname=specname)
        if (trim(spectype) == 'sulfate') then
          call aero_props%get(n, l, spec_mw=mw_so4)
          aqso4(:, n) = 0._kind_phys
          do k = 1, pver
            do i = 1, ncol
              aqso4(i, n) = aqso4(i, n) + dqdt_aqso4(i, k, mm)*mw_so4/mbar(i, k) &
                            *pdel(i, k)/gravit ! kg/m2/s
            end do
          end do
          aqh2so4(:, n) = 0._kind_phys
          do k = 1, pver
            do i = 1, ncol
              aqh2so4(i, n) = aqh2so4(i, n) + dqdt_aqh2so4(i, k, mm)*mw_so4/mbar(i, k) &
                              *pdel(i, k)/gravit ! kg/m2/s
            end do
          end do
        end if
      end do
    end do

    aqso4_h2o2(:) = 0._kind_phys
    do k = 1, pver
      do i = 1, ncol
        aqso4_h2o2(i) = aqso4_h2o2(i) + dqdt_aqhprxn(i, k)*mw_so4/mbar(i, k) &
                        *pdel(i, k)/gravit ! kg SO4 /m2/s
      end do
    end do

    if (present(aqso4_h2o2_3d)) then
      aqso4_h2o2_3d(:, :) = 0._kind_phys
      do k = 1, pver
        do i = 1, ncol
          aqso4_h2o2_3d(i, k) = dqdt_aqhprxn(i, k)*mw_so4/mbar(i, k) &
                                *pdel(i, k)/gravit ! kg SO4 /m2/s
        end do
      end do
    end if

    aqso4_o3(:) = 0._kind_phys
    do k = 1, pver
      do i = 1, ncol
        aqso4_o3(i) = aqso4_o3(i) + dqdt_aqo3rxn(i, k)*mw_so4/mbar(i, k) &
                      *pdel(i, k)/gravit ! kg SO4 /m2/s
      end do
    end do

    if (present(aqso4_o3_3d)) then
      aqso4_o3_3d(:, :) = 0._kind_phys
      do k = 1, pver
        do i = 1, ncol
          aqso4_o3_3d(i, k) = dqdt_aqo3rxn(i, k)*mw_so4/mbar(i, k) &
                              *pdel(i, k)/gravit ! kg SO4 /m2/s
        end do
      end do
    end if

  end subroutine sox_cldaero_update

  subroutine sox_cldaero_destroy_obj(conc_obj)
    type(cldaero_conc_t), pointer :: conc_obj

    call cldaero_deallocate(conc_obj)

  end subroutine sox_cldaero_destroy_obj

end module sox_cldaero_mod
