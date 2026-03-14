! Aerosol optical properties for RRTMGP radiation.
! Computes per-band SW/LW aerosol optics from aerosol model instances.
!
! Author: Haipeng Lin, NSF-NCAR/CGD/AMP, March 2026
module aerosol_optics
  implicit none
  private

  public :: aerosol_optics_init
  public :: aerosol_optics_run

contains

!> \section arg_table_aerosol_optics_init Argument Table
!! \htmlinclude aerosol_optics_init.html
  subroutine aerosol_optics_init(N_DIAG, active_calls, &
                                 num_bulk_aer, &
                                 errmsg, errflg)
    use aerosol_mmr_ccpp,                only: rad_aer_diag_init
    use radiative_aerosol_definitions,   only: aerlist_t
    use radiative_aerosol_definitions,   only: bulk_aerosol_list

    integer,            intent(in)  :: N_DIAG
    logical,            intent(in)  :: active_calls(:)
    !type(aerlist_t),    intent(in)  :: bulk_aerosol_list(:)
    ! does not work: errors with
    ! parse_source.CCPPError: No ddt_lib or ddt aerlist_t not in ddt_lib

    integer,            intent(out) :: num_bulk_aer
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: i

    errmsg = ''
    errflg = 0

    ! Set number of bulk aerosol constituents in climate list
    num_bulk_aer = bulk_aerosol_list(0)%numaerosols

    ! Register aerosol diagnostic history fields for active lists
    do i = 1, N_DIAG
      if (active_calls(i)) then
        call rad_aer_diag_init(bulk_aerosol_list(i))
      end if
    end do

  end subroutine aerosol_optics_init

!> \section arg_table_aerosol_optics_run Argument Table
!! \htmlinclude aerosol_optics_run.html
  subroutine aerosol_optics_run( &
    ncol, pver, &
    const_props, &
    nswbands, nlwbands, &
    top_lev, &
    rga, &
    idx_sw_diag, &
    num_bulk_aer, &
    relh, pdeldry, constituents, &
    crefwsw, crefwlw, &
    aer_tau, aer_tau_w, aer_tau_w_g, aer_lw_abs, &
    dustaod, sulfaod, bcaod, pomaod, soaaod, ssltaod, &
    aodabsbc, &
    burdendust, burdenso4, burdenbc, burdenpom, burdensoa, burdenseasalt, &
    ssavis, aodvis, &
    odv_col_aod, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    ! host-model dependency for aerosol objects:
    use aerosol_physical_properties,     only: nrh
    use aerosol_optics_core,             only: aerosol_optics_sw_bin, aerosol_optics_lw_bin
    use aerosol_instances_mod,           only: aerosol_instances_get_props, &
                                               aerosol_instances_get_state, &
                                               aerosol_instances_get_num_models
    use aerosol_properties_mod,          only: aerosol_properties
    use aerosol_state_mod,               only: aerosol_state
    use aerosol_physical_properties,     only: ot_length
    use aerosol_mmr_ccpp,                only: rad_aer_diag_out
    use radiative_aerosol_definitions,   only: N_DIAG, active_calls, bulk_aerosol_list

    ! Input arguments
    integer,          intent(in)  :: ncol                  ! number of columns [count]
    integer,          intent(in)  :: pver                  ! number of vertical layers [count]
    type(ccpp_constituent_prop_ptr_t), &
                      intent(in)  :: const_props(:)        ! ccpp constituent properties pointer
    integer,          intent(in)  :: nswbands              ! number of SW bands [count]
    integer,          intent(in)  :: nlwbands              ! number of LW bands [count]
    integer,          intent(in)  :: top_lev               ! top level for aerosol model [index]
    real(kind_phys),  intent(in)  :: rga
    integer,          intent(in)  :: idx_sw_diag           ! index of SW diagnostic band [index]
    integer,          intent(in)  :: num_bulk_aer          ! number of bulk aerosol constituents in climate list [count]
    real(kind_phys),  intent(in)  :: relh(:, :)            ! relative humidity [fraction]
    real(kind_phys),  intent(in)  :: pdeldry(:, :)         ! dry air pressure thickness [Pa]
    real(kind_phys),  intent(in)  :: constituents(:, :, :)
    complex(kind_phys), intent(in) :: crefwsw(:)           ! water refractive idx for SW rad [1]
    complex(kind_phys), intent(in) :: crefwlw(:)           ! water refractive idx for LW rad [1]

    ! Output arguments
    real(kind_phys),  intent(out) :: aer_tau(:, :, :)      ! SW extinction optical depth [1]
    real(kind_phys),  intent(out) :: aer_tau_w(:, :, :)    ! SW ssa * tau [1]
    real(kind_phys),  intent(out) :: aer_tau_w_g(:, :, :)  ! SW asy * ssa * tau [1]
    real(kind_phys),  intent(out) :: aer_lw_abs(:, :, :)   ! LW absorption optical depth [1]

    ! Per-species diagnostic outputs (all ncol)
    real(kind_phys),  intent(out) :: dustaod(:)             ! Dust AOD at vis [1]
    real(kind_phys),  intent(out) :: sulfaod(:)             ! Sulfate AOD at vis [1]
    real(kind_phys),  intent(out) :: bcaod(:)               ! BC AOD at vis [1]
    real(kind_phys),  intent(out) :: pomaod(:)              ! POM AOD at vis [1]
    real(kind_phys),  intent(out) :: soaaod(:)              ! SOA AOD at vis [1]
    real(kind_phys),  intent(out) :: ssltaod(:)             ! Seasalt AOD at vis [1]
    real(kind_phys),  intent(out) :: aodabsbc(:)            ! BC absorption OD at vis [1]
    real(kind_phys),  intent(out) :: burdendust(:)          ! Dust column burden [kg m-2]
    real(kind_phys),  intent(out) :: burdenso4(:)           ! Sulfate column burden [kg m-2]
    real(kind_phys),  intent(out) :: burdenbc(:)            ! BC column burden [kg m-2]
    real(kind_phys),  intent(out) :: burdenpom(:)           ! POM column burden [kg m-2]
    real(kind_phys),  intent(out) :: burdensoa(:)           ! SOA column burden [kg m-2]
    real(kind_phys),  intent(out) :: burdenseasalt(:)       ! Seasalt column burden [kg m-2]
    real(kind_phys),  intent(out) :: ssavis(:)              ! Column SSA at vis [1]
    real(kind_phys),  intent(out) :: aodvis(:)              ! Column AOD at vis [1]

    ! Per-constituent column visible OD for bulk aerosol (ncol, num_bulk_aer)
    real(kind_phys),  intent(out) :: odv_col_aod(:,:)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer  :: iaermod, ibin, nbins, iwav, ilev, icol, ilist
    integer  :: num_aero_models, numrh
    integer  :: bam_cnt
    real(kind_phys) :: mass(ncol, pver)

    ! Index into constituents array for VOLC_RAD_GEOM
    integer :: idx_volc_rad_geom = -1

    ! Aerosol object pointers
    class(aerosol_properties), pointer :: aeroprops
    class(aerosol_state),      pointer :: aerostate

    real(kind_phys) :: sulfwtpct(ncol, pver)

    ! Per-bin SW optics
    real(kind_phys) :: tau_bin(ncol, pver, nswbands)
    real(kind_phys) :: ssa_bin(ncol, pver, nswbands)
    real(kind_phys) :: asm_bin(ncol, pver, nswbands)
    real(kind_phys) :: pabs_vis(ncol, pver)
    real(kind_phys) :: dopaer0_vis(ncol, pver)

    ! Per-bin LW optics
    real(kind_phys) :: tau_lw_bin(ncol, pver, nlwbands)
    real(kind_phys) :: absorp_bin(ncol, pver, nlwbands)

    ! Volcanic geometric radius (pointer into constituents)
    real(kind_phys), target  :: geometric_radius(ncol, pver)
    logical :: has_volc_rad_geom

    character(len=ot_length) :: opticstype
    character(len=*), parameter :: subname = 'aerosol_optics_run'

    ! Per-species diagnostic locals
    integer :: ispec
    real(kind_phys) :: wetvol_d(ncol, pver), watervol_d(ncol, pver)
    real(kind_phys) :: vol_d(ncol)
    real(kind_phys) :: specdens_d, hygro_aer_d
    character(len=32) :: spectype_d
    real(kind_phys), pointer :: specmmr_d(:, :)
    complex(kind_phys), pointer :: specrefindex_d(:)
    real(kind_phys) :: scatdust_d(ncol), absdust_d(ncol), hygrodust_d(ncol)
    real(kind_phys) :: scatbc_d(ncol), absbc_d(ncol), hygrobc_d(ncol)
    real(kind_phys) :: scatsulf_d(ncol), abssulf_d(ncol), hygrosulf_d(ncol)
    real(kind_phys) :: scatpom_d(ncol), abspom_d(ncol), hygropom_d(ncol)
    real(kind_phys) :: scatsoa_d(ncol), abssoa_d(ncol), hygrosoa_d(ncol)
    real(kind_phys) :: scatsslt_d(ncol), abssslt_d(ncol), hygrosslt_d(ncol)
    real(kind_phys) :: scath2o_d, absh2o_d, sumscat_d, sumabs_d, sumhygro_d
    real(kind_phys) :: dopaer_d, palb_d, aodc_d

    errmsg = ''
    errflg = 0

    ! Initialize outputs
    aer_tau     = 0._kind_phys
    aer_tau_w   = 0._kind_phys
    aer_tau_w_g = 0._kind_phys
    aer_lw_abs  = 0._kind_phys

    ! Initialize per-constituent visible OD output for BAM
    bam_cnt = 0
    odv_col_aod(:,:) = 0._kind_phys

    ! Initialize per-species diagnostic outputs
    dustaod(:ncol)       = 0._kind_phys
    sulfaod(:ncol)       = 0._kind_phys
    bcaod(:ncol)         = 0._kind_phys
    pomaod(:ncol)        = 0._kind_phys
    soaaod(:ncol)        = 0._kind_phys
    ssltaod(:ncol)       = 0._kind_phys
    aodabsbc(:ncol)      = 0._kind_phys
    burdendust(:ncol)    = 0._kind_phys
    burdenso4(:ncol)     = 0._kind_phys
    burdenbc(:ncol)      = 0._kind_phys
    burdenpom(:ncol)     = 0._kind_phys
    burdensoa(:ncol)     = 0._kind_phys
    burdenseasalt(:ncol) = 0._kind_phys
    ssavis(:ncol)        = 0._kind_phys
    aodvis(:ncol)        = 0._kind_phys

    ! Layer mass [kg/m2]
    mass(:ncol, :) = pdeldry(:ncol, :) * rga

    ! Number of relative humidity bins for table lookup
    numrh = nrh

    ! Check if VOLC_RAD_GEOM is available
    call ccpp_const_get_idx(const_props, &
                            'VOLC_RAD_GEOM', idx_volc_rad_geom, errmsg, errflg)
    if(errflg /= 0) return

    has_volc_rad_geom = (idx_volc_rad_geom > 0)
    if (has_volc_rad_geom) then
      geometric_radius(:ncol, :pver) = constituents(:ncol, :pver, idx_volc_rad_geom)
    end if

    num_aero_models = aerosol_instances_get_num_models()

    ! Loop over aerosol models in the climate list.
    aerosol_model_loop: do iaermod = 1, num_aero_models

      aeroprops => aerosol_instances_get_props(iaermod, list_idx=0)
      aerostate => aerosol_instances_get_state(iaermod, list_idx=0)

      if (.not. associated(aeroprops)) cycle
      if (.not. associated(aerostate)) cycle

      nbins = aeroprops%nbins()

      ! Sulfate weight percent for hygroscopic_wtp optics
      sulfwtpct(:ncol, :pver) = aerostate%wgtpct(ncol, pver)

      ! Loop over aerosol bins/species
      bin_loop: do ibin = 1, nbins

        ! Check if volcanic_radius optics needs geometric_radius
        call aeroprops%optics_params(bin_ndx=ibin, opticstype=opticstype)
        if (index(opticstype, 'volcanic_radius') > 0 .and. .not. has_volc_rad_geom) then
          errflg = 1
          errmsg = 'VOLC_RAD_GEOM constituent required for volcanic_radius optics but not found'
          return
        end if

        !-------------------------------------------------
        ! Call the portable core for shortwave optics calculation per-bin:
        !-------------------------------------------------
        call aerosol_optics_sw_bin(aeroprops, aerostate, ibin, &
                                   ncol, pver, top_lev, nswbands, nlwbands, numrh, &
                                   idx_sw_diag, &
                                   relh, sulfwtpct, mass, crefwsw, crefwlw, &
                                   geometric_radius, &
                                   tau_bin, ssa_bin, asm_bin, &
                                   pabs_vis, dopaer0_vis, &
                                   errmsg, errflg)
        if (errflg /= 0) return

        ! Accumulate weighted products
        wavelength: do iwav = 1, nswbands
          vertical: do ilev = top_lev, pver
            column: do icol = 1, ncol
              aer_tau(icol, ilev, iwav)     = aer_tau(icol, ilev, iwav) &
                                              + tau_bin(icol, ilev, iwav)
              aer_tau_w(icol, ilev, iwav)   = aer_tau_w(icol, ilev, iwav) &
                                              + tau_bin(icol, ilev, iwav) * ssa_bin(icol, ilev, iwav)
              aer_tau_w_g(icol, ilev, iwav) = aer_tau_w_g(icol, ilev, iwav) &
                                              + tau_bin(icol, ilev, iwav) * ssa_bin(icol, ilev, iwav) &
                                                * asm_bin(icol, ilev, iwav)
            end do column
          end do vertical
        end do wavelength

        !-------------------------------------------------
        ! Per-species diagnostic accumulation at vis band
        !-------------------------------------------------
        wetvol_d(:ncol, :) = aerostate%wet_volume(aeroprops, ibin, ncol, pver)
        watervol_d(:ncol, :) = aerostate%water_volume(aeroprops, ibin, ncol, pver)

        do ilev = top_lev, pver
          ! Reset per-species accumulators for this bin/level
          scatdust_d(:ncol) = 0._kind_phys; absdust_d(:ncol) = 0._kind_phys; hygrodust_d(:ncol) = 0._kind_phys
          scatbc_d(:ncol) = 0._kind_phys; absbc_d(:ncol) = 0._kind_phys; hygrobc_d(:ncol) = 0._kind_phys
          scatsulf_d(:ncol) = 0._kind_phys; abssulf_d(:ncol) = 0._kind_phys; hygrosulf_d(:ncol) = 0._kind_phys
          scatpom_d(:ncol) = 0._kind_phys; abspom_d(:ncol) = 0._kind_phys; hygropom_d(:ncol) = 0._kind_phys
          scatsoa_d(:ncol) = 0._kind_phys; abssoa_d(:ncol) = 0._kind_phys; hygrosoa_d(:ncol) = 0._kind_phys
          scatsslt_d(:ncol) = 0._kind_phys; abssslt_d(:ncol) = 0._kind_phys; hygrosslt_d(:ncol) = 0._kind_phys

          do ispec = 1, aeroprops%nspecies(ibin)
            call aeroprops%get(ibin, ispec, density=specdens_d, &
                               spectype=spectype_d, refindex_sw=specrefindex_d, hygro=hygro_aer_d)
            call aerostate%get_ambient_mmr(species_ndx=ispec, bin_ndx=ibin, mmr=specmmr_d)

            do icol = 1, ncol
              vol_d(icol) = specmmr_d(icol, ilev) / specdens_d

              select case (trim(spectype_d))
              case('dust')
                burdendust(icol) = burdendust(icol) + specmmr_d(icol, ilev)*mass(icol, ilev)
                if (associated(specrefindex_d)) then
                  scatdust_d(icol) = vol_d(icol)*specrefindex_d(idx_sw_diag)%re
                  absdust_d(icol)  = -vol_d(icol)*specrefindex_d(idx_sw_diag)%im
                end if
                hygrodust_d(icol) = vol_d(icol)*hygro_aer_d
              case('black-c')
                burdenbc(icol) = burdenbc(icol) + specmmr_d(icol, ilev)*mass(icol, ilev)
                if (associated(specrefindex_d)) then
                  scatbc_d(icol) = vol_d(icol)*specrefindex_d(idx_sw_diag)%re
                  absbc_d(icol)  = -vol_d(icol)*specrefindex_d(idx_sw_diag)%im
                end if
                hygrobc_d(icol) = vol_d(icol)*hygro_aer_d
              case('sulfate')
                burdenso4(icol) = burdenso4(icol) + specmmr_d(icol, ilev)*mass(icol, ilev)
                if (associated(specrefindex_d)) then
                  scatsulf_d(icol) = vol_d(icol)*specrefindex_d(idx_sw_diag)%re
                  abssulf_d(icol)  = -vol_d(icol)*specrefindex_d(idx_sw_diag)%im
                end if
                hygrosulf_d(icol) = vol_d(icol)*hygro_aer_d
              case('p-organic')
                burdenpom(icol) = burdenpom(icol) + specmmr_d(icol, ilev)*mass(icol, ilev)
                if (associated(specrefindex_d)) then
                  scatpom_d(icol) = vol_d(icol)*specrefindex_d(idx_sw_diag)%re
                  abspom_d(icol)  = -vol_d(icol)*specrefindex_d(idx_sw_diag)%im
                end if
                hygropom_d(icol) = vol_d(icol)*hygro_aer_d
              case('s-organic')
                burdensoa(icol) = burdensoa(icol) + specmmr_d(icol, ilev)*mass(icol, ilev)
                if (associated(specrefindex_d)) then
                  scatsoa_d(icol) = vol_d(icol)*specrefindex_d(idx_sw_diag)%re
                  abssoa_d(icol)  = -vol_d(icol)*specrefindex_d(idx_sw_diag)%im
                end if
                hygrosoa_d(icol) = vol_d(icol)*hygro_aer_d
              case('seasalt')
                burdenseasalt(icol) = burdenseasalt(icol) + specmmr_d(icol, ilev)*mass(icol, ilev)
                if (associated(specrefindex_d)) then
                  scatsslt_d(icol) = vol_d(icol)*specrefindex_d(idx_sw_diag)%re
                  abssslt_d(icol)  = -vol_d(icol)*specrefindex_d(idx_sw_diag)%im
                end if
                hygrosslt_d(icol) = vol_d(icol)*hygro_aer_d
              end select
            end do ! icol
          end do ! ispec

          ! Volume-based partitioning and per-species AOD accumulation
          do icol = 1, ncol
            if (wetvol_d(icol, ilev) > 1.e-40_kind_phys .and. vol_d(icol) > 0._kind_phys) then
              scath2o_d  = watervol_d(icol, ilev)*crefwsw(idx_sw_diag)%re
              absh2o_d   = -watervol_d(icol, ilev)*crefwsw(idx_sw_diag)%im
              sumscat_d  = scatsulf_d(icol) + scatpom_d(icol) + scatsoa_d(icol) + scatbc_d(icol) + &
                           scatdust_d(icol) + scatsslt_d(icol) + scath2o_d
              sumabs_d   = abssulf_d(icol) + abspom_d(icol) + abssoa_d(icol) + absbc_d(icol) + &
                           absdust_d(icol) + abssslt_d(icol) + absh2o_d
              sumhygro_d = hygrosulf_d(icol) + hygropom_d(icol) + hygrosoa_d(icol) + hygrobc_d(icol) + &
                           hygrodust_d(icol) + hygrosslt_d(icol)

              ! Partition water contribution proportionally to hygroscopicity
              scatdust_d(icol) = (scatdust_d(icol) + scath2o_d*hygrodust_d(icol)/sumhygro_d)/sumscat_d
              absdust_d(icol)  = (absdust_d(icol) + absh2o_d*hygrodust_d(icol)/sumhygro_d)/sumabs_d

              scatsulf_d(icol) = (scatsulf_d(icol) + scath2o_d*hygrosulf_d(icol)/sumhygro_d)/sumscat_d
              abssulf_d(icol)  = (abssulf_d(icol) + absh2o_d*hygrosulf_d(icol)/sumhygro_d)/sumabs_d

              scatpom_d(icol) = (scatpom_d(icol) + scath2o_d*hygropom_d(icol)/sumhygro_d)/sumscat_d
              abspom_d(icol)  = (abspom_d(icol) + absh2o_d*hygropom_d(icol)/sumhygro_d)/sumabs_d

              scatsoa_d(icol) = (scatsoa_d(icol) + scath2o_d*hygrosoa_d(icol)/sumhygro_d)/sumscat_d
              abssoa_d(icol)  = (abssoa_d(icol) + absh2o_d*hygrosoa_d(icol)/sumhygro_d)/sumabs_d

              scatbc_d(icol) = (scatbc_d(icol) + scath2o_d*hygrobc_d(icol)/sumhygro_d)/sumscat_d
              absbc_d(icol)  = (absbc_d(icol) + absh2o_d*hygrobc_d(icol)/sumhygro_d)/sumabs_d

              scatsslt_d(icol) = (scatsslt_d(icol) + scath2o_d*hygrosslt_d(icol)/sumhygro_d)/sumscat_d
              abssslt_d(icol)  = (abssslt_d(icol) + absh2o_d*hygrosslt_d(icol)/sumhygro_d)/sumabs_d

              dopaer_d = tau_bin(icol, ilev, idx_sw_diag)
              palb_d   = ssa_bin(icol, ilev, idx_sw_diag)

              ! Accumulate per-species AODs
              aodabsbc(icol) = aodabsbc(icol) + absbc_d(icol)*dopaer_d*(1.0_kind_phys - palb_d)

              aodc_d = (abssulf_d(icol)*(1.0_kind_phys - palb_d) + palb_d*scatsulf_d(icol))*dopaer_d
              sulfaod(icol) = sulfaod(icol) + aodc_d

              aodc_d = (abspom_d(icol)*(1.0_kind_phys - palb_d) + palb_d*scatpom_d(icol))*dopaer_d
              pomaod(icol) = pomaod(icol) + aodc_d

              aodc_d = (abssoa_d(icol)*(1.0_kind_phys - palb_d) + palb_d*scatsoa_d(icol))*dopaer_d
              soaaod(icol) = soaaod(icol) + aodc_d

              aodc_d = (absbc_d(icol)*(1.0_kind_phys - palb_d) + palb_d*scatbc_d(icol))*dopaer_d
              bcaod(icol) = bcaod(icol) + aodc_d

              aodc_d = (abssslt_d(icol)*(1.0_kind_phys - palb_d) + palb_d*scatsslt_d(icol))*dopaer_d
              ssltaod(icol) = ssltaod(icol) + aodc_d

              aodc_d = (absdust_d(icol)*(1.0_kind_phys - palb_d) + palb_d*scatdust_d(icol))*dopaer_d
              dustaod(icol) = dustaod(icol) + aodc_d
            end if
          end do ! icol
        end do ! ilev (per-species diagnostics)

        !-------------------------------------------------
        ! Per-bin column visible tau for bulk aerosol (ODV_ diagnostics)
        ! (was in aer_vis_diag_mod in CAM)
        !-------------------------------------------------
        if (aeroprops%model_is('BAM')) then
          bam_cnt = bam_cnt + 1
          if (bam_cnt <= num_bulk_aer) then
            ! TODO: Use hist_fld_active to gate this when available in CAM-SIMA
            do icol = 1, ncol
              odv_col_aod(icol, bam_cnt) = sum(tau_bin(icol, top_lev:pver, idx_sw_diag))
            end do
          end if
        end if

        !-------------------------------------------------
        ! Call the portable core for longwave optics calculation per-bin:
        !-------------------------------------------------
        call aerosol_optics_lw_bin(aeroprops, aerostate, ibin, &
                                   ncol, pver, nswbands, nlwbands, numrh, &
                                   relh, sulfwtpct, mass, crefwsw, crefwlw, &
                                   geometric_radius, &
                                   tau_lw_bin, absorp_bin, &
                                   errmsg, errflg)
        if (errflg /= 0) return

        ! Accumulate LW absorption optical depth
        do iwav = 1, nlwbands
          do ilev = 1, pver
            do icol = 1, ncol
              aer_lw_abs(icol, ilev, iwav) = aer_lw_abs(icol, ilev, iwav) &
                                             + tau_lw_bin(icol, ilev, iwav)
            end do
          end do
        end do

      end do bin_loop ! ibin
    end do aerosol_model_loop ! iaermod

    ! Compute column-level SSA and AOD at vis from total aer_tau/aer_tau_w
    do icol = 1, ncol
      aodvis(icol) = sum(aer_tau(icol, :, idx_sw_diag))
      if (aodvis(icol) > 1.e-10_kind_phys) then
        ssavis(icol) = sum(aer_tau_w(icol, :, idx_sw_diag)) / aodvis(icol)
      else
        ssavis(icol) = 0.925_kind_phys
      end if
    end do

    ! Output diagnostics for active climate and diagnostic lists
    do ilist = 0, N_DIAG
      if (active_calls(ilist)) then
        call rad_aer_diag_out(ilist, constituents, pdeldry, ncol)
      end if
    end do

  end subroutine aerosol_optics_run

end module aerosol_optics
