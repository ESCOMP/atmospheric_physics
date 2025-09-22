module rrtmgp_inputs
 use cam_logfile, only: iulog

 implicit none
 private

 public :: rrtmgp_inputs_run

 contains

!> \section arg_table_rrtmgp_inputs_run Argument Table
!! \htmlinclude rrtmgp_inputs_run.html
!!
  subroutine rrtmgp_inputs_run(dosw, dolw, snow_associated, graupel_associated, &
                  is_root, iulog, is_mpas, pmid, pint, t, nday, idxday,          &
                  cldfprime, coszrs, kdist_sw, t_sfc, emis_sfc, t_rad,          &
                  pmid_rad, pint_rad, t_day, pmid_day, pint_day, coszrs_day,    &
                  alb_dir, alb_dif, lwup, stebol, ncol, ktopcam, ktoprad, &
                  nswbands, asdir, asdif, sw_low_bounds, sw_high_bounds,  &
                  aldir, aldif, nlay, pverp, pver, cld, cldfsnow,         &
                  cldfgrau, graupel_in_rad, gasnamelength, gaslist_lc,    &
                  gas_concs_lw, aer_lw, atm_optics_lw, kdist_lw,          &
                  sources_lw, aer_sw, atm_optics_sw, gas_concs_sw,        &
                  errmsg, errflg)
     use ccpp_kinds,              only: kind_phys
     use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
     use ccpp_optical_props,      only: ty_optical_props_1scl_ccpp, ty_optical_props_2str_ccpp
     use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
     use ccpp_source_functions,   only: ty_source_func_lw_ccpp
     use atmos_phys_rad_utils,    only: is_visible

     ! Inputs
     logical,                              intent(in) :: graupel_in_rad        ! Flag to include graupel in radiation calculation
     integer,                              intent(in) :: ncol                  ! Number of columns
     integer,                              intent(in) :: pver                  ! Number of vertical layers
     integer,                              intent(in) :: pverp                 ! Number of vertical interfaces
     integer,                              intent(in) :: nlay                  ! Number of vertical layers used in radiation calculation
     integer,                              intent(in) :: nswbands              ! Number of shortwave bands
     integer,                              intent(in) :: ktopcam               ! Index in CAM arrays of top level (layer or interface) at which RRTMGP is active
     integer,                              intent(in) :: ktoprad               ! Index in RRTMGP array corresponding to top layer or interface of CAM arrays
     integer,                              intent(in) :: gasnamelength         ! Length of gases in gas_list
     integer,                              intent(in) :: nday                  ! Number of daylight columns
     logical,                              intent(in) :: dosw                  ! Flag for performing the shortwave calculation
     logical,                              intent(in) :: dolw                  ! Flag for performing the longwave calculation
     logical,                              intent(in) :: snow_associated       ! Flag for whether the cloud snow fraction argument should be used
     logical,                              intent(in) :: graupel_associated    ! Flag for whether the cloud graupel fraction argument should be used
     logical,                              intent(in) :: is_root
     logical,                              intent(in) :: is_mpas
     integer,                              intent(in) :: iulog
     integer,         dimension(:),        intent(in) :: idxday                ! Indices of daylight columns
     real(kind_phys), dimension(:,:),      intent(in) :: pmid                  ! Air pressure at midpoint (Pa)
     real(kind_phys), dimension(:,:),      intent(in) :: pint                  ! Air pressure at interface (Pa)
     real(kind_phys), dimension(:,:),      intent(in) :: t                     ! Air temperature (K)
     real(kind_phys), dimension(:,:),      intent(in) :: cldfsnow              ! Cloud fraction of just "snow clouds"
     real(kind_phys), dimension(:,:),      intent(in) :: cldfgrau              ! Cloud fraction of just "graupel clouds"
     real(kind_phys), dimension(:,:),      intent(in) :: cld                   ! Cloud fraction (liq+ice)
     real(kind_phys), dimension(:),        intent(in) :: sw_low_bounds         ! Lower bounds for shortwave bands
     real(kind_phys), dimension(:),        intent(in) :: sw_high_bounds        ! Upper bounds for shortwave bands
     real(kind_phys), dimension(:),        intent(in) :: coszrs                ! Cosine of solar senith angle (radians)
     real(kind_phys), dimension(:),        intent(in) :: lwup                  ! Longwave up flux (W m-2)
     real(kind_phys), dimension(:),        intent(in) :: asdir                 ! Shortwave direct albedo (fraction)
     real(kind_phys), dimension(:),        intent(in) :: asdif                 ! Shortwave diffuse albedo (fraction)
     real(kind_phys), dimension(:),        intent(in) :: aldir                 ! Longwave direct albedo (fraction)
     real(kind_phys), dimension(:),        intent(in) :: aldif                 ! Longwave diffuse albedo (fraction)
     real(kind_phys),                      intent(in) :: stebol                ! Stefan-Boltzmann constant (W m-2 K-4)
     type(ty_gas_optics_rrtmgp_ccpp),      intent(in) :: kdist_sw              ! Shortwave gas optics object
     type(ty_gas_optics_rrtmgp_ccpp),      intent(in) :: kdist_lw              ! Longwave gas optics object
     character(len=5), dimension(:),       intent(in) :: gaslist_lc            ! Radiatively active gases
     ! Outputs
     real(kind_phys), dimension(:,:),      intent(out) :: t_rad                ! Air temperature with radiation indexing (K)
     real(kind_phys), dimension(:,:),      intent(out) :: pmid_rad             ! Midpoint pressure with radiation indexing (Pa)
     real(kind_phys), dimension(:,:),      intent(out) :: pint_rad             ! Interface pressure with radiation indexing (Pa)
     real(kind_phys), dimension(:,:),      intent(out) :: t_day                ! Air temperature of daylight columns (K)
     real(kind_phys), dimension(:,:),      intent(out) :: pint_day             ! Interface pressure of daylight columns (Pa)
     real(kind_phys), dimension(:,:),      intent(out) :: pmid_day             ! Midpoint pressure of daylight columns (Pa)
     real(kind_phys), dimension(:,:),      intent(out) :: emis_sfc             ! Surface emissivity (fraction)
     real(kind_phys), dimension(:,:),      intent(out) :: alb_dir              ! Surface albedo due to UV and VIS direct (fraction)
     real(kind_phys), dimension(:,:),      intent(out) :: alb_dif              ! Surface albedo due to IR diffused (fraction)
     real(kind_phys), dimension(:,:),      intent(out) :: cldfprime            ! modified cloud fraciton

     real(kind_phys), dimension(:),        intent(out) :: t_sfc                ! Surface temperature (K)
     real(kind_phys), dimension(:),        intent(out) :: coszrs_day           ! Cosine of solar zenith angle for daylight columns (radians)
     type(ty_gas_concs_ccpp),              intent(out) :: gas_concs_lw         ! Gas concentrations object for longwave radiation
     type(ty_optical_props_1scl_ccpp),     intent(out) :: atm_optics_lw        ! Atmosphere optical properties object for longwave radiation
     type(ty_optical_props_1scl_ccpp),     intent(out) :: aer_lw               ! Aerosol optical properties object for longwave radiation
     type(ty_source_func_lw_ccpp),         intent(out) :: sources_lw           ! Longwave sources object for longwave radiation
     type(ty_gas_concs_ccpp),              intent(out) :: gas_concs_sw         ! Gas concentrations object for shortwave radiation
     type(ty_optical_props_2str_ccpp),     intent(out) :: atm_optics_sw        ! Atmosphere optical properties object for shortwave radiation
     type(ty_optical_props_2str_ccpp),     intent(out) :: aer_sw               ! Aerosol optical properties object for shortwave radiation
     character(len=512),                   intent(out) :: errmsg
     integer,                              intent(out) :: errflg

     ! Local variables
     real(kind_phys) :: tref_min
     real(kind_phys) :: tref_max
     integer :: idx, kdx, iband
     logical :: ltrick_rrtmgp

     ! Set error variables
     errmsg = ''
     errflg = 0

     if (.not. dosw .and. .not. dolw) then
        return
     end if

     !------------------------------------------------------------------------------
     ! Compile logic to determine whether it is necessary AND possible to 'trick'
     ! rrtmgp to violate its own validity limits by hacking its vertical grid.
     ! Conditions:
     !   1) top CAM interface (k=1) is above 1 Pa
     !   2) next interface (k=2) is below 1 Pa
     !   3) RRTMGP is asked to active over ALL CAM layers
     !      (nlay=pverp, e.g., set by spec p_top_for_equil_rad=0.)
     !   4) dycore is NOT MPAS
     !
     ! These conditions are generally only satisfied in a non-MPAS MT configuration
     !------------------------------------------------------------------------------
     if (( .not. is_mpas ) .and. &
          (nlay==pverp) .and. &
          (minval(pint(:,1)) < 1._kind_phys) .and. &
          (minval(pint(:,2)) > 1._kind_phys) ) then
        ! we can and need to trick rrtmgp
        ltrick_rrtmgp = .true.
     else
        ! we cannot or don't need to trick rrtmgp
        ltrick_rrtmgp = .false.
     end if

     if (is_root) then
        if (ltrick_rrtmgp) then
           write(iulog,*) ' *** TRICKING RRTMGP INTO GOING AN EXTRA LEVEL  ',nlay,pverp
        else
           write(iulog,*) ' *** CANT or WONT trick RRTMGP ',nlay,pverp
        end if
     end if

     ! RRTMGP set state
     t_sfc = sqrt(sqrt(lwup(:)/stebol))  ! Surface temp set based on longwave up flux.

     ! Set surface emissivity to 1.0.
     ! The land model *does* have its own surface emissivity, but is not spectrally resolved.
     ! The LW upward flux is calculated with that land emissivity, and the "radiative temperature"
     ! t_sfc is derived from that flux. We assume, therefore, that the emissivity is unity
     ! to be consistent with t_sfc.
     emis_sfc(:,:) = 1._kind_phys

     !-------------------------------------------------------------------------
     ! RRTMGP enforces  P > 1 Pa for validity.
     ! Actual range of RRTMGP in CAM is set with namelist variable p_top_for_equil_rad.
     ! In rrtmg_inputs_setup.F90, active layers for RRTMGP are counted based on
     ! the logical P_ref > p_top_for_equil_rad. Returned as nlay.
     !
     ! If 
     ! 1) entire vertical domain has P_ref> p_top_for_equil_rad then
     !    nlay = pverp
     !    ktoprad = 2
     !    ktopcam = 1
     ! 2) min(P_ref) < p_top_for_equil_rad (e.g. MPAS MT) then
     !    nlay < pverp
     !    ktoprad = 1
     !    ktopcam = pver - nlay + 1
     !
     ! Level ordering is the same for both CAM and RRTMGP (top to bottom)
     ! Note in Case 2, tops of {t,pint,pmid}_rad start at a 'valid' level,
     ! i.e. pref > p_top_for_rrtmgp
     !-------------------------------------------------------
     t_rad(:,ktoprad:) = t(:,ktopcam:)
     pmid_rad(:,ktoprad:) = pmid(:,ktopcam:)
     pint_rad(:,ktoprad:) = pint(:,ktopcam:)

     ! Deal with vertical grid for RRTMGP 
     if (nlay == pverp) then
        ! All model layers are within RRTMGP's range of valid pressures
        ! as specified by p_top_for_equil_rad - (Case 1 above)
        t_rad(:,1)      = t(:,1)
        if (ltrick_rrtmgp) then
           ! The top reference pressure from the RRTMGP coefficients datasets is 1.005183574463 Pa
           ! Set the top of the extra layer just below that.
           pint_rad(:,1) = 1.01_kind_phys
           pint_rad(:,2) = 1.02_kind_phys
           ! set the highest pmid (in the "extra layer") to the midpoint (guarantees > 1Pa)
           pmid_rad(:,1)   = 1.015_kind_phys ! pint_rad(:,1) + 0.5_kind_phys * (pint_rad(:,2) - pint_rad(:,1))
           pmid_rad(:,2)   = 0.5*(pint_rad(:,2) + pint_rad(:,3))
        else
           ! The top reference pressure from the RRTMGP coefficients datasets is 1.005183574463 Pa
           ! Set the top of the extra layer just below that.
           pint_rad(:,1) = 1.01_kind_phys
           ! set the highest pmid (in the "extra layer") to the midpoint (guarantees > 1Pa)
           pmid_rad(:,1)   = pint_rad(:,1) + 0.5_kind_phys * (pint_rad(:,2) - pint_rad(:,1))
        end if
     else
        ! nlay < pverp : model min(pref) < p_top_for_rrtmgp  (Case 2 above)
        ! Not sure why is this needed since pint_rad should have been specified
        ! at RRTMGP valid values above. But this is the way it was done in
        ! original RRTMG implementation.
        pint_rad(:,1) = 1.01_kind_phys
        ! The following *should* work since pint_rad is all in valid range.
        pmid_rad(:,1) = 0.5_kind_phys * (pint_rad(:,1) + pint_rad(:,2))
     end if

     ! Limit temperatures to be within the limits of RRTMGP validity.
     tref_min = kdist_sw%gas_props%get_temp_min()
     tref_max = kdist_sw%gas_props%get_temp_max()
     t_rad = merge(t_rad, tref_min, t_rad > tref_min)
     t_rad = merge(t_rad, tref_max, t_rad < tref_max)

     ! Construct arrays containing only daylight columns
     do idx = 1, nday
        t_day(idx,:)    = t_rad(idxday(idx),:)
        pmid_day(idx,:) = pmid_rad(idxday(idx),:)
        pint_day(idx,:) = pint_rad(idxday(idx),:)
        coszrs_day(idx) = coszrs(idxday(idx))
     end do
     ! Assign albedos to the daylight columns (from E3SM implementation)
     ! Albedos are imported from the surface models as broadband (visible, and near-IR),
     ! and we need to map these to appropriate narrower bands used in RRTMGP. Bands
     ! are categorized broadly as "visible/UV" or "infrared" based on wavenumber.
     ! Loop over bands, and determine for each band whether it is broadly in the
     ! visible or infrared part of the spectrum based on a dividing line of
     ! 0.7 micron, or 14286 cm^-1
     do iband = 1,nswbands
        if (is_visible(sw_low_bounds(iband)) .and. &
           is_visible(sw_high_bounds(iband))) then

           ! Entire band is in the visible
           do idx = 1, nday
              alb_dir(iband,idx) = asdir(idxday(idx))
              alb_dif(iband,idx) = asdif(idxday(idx))
           end do

        else if (.not.is_visible(sw_low_bounds(iband)) .and. &
                 .not.is_visible(sw_high_bounds(iband))) then
           ! Entire band is in the longwave (near-infrared)
           do idx = 1, nday
              alb_dir(iband,idx) = aldir(idxday(idx))
              alb_dif(iband,idx) = aldif(idxday(idx))
           end do
        else
           ! Band straddles the visible to near-infrared transition, so we take
           ! the albedo to be the average of the visible and near-infrared
           ! broadband albedos
           do idx = 1, nday
              alb_dir(iband,idx) = 0.5_kind_phys * (aldir(idxday(idx)) + asdir(idxday(idx)))
              alb_dif(iband,idx) = 0.5_kind_phys * (aldif(idxday(idx)) + asdif(idxday(idx)))
           end do
        end if
     end do
     ! Strictly enforce albedo bounds
     where (alb_dir < 0)
        alb_dir = 0.0_kind_phys
     end where
     where (alb_dir > 1)
        alb_dir = 1.0_kind_phys
     end where
     where (alb_dif < 0)
        alb_dif = 0.0_kind_phys
     end where
     where (alb_dif > 1)
        alb_dif = 1.0_kind_phys
     end where

     ! modified cloud fraction
     ! Compute modified cloud fraction, cldfprime.
     ! 1. initialize as cld
     ! 2. modify for snow. use max(cld, cldfsnow)
     ! 3. modify for graupel if graupel_in_rad is true.
     !    use max(cldfprime, cldfgrau)
     if (snow_associated) then
        do kdx = 1, pver
           do idx = 1, ncol
              cldfprime(idx,kdx) = max(cld(idx,kdx), cldfsnow(idx,kdx))
           end do
        end do
     else
        cldfprime(:,:) = cld(:,:)
     end if

     if (graupel_associated .and. graupel_in_rad) then
        do kdx = 1, pver
           do idx = 1, ncol
              cldfprime(idx,kdx) = max(cldfprime(idx,kdx), cldfgrau(idx,kdx))
           end do
        end do
     end if

     ! If no daylight columns, can't create empty RRTMGP objects
     if (dosw .and. nday > 0) then
        ! Initialize object for gas concentrations.
         errmsg = gas_concs_sw%gas_concs%init(gaslist_lc)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        write(iulog,*) 'peverwhee - nlay before alloc'
        write(iulog,*) nlay

        ! Initialize object for combined gas + aerosol + cloud optics.
        ! Allocates arrays for properties represented on g-points.
        errmsg = atm_optics_sw%optical_props%alloc_2str(nday, nlay, kdist_sw%gas_props)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for SW aerosol optics.  Allocates arrays
        ! for properties represented by band.
        errmsg = aer_sw%optical_props%alloc_2str(nday, nlay, kdist_sw%gas_props%get_band_lims_wavenumber())
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
     end if

     if (dolw) then
        ! Initialize object for gas concentrations
        errmsg = gas_concs_lw%gas_concs%init(gaslist_lc)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for combined gas + aerosol + cloud optics.
        errmsg = atm_optics_lw%optical_props%alloc_1scl(ncol, nlay, kdist_lw%gas_props)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for LW aerosol optics.
        errmsg = aer_lw%optical_props%alloc_1scl(ncol, nlay, kdist_lw%gas_props%get_band_lims_wavenumber())
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
        ! PEVERWHEE - ZERO AEROSOLS FOR TESTING!
        aer_lw%optical_props%tau = 0.0_kind_phys

        ! Initialize object for Planck sources.
        errmsg = sources_lw%sources%alloc(ncol, nlay, kdist_lw%gas_props)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
     end if

  end subroutine rrtmgp_inputs_run

end module rrtmgp_inputs
