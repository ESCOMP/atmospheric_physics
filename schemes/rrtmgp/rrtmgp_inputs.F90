module rrtmgp_inputs

 implicit none
 private

 public :: rrtmgp_inputs_run

 contains

!> \section arg_table_rrtmgp_inputs_run Argument Table
!! \htmlinclude rrtmgp_inputs_run.html
!!
  subroutine rrtmgp_inputs_run(dosw, dolw, snow_associated, graupel_associated, &
                  pmid, pint, t, nday, idxday, cldfprime, &
                  coszrs, kdist_sw, t_sfc, emis_sfc, t_rad, pmid_rad,     &
                  pint_rad, t_day, pmid_day, pint_day, coszrs_day,        &
                  alb_dir, alb_dif, lwup, stebol, ncol, ktopcam, ktoprad, &
                  nswbands, asdir, asdif, sw_low_bounds, sw_high_bounds,  &
                  aldir, aldif, nlay, pverp, pver, cld, cldfsnow,         &
                  cldfgrau, graupel_in_rad, gasnamelength, gaslist,       &
                  gas_concs_lw, aer_lw, atm_optics_lw, kdist_lw,          &
                  sources_lw, aer_sw, atm_optics_sw, gas_concs_sw,        &
                  errmsg, errflg)
     use ccpp_kinds,              only: kind_phys
     use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
     use ccpp_optical_props,      only: ty_optical_props_1scl_ccpp, ty_optical_props_2str_ccpp
     use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
     use ccpp_source_functions,   only: ty_source_func_lw_ccpp
     use atmos_phys_string_utils, only: to_lower
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
     character(len=*), dimension(:),       intent(in) :: gaslist               ! Radiatively active gases
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
     character(len=*),                     intent(out) :: errmsg
     integer,                              intent(out) :: errflg

     ! Local variables
     real(kind_phys) :: tref_min
     real(kind_phys) :: tref_max
     integer :: idx, kdx, iband
     character(len=gasnamelength) :: gaslist_lc(size(gaslist))

     ! Set error variables
     errmsg = ''
     errflg = 0

     if (.not. dosw .and. .not. dolw) then
        return
     end if

     ! RRTMGP set state
     t_sfc = sqrt(sqrt(lwup(:)/stebol))  ! Surface temp set based on longwave up flux.

     ! Set surface emissivity to 1.0.
     ! The land model *does* have its own surface emissivity, but is not spectrally resolved.
     ! The LW upward flux is calculated with that land emissivity, and the "radiative temperature"
     ! t_sfc is derived from that flux. We assume, therefore, that the emissivity is unity
     ! to be consistent with t_sfc.
     emis_sfc(:,:) = 1._kind_phys

     ! Level ordering is the same for both CAM and RRTMGP (top to bottom)
     t_rad(:,ktoprad:) = t(:,ktopcam:)
     pmid_rad(:,ktoprad:) = pmid(:,ktopcam:)
     pint_rad(:,ktoprad:) = pint(:,ktopcam:)

     ! Add extra layer values if needed.
     if (nlay == pverp) then
        t_rad(:,1) = t(:,1)
        ! The top reference pressure from the RRTMGP coefficients datasets is 1.005183574463 Pa
        ! Set the top of the extra layer just below that.
        pint_rad(:,1) = 1.01_kind_phys

        ! next interface down in LT will always be > 1Pa
        ! but in MT we apply adjustment to have it be 1.02 Pa if it was too high
        where (pint_rad(:,2) <= pint_rad(:,1)) pint_rad(:,2) = pint_rad(:,1)+0.01_kind_phys

        ! set the highest pmid (in the "extra layer") to the midpoint (guarantees > 1Pa)
        pmid_rad(:,1)   = pint_rad(:,1) + 0.5_kind_phys * (pint_rad(:,2) - pint_rad(:,1))

        ! For case of CAM MT, also ensure pint_rad(:,2) > pint_rad(:,1) & pmid_rad(:,2) > max(pmid_rad(:,1), min_pressure)
        where (pmid_rad(:,2) <= kdist_sw%gas_props%get_press_min()) pmid_rad(:,2) = pint_rad(:,2) + 0.01_kind_phys
     else
        ! nlay < pverp, thus the 1 Pa level is within a CAM layer.  Assuming the top interface of
        ! this layer is at a pressure < 1 Pa, we need to adjust the top of this layer so that it
        ! is within the valid pressure range of RRTMGP (otherwise RRTMGP issues an error).  Then
        ! set the midpoint pressure halfway between the interfaces.
        pint_rad(:,1) = 1.01_kind_phys
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

     ! Create lowercase version of the gaslist for RRTMGP.  The ty_gas_concs_ccpp objects
     ! work with CAM's uppercase names, but other objects that get input from the gas
     ! concs objects don't work.
     do idx = 1, size(gaslist)
        gaslist_lc(idx) = to_lower(gaslist(idx))
     end do

     ! If no daylight columns, can't create empty RRTMGP objects
     if (dosw .and. nday > 0) then
        ! Initialize object for gas concentrations.
         errmsg = gas_concs_sw%gas_concs%init(gaslist_lc)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

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

        ! Initialize object for Planck sources.
        errmsg = sources_lw%sources%alloc(ncol, nlay, kdist_lw%gas_props)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
     end if

  end subroutine rrtmgp_inputs_run

end module rrtmgp_inputs
