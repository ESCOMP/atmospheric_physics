module rrtmgp_inputs_setup

 implicit none
 private

 public :: rrtmgp_inputs_setup_init

 contains
!> \section arg_table_rrtmgp_inputs_setup_init Argument Table
!! \htmlinclude rrtmgp_inputs_setup_init.html
!!
  subroutine rrtmgp_inputs_setup_init(nswbands, nlwbands, pref_edge, pver, pverp, kdist_sw, &
                   kdist_lw, qrl, is_first_step, use_rad_dt_cosz, timestep_size, nstep, iradsw, dt_avg,  &
                   irad_always, is_first_restart_step, p_top_for_rrtmgp, nradgas, gasnamelength, current_cal_day,  &
                   ktopcam, ktoprad, nlaycam, sw_low_bounds, sw_high_bounds, idx_sw_diag, idx_nir_diag,        &
                   idx_uv_diag, idx_sw_cloudsim, idx_lw_diag, idx_lw_cloudsim, nswgpts, nlwgpts, changeseed,   &
                   nlay, nlayp, nextsw_cday, band2gpt_sw, irad_always_out, errmsg, errflg)
     use ccpp_kinds,             only: kind_phys
     use ccpp_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp_ccpp
     use radiation_utils,        only: radiation_utils_init, get_sw_spectral_boundaries_ccpp

     ! Inputs
     integer,                         intent(in) :: nswbands
     integer,                         intent(in) :: nlwbands               ! Number of longwave bands
     integer,                         intent(in) :: nradgas                ! Number of radiatively active gases
     integer,                         intent(in) :: pverp                  ! Number of vertical interfaces
     integer,                         intent(in) :: pver                   ! Number of vertical layers
     integer,                         intent(in) :: iradsw                 ! Freq. of shortwave radiation calc in time steps (positive) or hours (negative).
     integer,                         intent(in) :: nstep                  ! Current timestep number
     integer,                         intent(in) :: gasnamelength          ! Length of all of the gas_list entries
     integer,                         intent(in) :: irad_always            ! Number of time steps to execute radiation continuously
     real(kind_phys),                 intent(in) :: timestep_size          ! Timestep size (s)
     real(kind_phys),                 intent(in) :: current_cal_day        ! Current calendar day
     real(kind_phys), dimension(:),   intent(in) :: pref_edge              ! Reference pressures (interfaces) (Pa)
     type(ty_gas_optics_rrtmgp_ccpp), intent(in) :: kdist_sw               ! Shortwave gas optics object
     type(ty_gas_optics_rrtmgp_ccpp), intent(in) :: kdist_lw               ! Longwave gas optics object
     logical,                         intent(in) :: is_first_step          ! Flag for whether this is the first timestep (.true. = yes)
     logical,                         intent(in) :: is_first_restart_step  ! Flag for whether this is the first restart step (.true. = yes)
     logical,                         intent(in) :: use_rad_dt_cosz        ! Use adjusted radiation timestep for cosz calculation
     real(kind_phys),                 intent(in) :: p_top_for_rrtmgp       ! Top pressure to use for RRTMGP (Pa)

     ! Outputs
     integer,                         intent(out) :: ktopcam               ! Index in CAM arrays of top level (layer or interface) at which RRTMGP is active
     integer,                         intent(out) :: ktoprad               ! Index in RRTMGP array corresponding to top layer or interface of CAM arrays
     integer,                         intent(out) :: nlaycam               ! Number of vertical layers in CAM. Is either equal to nlay
                                                                           !  or is 1 less than nlay if "extra layer" is used in the radiation calculations
     integer,                         intent(out) :: nlay                  ! Number of vertical layers in radiation calculation
     integer,                         intent(out) :: nlayp                 ! Number of vertical interfaces in radiation calculations (nlay + 1)
     ! Indices to specific bands for diagnostic output and COSP input
     integer,                         intent(out) :: idx_sw_diag           ! Index of band containing 500-nm wave
     integer,                         intent(out) :: idx_nir_diag          ! Index of band containing 1000-nm wave
     integer,                         intent(out) :: idx_uv_diag           ! Index of band containing 400-nm wave
     integer,                         intent(out) :: idx_sw_cloudsim       ! Index of band for shortwave cosp diagnostics
     integer,                         intent(out) :: idx_lw_diag           ! Index of band containing 1000-cm-1 wave (H2O window)
     integer,                         intent(out) :: idx_lw_cloudsim       ! Index of band for longwave cosp diagnostics

     integer,                         intent(out) :: nswgpts               ! Number of shortwave g-points
     integer,                         intent(out) :: nlwgpts               ! Number of longwave g-points
     integer,                         intent(out) :: changeseed            ! Random number seed for mcica longwave
     integer, dimension(:,:),         intent(out) :: band2gpt_sw           ! Array for converting shortwave band limits to g-points
     real(kind_phys),                 intent(out) :: nextsw_cday           ! The next calendar day during which the shortwave radiation calculation will be performed
     real(kind_phys), dimension(:),   intent(out) :: sw_low_bounds         ! Lower bounds of shortwave bands
     real(kind_phys), dimension(:),   intent(out) :: sw_high_bounds        ! Upper bounds of shortwave bands
     real(kind_phys), dimension(:,:), intent(inout) :: qrl                 ! Longwave radiative heating
     character(len=*),                intent(out) :: errmsg
     integer,                         intent(out) :: errflg
     integer,                         intent(out) :: irad_always_out       ! Number of time steps to execute radiation continuously
     real(kind_phys),                 intent(out) :: dt_avg                ! averaging time interval for zenith angle

     ! Local variables
     real(kind_phys), target :: wavenumber_low_shortwave(nswbands)
     real(kind_phys), target :: wavenumber_high_shortwave(nswbands)
     real(kind_phys), target :: wavenumber_low_longwave(nlwbands)
     real(kind_phys), target :: wavenumber_high_longwave(nlwbands)

     ! Set error variables
     errflg = 0
     errmsg = ''

     ! Number of layers in radiation calculation is capped by the number of
     ! pressure interfaces below p_top_for_rrtmgp.  When the entire model atmosphere is
     ! below p_top_for_rrtmgp then an extra layer is added to the top of the model for
     ! the purpose of the radiation calculation.
     nlay = count( pref_edge(:) > p_top_for_rrtmgp )
     nlayp = nlay + 1

     if (nlay == pverp) then
        ! Top model interface is below p_top_for_rrtmgp.  RRTMGP is active in all model
        ! layers plus 1 extra layer between model top and p_top_for_rrtmgp.
        ktopcam = 1
        ktoprad = 2
        nlaycam = pver
     else
        ! nlay < pverp.  nlay layers are used in radiation calcs, and they are
        ! all CAM layers.
        ktopcam = pver - nlay + 1
        ktoprad = 1
        nlaycam = nlay
     end if

     ! Set the sw/lw band boundaries in radconstants.  Also sets
     ! indicies of specific bands for diagnostic output and COSP input.
     call set_wavenumber_bands(kdist_sw, kdist_lw, nswbands, nlwbands, idx_sw_diag, idx_nir_diag, &
                 idx_uv_diag, idx_sw_cloudsim, idx_lw_diag, idx_lw_cloudsim, nswgpts, nlwgpts,    &
                 wavenumber_low_shortwave, wavenumber_high_shortwave, wavenumber_low_longwave,    &
                 wavenumber_high_longwave, band2gpt_sw, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     call radiation_utils_init(nswbands, nlwbands, wavenumber_low_shortwave, wavenumber_high_shortwave, &
             wavenumber_low_longwave, wavenumber_high_longwave, errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     ! Initialize the SW band boundaries
     call get_sw_spectral_boundaries_ccpp(sw_low_bounds, sw_high_bounds, 'cm-1', errmsg, errflg)
     if (errflg /= 0) then
        return
     end if

     if (is_first_step) then
        qrl = 0._kind_phys
     end if

     ! Set the radiation timestep for cosz calculations if requested using
     ! the adjusted iradsw value from radiation
     if (use_rad_dt_cosz)  then
        dt_avg = iradsw*timestep_size
     else
        dt_avg = 0._kind_phys
     end if

     ! "irad_always" is number of time steps to execute radiation continuously from
     ! start of initial OR restart run
     irad_always_out = irad_always
     if (irad_always > 0) then
        irad_always_out = irad_always + nstep
     end if

     ! Surface components to get radiation computed today
     if (.not. is_first_restart_step) then
        nextsw_cday = current_cal_day
     end if

     changeseed = nlwgpts

  end subroutine rrtmgp_inputs_setup_init

!=========================================================================================
!                                     HELPER FUNCTIONS                                   !
!=========================================================================================
  subroutine set_wavenumber_bands(kdist_sw, kdist_lw, nswbands, nlwbands, idx_sw_diag, idx_nir_diag, &
                  idx_uv_diag, idx_sw_cloudsim, idx_lw_diag, idx_lw_cloudsim, nswgpts, nlwgpts,      &
                  wavenumber_low_shortwave, wavenumber_high_shortwave, wavenumber_low_longwave,      &
                  wavenumber_high_longwave, band2gpt_sw, errmsg, errflg)
   use ccpp_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_kinds,             only: kind_phys
   ! Set the low and high limits of the wavenumber grid for sw and lw.
   ! Values come from RRTMGP coefficients datasets, and are stored in the
   ! kdist objects.
   !
   ! Set band indices for bands containing specific wavelengths.

   ! Arguments
   type(ty_gas_optics_rrtmgp_ccpp), intent(in) :: kdist_sw
   type(ty_gas_optics_rrtmgp_ccpp), intent(in) :: kdist_lw
   integer,                         intent(in) :: nswbands
   integer,                         intent(in) :: nlwbands

   integer,                        intent(out) :: idx_sw_diag
   integer,                        intent(out) :: idx_nir_diag
   integer,                        intent(out) :: idx_uv_diag
   integer,                        intent(out) :: idx_sw_cloudsim
   integer,                        intent(out) :: idx_lw_diag
   integer,                        intent(out) :: idx_lw_cloudsim
   integer,                        intent(out) :: nswgpts
   integer,                        intent(out) :: nlwgpts
   integer, dimension(:,:),        intent(out) :: band2gpt_sw
   real(kind_phys), dimension(:),  intent(out) :: wavenumber_low_shortwave
   real(kind_phys), dimension(:),  intent(out) :: wavenumber_high_shortwave
   real(kind_phys), dimension(:),  intent(out) :: wavenumber_low_longwave
   real(kind_phys), dimension(:),  intent(out) :: wavenumber_high_longwave
   character(len=*),               intent(out) :: errmsg
   integer,                        intent(out) :: errflg

   ! Local variables
   integer :: istat
   real(kind_phys), allocatable :: values(:,:)
   character(len=256) :: alloc_errmsg

   character(len=*), parameter :: sub = 'set_wavenumber_bands'
   !----------------------------------------------------------------------------

   ! Initialize error variables
   errflg = 0
   errmsg = ''
   ! Check that number of sw/lw bands in gas optics files matches the parameters.
   if (kdist_sw%gas_props%get_nband() /= nswbands) then
      write(errmsg,'(a, a,i4,a,i4)') sub, ': ERROR: number of sw bands in file, ', kdist_sw%gas_props%get_nband(), &
         ", doesn't match parameter nswbands= ", nswbands
      errflg = 1
      return
   end if
   if (kdist_lw%gas_props%get_nband() /= nlwbands) then
      write(errmsg,'(a, a,i4,a,i4)') sub, ': ERROR: number of lw bands in file, ', kdist_lw%gas_props%get_nband(), &
         ", doesn't match parameter nlwbands= ", nlwbands
      errflg = 1
      return
   end if

   nswgpts = kdist_sw%gas_props%get_ngpt()
   nlwgpts = kdist_lw%gas_props%get_ngpt()

   ! SW band bounds in cm^-1
   allocate( values(2,nswbands), stat=istat, errmsg=alloc_errmsg )
   if (istat/=0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR allocating array: values(2,nswbands); message - ', alloc_errmsg
      errflg = 1
      return
   end if
   values = kdist_sw%gas_props%get_band_lims_wavenumber()
   wavenumber_low_shortwave = values(1,:)
   wavenumber_high_shortwave = values(2,:)

   ! First and last g-point for each SW band:
   band2gpt_sw = kdist_sw%gas_props%get_band_lims_gpoint()

   ! Indices into specific bands
   call get_band_index_by_value('sw', 500.0_kind_phys, 'nm', nswbands, &
        wavenumber_low_shortwave, wavenumber_high_shortwave, idx_sw_diag, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call get_band_index_by_value('sw', 1000.0_kind_phys, 'nm', nswbands, &
        wavenumber_low_shortwave, wavenumber_high_shortwave, idx_nir_diag, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call get_band_index_by_value('sw', 400._kind_phys, 'nm', nswbands, &
        wavenumber_low_shortwave, wavenumber_high_shortwave, idx_uv_diag, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call get_band_index_by_value('sw', 0.67_kind_phys, 'micron', nswbands, &
        wavenumber_low_shortwave, wavenumber_high_shortwave, idx_sw_cloudsim, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if

   deallocate(values)

   ! LW band bounds in cm^-1
   allocate( values(2,nlwbands), stat=istat, errmsg=alloc_errmsg )
   if (istat/=0) then
      write(errmsg, '(a,a,a)') sub, ': ERROR allocating array: values(2,nlwbands); message - ', alloc_errmsg
      errflg = 1
      return
   end if
   values = kdist_lw%gas_props%get_band_lims_wavenumber()
   wavenumber_low_longwave = values(1,:)
   wavenumber_high_longwave = values(2,:)

   ! Indices into specific bands
   call get_band_index_by_value('lw', 1000.0_kind_phys, 'cm^-1', nlwbands, &
        wavenumber_low_longwave, wavenumber_high_longwave, idx_lw_diag, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call get_band_index_by_value('lw', 10.5_kind_phys, 'micron', nlwbands, &
        wavenumber_low_longwave, wavenumber_high_longwave, idx_lw_cloudsim, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if

 end subroutine set_wavenumber_bands

!=========================================================================================

 subroutine get_band_index_by_value(swlw, targetvalue, units, nbnds, wavenumber_low, &
                wavenumber_high, ans, errmsg, errflg)
   use ccpp_kinds, only: kind_phys

   ! Find band index for requested wavelength/wavenumber.

   character(len=*),                      intent(in) :: swlw        ! sw or lw bands
   real(kind_phys),                       intent(in) :: targetvalue  
   character(len=*),                      intent(in) :: units       ! units of targetvalue
   integer,                               intent(in) :: nbnds
   real(kind_phys),  dimension(:),        intent(in) :: wavenumber_low
   real(kind_phys),  dimension(:),        intent(in) :: wavenumber_high
   character(len=*),                     intent(out) :: errmsg
   integer,                              intent(out) :: errflg
   integer,                              intent(out) :: ans

   ! local
   real(kind_phys) :: tgt
   integer  :: idx

   character(len=*), parameter :: sub = 'get_band_index_by_value'
   !----------------------------------------------------------------------------

   ! Initialize error variables
   errflg = 0
   errmsg = ''
   if (trim(swlw) /= 'sw' .and. trim(swlw) /= 'lw') then
      write(errmsg,'(a,a)') 'rrtmgp_inputs_setup: get_band_index_by_value: type of bands not recognized: ', swlw
      errflg = 1
      return
   end if

   ! band info is in cm^-1 but target value may be other units,
   ! so convert targetvalue to cm^-1
   select case (units)
   case ('inv_cm','cm^-1','cm-1')
      tgt = targetvalue
   case('m','meter','meters')
      tgt = 1.0_kind_phys / (targetvalue * 1.e2_kind_phys)
   case('nm','nanometer','nanometers')
      tgt = 1.0_kind_phys / (targetvalue * 1.e-7_kind_phys)
   case('um','micrometer','micrometers','micron','microns')
      tgt = 1.0_kind_phys / (targetvalue * 1.e-4_kind_phys)
   case('cm','centimeter','centimeters')
      tgt = 1._kind_phys/targetvalue
   case default
      write(errmsg,'(a,a)') 'rrtmgp_inputs_setup: get_band_index_by_value: units not recognized: ', units
      errflg = 1
   end select

   ! now just loop through the array
   ans = 0
   do idx = 1,nbnds
      if ((tgt > wavenumber_low(idx)) .and. (tgt <= wavenumber_high(idx))) then
         ans = idx
         exit
      end if
   end do

   if (ans == 0) then
      write(errmsg,'(a,f10.3,a,a)') 'rrtmgp_inputs_setup: get_band_index_by_value: no index found for wavelength ', targetvalue, ' ', trim(units)
      errflg = 1
   end if
   
 end subroutine get_band_index_by_value

end module rrtmgp_inputs_setup
