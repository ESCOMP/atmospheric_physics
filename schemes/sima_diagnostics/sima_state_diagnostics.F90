module sima_state_diagnostics

   use ccpp_kinds, only: kind_phys
   use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
   use cam_history_support,       only: fieldname_len

   implicit none
   private

   public :: sima_state_diagnostics_init ! init routine
   public :: sima_state_diagnostics_run  ! main routine

   character(len=*), parameter :: wv_std_name = 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water'

   character(len=65) :: const_std_names(6) = &
    [character(len=65) :: wv_std_name, &
     'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', &
     'rain_mixing_ratio_wrt_moist_air_and_condensed_water              ', &
     'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water         ', &
     'snow_mixing_ratio_wrt_moist_air_and_condensed_water              ', &
     'graupel_water_mixing_ratio_wrt_moist_air_and_condensed_water     ']

   character(len=6) :: const_diag_names(6) =  ['Q     ', &
                                               'CLDLIQ', &
                                               'RAINQM', &
                                               'CLDICE', &
                                               'SNOWQM', &
                                               'GRAUQM']

   ! Index of water vapor in the constituent array (0 if absent)
   ! QFLX and TMQ are only written when it is present.
   integer :: wv_idx = 0

CONTAINS

   !> \section arg_table_sima_state_diagnostics_init  Argument Table
   !! \htmlinclude sima_state_diagnostics_init.html
   subroutine sima_state_diagnostics_init(const_props, errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      integer :: const_idx, name_idx
      integer :: const_num_found
      character(len=64) :: units
      logical :: is_wet
      character(len=512) :: standard_name
      character(len=fieldname_len) :: diagnostic_name
      character(len=512) :: long_name
      character(len=3) :: mixing_ratio_type
      logical :: const_found(size(const_props))

      errmsg = ''
      errflg = 0

      ! Add state fields
      call history_add_field('PS',        'surface_pressure',                                              horiz_only,  'avg', 'Pa')
      call history_add_field('PSDRY',     'surface_pressure_of_dry_air',                                   horiz_only,  'avg', 'Pa')
      call history_add_field('PHIS',      'surface_geopotential',                                          horiz_only, 'inst', 'Pa')
      call history_add_field('T',         'air_temperature',                                                    'lev',  'avg', 'K')
      call history_add_field('U',         'eastward_wind',                                                      'lev',  'avg', 'm s-1')
      call history_add_field('V',         'northward_wind',                                                     'lev',  'avg', 'm s-1')
      call history_add_field('DSE',       'dry_static_energy',                                                  'lev',  'avg', 'm s-1')
      call history_add_field('OMEGA',     'lagrangian_tendency_of_air_pressure',                                'lev',  'avg', 'Pa s-1')
      call history_add_field('PMID',      'air_pressure',                                                       'lev',  'avg', 'Pa')
      call history_add_field('PMIDDRY',   'air_pressure_of_dry_air',                                            'lev',  'avg', 'Pa')
      call history_add_field('PDEL',      'air_pressure_thickness',                                             'lev',  'avg', 'Pa')
      call history_add_field('PDELDRY',   'air_pressure_thickness_of_dry_air',                                  'lev',  'avg', 'Pa')
      call history_add_field('RPDEL',     'reciprocal_of_air_pressure_thickness',                               'lev',  'avg', 'Pa-1')
      call history_add_field('RPDELDRY',  'reciprocal_of_air_pressure_thickness_of_dry_air',                    'lev',  'avg', 'Pa-1')
      call history_add_field('LNPMID',    'ln_air_pressure',                                                    'lev',  'avg', '1')
      call history_add_field('LNPMIDDRY', 'ln_air_pressure_of_dry_air',                                         'lev',  'avg', '1')
      call history_add_field('EXNER',     'reciprocal_of_dimensionless_exner_function_wrt_surface_air_pressure','lev',  'avg', '1')
      call history_add_field('ZM',        'geopotential_height_wrt_surface',                                    'lev',  'avg', 'm')
      call history_add_field('PINT',      'air_pressure_at_interfaces',                                         'ilev', 'avg', 'Pa')
      call history_add_field('PINTDRY',   'air_pressure_of_dry_air_at_interfaces',                              'ilev', 'avg', 'Pa')
      call history_add_field('LNPINT',    'ln_air_pressure_at_interfaces',                                      'ilev', 'avg', '1')
      call history_add_field('LNPINTDRY', 'ln_air_pressure_of_dry_air_at_interfaces',                           'ilev', 'avg', '1')
      call history_add_field('ZI',        'geopotential_height_wrt_surface_at_interfaces',                      'ilev', 'avg', 'm')

      ! Add surface fluxes received from the coupler
      call history_add_field('SHFLX',     'surface_upward_sensible_heat_flux_from_coupler',                     horiz_only, 'avg', 'W m-2')
      call history_add_field('LHFLX',     'surface_upward_latent_heat_flux_from_coupler',                       horiz_only, 'avg', 'W m-2')
      call history_add_field('TAUX',      'surface_eastward_wind_stress_from_coupler',                          horiz_only, 'avg', 'N m-2')
      call history_add_field('TAUY',      'surface_northward_wind_stress_from_coupler',                         horiz_only, 'avg', 'N m-2')

      ! Add expected constituent fields
      const_num_found = 0
      const_found = .false.
      wv_idx = 0
      do const_idx = 1, size(const_props)
         call const_props(const_idx)%standard_name(standard_name, errflg, errmsg)
         if (errflg /= 0) then
            return
         end if
         if (trim(standard_name) == wv_std_name) then
            ! Water vapor: also its surface flux and column integral
            wv_idx = const_idx
            call history_add_field('QFLX', 'surface_upward_water_vapor_flux_from_coupler', horiz_only, 'avg', 'kg m-2 s-1')
            call history_add_field('TMQ',  'vertically_integrated_water_vapor',            horiz_only, 'avg', 'kg m-2')
         end if
         do name_idx = 1, size(const_std_names)
            if (trim(standard_name) == trim(const_std_names(name_idx))) then
               call history_add_field(trim(const_diag_names(name_idx)), trim(const_std_names(name_idx)), 'lev', 'avg', 'kg kg-1', mixing_ratio='wet')
               const_num_found = const_num_found + 1
               const_found(const_idx) = .true.
            end if
         end do
         if (const_num_found == size(const_std_names)) then
            exit
         end if
      end do

      ! Add fields for all other constituents
      do const_idx = 1, size(const_props)
         if (.not. const_found(const_idx)) then
            call const_props(const_idx)%diagnostic_name(diagnostic_name, errflg, errmsg)
            if (errflg /= 0) then
               return
            end if
            call const_props(const_idx)%units(units, errflg, errmsg)
            if (errflg /= 0) then
               return
            end if
            call const_props(const_idx)%is_wet(is_wet, errflg, errmsg)
            if (errflg /= 0) then
               return
            end if
            call const_props(const_idx)%long_name(long_name, errflg, errmsg)
            if (errflg /= 0) then
               return
            end if
            if (is_wet) then
               mixing_ratio_type = 'wet'
            else
               mixing_ratio_type = 'dry'
            end if
            call history_add_field(trim(diagnostic_name), trim(long_name), 'lev', 'avg', trim(units), mixing_ratio=mixing_ratio_type)
         end if
      end do

   end subroutine sima_state_diagnostics_init

   !> \section arg_table_sima_state_diagnostics_run  Argument Table
   !! \htmlinclude sima_state_diagnostics_run.html
   subroutine sima_state_diagnostics_run(ps, psdry, phis, T, u, v, dse, omega, &
        pmid, pmiddry, pdel, pdeldry, rpdel, rpdeldry, lnpmid, lnpmiddry,     &
        inv_exner, zm, pint, pintdry, lnpint, lnpintdry, zi,                  &
        shf, lhf, cflx, wsx, wsy, gravit, const_array,                        &
        const_props, errmsg, errflg)

      use cam_history, only: history_out_field
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      ! State variables
      real(kind_phys), intent(in) :: ps(:)          ! surface pressure
      real(kind_phys), intent(in) :: psdry(:)       ! surface pressure of dry air
      real(kind_phys), intent(in) :: phis(:)        ! surface geopotential
      real(kind_phys), intent(in) :: T(:,:)         ! air temperature
      real(kind_phys), intent(in) :: u(:,:)         ! eastward wind (x wind)
      real(kind_phys), intent(in) :: v(:,:)         ! northward wind (y wind)
      real(kind_phys), intent(in) :: dse(:,:)       ! dry static energy
      real(kind_phys), intent(in) :: omega(:,:)     ! lagrangian tendency of air pressure
      real(kind_phys), intent(in) :: pmid(:,:)      ! air pressure
      real(kind_phys), intent(in) :: pmiddry(:,:)   ! air pressure of dry air
      real(kind_phys), intent(in) :: pdel(:,:)      ! air pressure thickness
      real(kind_phys), intent(in) :: pdeldry(:,:)   ! air pressure thickness of dry air
      real(kind_phys), intent(in) :: rpdel(:,:)     ! reciprocal of air pressure thickness
      real(kind_phys), intent(in) :: rpdeldry(:,:)  ! reciprocal of air pressure thickness of dry air
      real(kind_phys), intent(in) :: lnpmid(:,:)    ! ln air pressure
      real(kind_phys), intent(in) :: lnpmiddry(:,:) ! ln air pressure of dry air
      real(kind_phys), intent(in) :: inv_exner(:,:) ! inverse exner function wrt surface pressure
      real(kind_phys), intent(in) :: zm(:,:)        ! geopotential height wrt surface
      real(kind_phys), intent(in) :: pint(:,:)      ! air pressure at interfaces
      real(kind_phys), intent(in) :: pintdry(:,:)   ! air pressure of dry air at interfaces
      real(kind_phys), intent(in) :: lnpint(:,:)    ! ln air pressure at interfaces
      real(kind_phys), intent(in) :: lnpintdry(:,:) ! ln air pressure of dry air at interfaces
      real(kind_phys), intent(in) :: zi(:,:)        ! geopotential height wrt surface at interfaces
      ! Surface fluxes from the coupler
      real(kind_phys), intent(in) :: shf(:)         ! surface upward sensible heat flux
      real(kind_phys), intent(in) :: lhf(:)         ! surface upward latent heat flux
      real(kind_phys), intent(in) :: cflx(:,:)      ! surface upward constituent fluxes
      real(kind_phys), intent(in) :: wsx(:)         ! surface eastward wind stress
      real(kind_phys), intent(in) :: wsy(:)         ! surface northward wind stress
      real(kind_phys), intent(in) :: gravit         ! gravitational acceleration
      ! Constituent variables
      real(kind_phys), intent(in) :: const_array(:,:,:)
      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      character(len=512) :: standard_name
      character(len=fieldname_len)  :: diagnostic_name
      integer :: const_idx, name_idx
      integer :: const_num_found
      logical :: const_found(size(const_props))
      integer :: k
      real(kind_phys) :: tmq(size(pdel, 1))         ! vertically integrated water vapor

      errmsg = ''
      errflg = 0

      ! Capture state fields
      call history_out_field('PS'       , ps)
      call history_out_field('PSDRY'    , psdry)
      call history_out_field('PHIS'     , phis)
      call history_out_field('T'        , T)
      call history_out_field('U'        , u)
      call history_out_field('V'        , v)
      call history_out_field('DSE'      , dse)
      call history_out_field('OMEGA'    , omega)
      call history_out_field('PMID'     , pmid)
      call history_out_field('PMIDDRY'  , pmiddry)
      call history_out_field('PDEL'     , pdel)
      call history_out_field('PDELDRY'  , pdeldry)
      call history_out_field('RPDEL'    , rpdel)
      call history_out_field('RPDELDRY' , rpdeldry)
      call history_out_field('LNPMID'   , lnpmid)
      call history_out_field('LNPMIDDRY', lnpmiddry)
      call history_out_field('EXNER'    , inv_exner)
      call history_out_field('ZM'       , zm)
      call history_out_field('PINT'     , pint)
      call history_out_field('PINTDRY'  , pintdry)
      call history_out_field('LNPINT'   , lnpint)
      call history_out_field('LNPINTDRY', lnpintdry)
      call history_out_field('ZI'       , zi)

      ! Capture surface fluxes from the coupler
      call history_out_field('SHFLX'    , shf)
      call history_out_field('LHFLX'    , lhf)
      call history_out_field('TAUX'     , wsx)
      call history_out_field('TAUY'     , wsy)
      if (wv_idx > 0) then
         call history_out_field('QFLX', cflx(:, wv_idx))
         ! Column water vapor, as in CAM: sum of q * pdel / g over the layers
         tmq(:) = 0._kind_phys
         do k = 1, size(pdel, 2)
            tmq(:) = tmq(:) + const_array(:, k, wv_idx) * pdel(:, k)
         end do
         tmq(:) = tmq(:) / gravit
         call history_out_field('TMQ', tmq)
      end if

      ! Capture expected constituent fields
      const_num_found = 0
      const_found = .false.
      do const_idx = 1, size(const_props)
         call const_props(const_idx)%standard_name(standard_name, errflg, errmsg)
         if (errflg /= 0) then
            return
         end if
         do name_idx = 1, size(const_std_names)
            if (trim(standard_name) == trim(const_std_names(name_idx))) then
               call history_out_field(trim(const_diag_names(name_idx)), const_array(:,:,const_idx))
               const_num_found = const_num_found + 1
               const_found(const_idx) = .true.
            end if
         end do
         if (const_num_found == size(const_std_names)) then
            exit
         end if
      end do

      ! Capture all other constituent fields
      do const_idx = 1, size(const_props)
         if (.not. const_found(const_idx)) then
            ! Must match the name registered in sima_state_diagnostics_init.
            call const_props(const_idx)%diagnostic_name(diagnostic_name, errflg, errmsg)
            if (errflg /= 0) then
               return
            end if
            call history_out_field(trim(diagnostic_name), const_array(:,:,const_idx))
         end if
      end do

   end subroutine sima_state_diagnostics_run
end module sima_state_diagnostics
