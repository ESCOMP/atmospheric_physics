module cam_diagnostics

   use ccpp_kinds, only:  kind_phys
   use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

   implicit none
   private
   save

   public :: cam_state_diagnostics_init ! init routine
   public :: cam_state_diagnostics_run  ! main routine
   public :: cam_tend_diagnostics_init ! init routine
   public :: cam_tend_diagnostics_run  ! main routine

   character(len=65) :: const_std_names(4) = &
   (/'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water       ', &
     'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', &
     'rain_mixing_ratio_wrt_moist_air_and_condensed_water              ', &
     'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water         '/)

   character(len=6) :: const_diag_names(4) = (/'Q     ', &
                                               'CLDLIQ', &
                                               'CLDICE', &
                                               'RAINQM'/)

CONTAINS

   !> \section arg_table_cam_state_diagnostics_init  Argument Table
   !! \htmlinclude cam_state_diagnostics_init.html
   subroutine cam_state_diagnostics_init(const_props, errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      integer :: const_idx, name_idx
      character(len=512) :: standard_name

      errmsg = ''
      errflg = 0

      ! Add state fields
      call history_add_field('PS',        'surface_pressure',                           horiz_only,  'avg', 'Pa')
      call history_add_field('PSDRY',     'surface_pressure_of_dry_air',                horiz_only,  'avg', 'Pa')
      call history_add_field('PHIS',      'surface_geopotential',                       horiz_only,  'lst', 'Pa')
      call history_add_field('T',         'air_temperature',                                 'lev',  'avg', 'K')
      call history_add_field('U',         'eastward_wind',                                   'lev',  'avg', 'm s-1')
      call history_add_field('V',         'northward_wind',                                  'lev',  'avg', 'm s-1')
      call history_add_field('DSE',       'dry_static_energy',                               'lev',  'avg', 'm s-1')
      call history_add_field('OMEGA',     'lagrangian_tendency_of_air_pressure',             'lev',  'avg', 'Pa s-1')
      call history_add_field('PMID',      'air_pressure',                                    'lev',  'avg', 'Pa')
      call history_add_field('PMIDDRY',   'air_pressure_of_dry_air',                         'lev',  'avg', 'Pa')
      call history_add_field('PDEL',      'air_pressure_thickness',                          'lev',  'avg', 'Pa')
      call history_add_field('PDELDRY',   'air_pressure_thickness_of_dry_air',               'lev',  'avg', 'Pa')
      call history_add_field('RPDEL',     'reciprocal_of_air_pressure_thickness',            'lev',  'avg', 'Pa-1')
      call history_add_field('RPDELDRY',  'reciprocal_of_air_pressure_thickness_of_dry_air', 'lev',  'avg', 'Pa-1')
      call history_add_field('LNPMID',    'ln_air_pressure',                                 'lev',  'avg', '1')
      call history_add_field('LNPMIDDRY', 'ln_air_pressure_of_dry_air',                      'lev',  'avg', '1')
      call history_add_field('EXNER',     'inverse_exner_function_wrt_surface_pressure',     'lev',  'avg', '1')
      call history_add_field('ZM',        'geopotential_height_wrt_surface',                 'lev',  'avg', 'm')
      call history_add_field('PINT',      'air_pressure_at_interface',                       'ilev', 'avg', 'Pa')
      call history_add_field('PINTDRY',   'air_pressure_of_dry_air_at_interface',            'ilev', 'avg', 'Pa')
      call history_add_field('LNPINT',    'ln_air_pressure_at_interface',                    'ilev', 'avg', '1')
      call history_add_field('LNPINTDRY', 'ln_air_pressure_of_dry_air_at_interface',         'ilev', 'avg', '1')
      call history_add_field('ZI',        'geopotential_heigh_wrt_surface_at_interface',     'ilev', 'avg', 'm')
      ! Add constituent fields
      do const_idx = 1, size(const_props)
         call const_props(const_idx)%standard_name(standard_name, errflg, errmsg)
         do name_idx = 1, size(const_std_names)
            if (trim(standard_name) == trim(const_std_names(name_idx))) then
               call history_add_field(trim(const_diag_names(name_idx)), trim(const_std_names(name_idx)), 'lev', 'avg', 'kg kg-1', mixing_ratio='wet')
            end if
         end do
      end do

   end subroutine cam_state_diagnostics_init

   !> \section arg_table_cam_state_diagnostics_run  Argument Table
   !! \htmlinclude cam_state_diagnostics_run.html
   subroutine cam_state_diagnostics_run(ps, psdry, phis, T, u, v, dse, omega, &
        pmid, pmiddry, pdel, pdeldry, rpdel, rpdeldry, lnpmid, lnpmiddry,     &
        inv_exner, zm, pint, pintdry, lnpint, lnpintdry, zi, const_array,     &
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
      real(kind_phys), intent(in) :: pint(:,:)      ! air pressure at interface
      real(kind_phys), intent(in) :: pintdry(:,:)   ! air pressure of dry air at interface
      real(kind_phys), intent(in) :: lnpint(:,:)    ! ln air pressure at interface
      real(kind_phys), intent(in) :: lnpintdry(:,:) ! ln air pressure of dry air at interface
      real(kind_phys), intent(in) :: zi(:,:)        ! geopotential height wrt surface at interface
      ! Constituent variables
      real(kind_phys), intent(in) :: const_array(:,:,:)
      type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      character(len=512) :: standard_name
      integer :: const_idx, name_idx

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

      ! Capture constituent fields
      do const_idx = 1, size(const_props)
         call const_props(const_idx)%standard_name(standard_name, errflg, errmsg)
         do name_idx = 1, size(const_std_names)
            if (trim(standard_name) == trim(const_std_names(name_idx))) then
               call history_out_field(trim(const_diag_names(name_idx)), const_array(:,:,const_idx))
            end if
         end do
      end do

   end subroutine cam_state_diagnostics_run

   !> \section arg_table_cam_tend_diagnostics_init  Argument Table
   !! \htmlinclude cam_tend_diagnostics_init.html
   subroutine cam_tend_diagnostics_init(errmsg, errflg)
      use cam_history, only: history_add_field
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Add tendency fields
      call history_add_field('TTEND', 'tendency_of_air_temperature_due_to_model_physics', 'lev', 'avg', 'K s-1')
      call history_add_field('UTEND', 'tendency_of_eastward_wind_due_to_model_physics',   'lev', 'avg', 'm s-2')
      call history_add_field('VTEND', 'tendency_of_northward_wind_due_to_model_physics',  'lev', 'avg', 'm s-2')

   end subroutine cam_tend_diagnostics_init

   !> \section arg_table_cam_tend_diagnostics_run  Argument Table
   !! \htmlinclude cam_tend_diagnostics_run.html
   subroutine cam_tend_diagnostics_run(dTdt_total, dudt_total, dvdt_total, errmsg, errflg)
      use cam_history, only: history_out_field
      ! Tendency variables
      real(kind_phys), intent(in) :: dTdt_total(:,:) ! tendency of air temperature due to model physics
      real(kind_phys), intent(in) :: dudt_total(:,:) ! tendency of eastward wind due to model physics
      real(kind_phys), intent(in) :: dvdt_total(:,:) ! tendency of northward wind due to model physics
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Capture tendency fields
      call history_out_field('TTEND', dTdt_total)
      call history_out_field('UTEND', dudt_total)
      call history_out_field('VTEND', dvdt_total)

   end subroutine cam_tend_diagnostics_run
   !=======================================================================
end module cam_diagnostics
