! Diagnostics for RK stratiform - cloud particle sedimentation
module cloud_particle_sedimentation_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: cloud_particle_sedimentation_diagnostics_init
   public :: cloud_particle_sedimentation_diagnostics_run

contains

   !> \section arg_table_cloud_particle_sedimentation_diagnostics_init  Argument Table
   !! \htmlinclude cloud_particle_sedimentation_diagnostics_init.html
   subroutine cloud_particle_sedimentation_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('DQSED', 'tendency_of_water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_sedimentation', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('DISED', 'tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_sedimentation', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('DLSED', 'tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_sedimentation', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('HSED', 'tendency_of_dry_air_enthalpy_at_constant_pressure_due_to_sedimentation', 'lev', 'avg', 'J kg-1 s-1')
      call history_add_field('PRECSED', 'stratiform_cloud_water_surface_flux_due_to_sedimentation', horiz_only, 'inst', 'm s-1')
      call history_add_field('SNOWSED', 'stratiform_lwe_cloud_ice_surface_flux_due_to_sedimentation', horiz_only, 'inst', 'm s-1')
      call history_add_field('RAINSED', 'stratiform_rain_flux_at_surface_due_to_sedimentation', horiz_only, 'inst', 'm s-1')

   end subroutine cloud_particle_sedimentation_diagnostics_init

   !> \section arg_table_cloud_particle_sedimentation_diagnostics_run  Argument Table
   !! \htmlinclude cloud_particle_sedimentation_diagnostics_run.html
   subroutine cloud_particle_sedimentation_diagnostics_run( &
      ncol, &
      wvtend, icetend, liqtend, htend, &
      snow_sed, sfliq, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      integer,            intent(in)  :: ncol
      real(kind_phys),    intent(in)  :: wvtend(:,:)    ! water vapor tendency -- to apply wv tendency
      real(kind_phys),    intent(in)  :: icetend(:,:)   ! ice condensate tendency -- to apply cldice tendency
      real(kind_phys),    intent(in)  :: liqtend(:,:)   ! liquid condensate tendency -- to apply cldliq tendency
      real(kind_phys),    intent(in)  :: htend(:,:)     ! heating rate [J kg-1 s-1] -- to apply s tendency

      real(kind_phys),    intent(in)  :: snow_sed(:)    ! stratiform_lwe_cloud_ice_surface_flux_due_to_sedimentation [m s-1]
      real(kind_phys),    intent(in)  :: sfliq(:)       ! stratiform_rain_flux_at_surface_due_to_sedimentation [kg m-2 s-1]

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      real(kind_phys),    parameter   :: rhofw = 1000._kind_phys  ! density of fresh water [kg m-3]
      real(kind_phys) :: prec_sed(ncol)

      ! repeat computation of prec_sed here for diagnostics [m s-1]
      prec_sed(:ncol) = sfliq(:ncol)/rhofw + snow_sed(:ncol)

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('DQSED'  , wvtend)
      call history_out_field('DISED'  , icetend)
      call history_out_field('DLSED'  , liqtend)
      call history_out_field('HSED'   , htend)

      call history_out_field('PRECSED', prec_sed) ! calculated as m s-1
      call history_out_field('SNOWSED', snow_sed) ! already in m s-1
      call history_out_field('RAINSED', sfliq/rhofw) ! convert from kg m-2 s-1 to m s-1 (precip units) for output

   end subroutine cloud_particle_sedimentation_diagnostics_run

end module cloud_particle_sedimentation_diagnostics
