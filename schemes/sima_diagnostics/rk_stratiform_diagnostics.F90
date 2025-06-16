! Diagnostics for RK stratiform - miscellaneous interstitial schemes
module rk_stratiform_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: rk_stratiform_diagnostics_init
   public :: rk_stratiform_cloud_fraction_perturbation_diagnostics_run
   public :: rk_stratiform_condensate_repartioning_diagnostics_run
   public :: rk_stratiform_prognostic_cloud_water_tendencies_diagnostics_run
   public :: rk_stratiform_cloud_optical_properties_diagnostics_run

contains

   !> \section arg_table_rk_stratiform_diagnostics_init  Argument Table
   !! \htmlinclude rk_stratiform_diagnostics_init.html
   subroutine rk_stratiform_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! There is one initialization scheme for all RK diagnostics
      ! but there are separate run phases for diagnostics
      ! pertaining to each interstitial scheme. See RK SDF file.

      errmsg = ''
      errflg = 0

      ! rk_stratiform_cloud_fraction_perturbation_diagnostics
      call history_add_field('AST', 'cloud_area_fraction', 'lev', 'avg', 'fraction')

      ! rk_stratiform_condensate_repartioning_diagnostics
      call history_add_field('FICE', 'mass_fraction_of_ice_content_within_stratiform_cloud', 'lev', 'avg', 'fraction')
      call history_add_field('REPARTICE', 'tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_cloud_ice_and_cloud_liquid_repartitioning', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('REPARTLIQ', 'tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water_due_to_cloud_ice_and_cloud_liquid_repartitioning', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('HREPART', 'tendency_of_dry_air_enthalpy_at_constant_pressure_due_to_cloud_ice_and_cloud_liquid_repartitioning', 'lev', 'avg', 'J kg-1 s-1')

      ! rk_stratiform_prognostic_cloud_water_tendencies_diagnostics
      call history_add_field('FWAUT', 'relative_importance_of_cloud_liquid_water_autoconversion', 'lev', 'avg', 'fraction')
      call history_add_field('FSAUT', 'relative_importance_of_cloud_ice_autoconversion', 'lev', 'avg', 'fraction')
      call history_add_field('FRACW', 'relative_importance_of_rain_accreting_cloud_liquid_water', 'lev', 'avg', 'fraction')
      call history_add_field('FSACW', 'relative_importance_of_snow_accreting_cloud_liquid_water', 'lev', 'avg', 'fraction')
      call history_add_field('FSACI', 'relative_importance_of_snow_accreting_cloud_ice', 'lev', 'avg', 'fraction')
      call history_add_field('PCSNOW', 'lwe_snow_precipitation_rate_at_surface_due_to_microphysics', horiz_only, 'avg', 'fraction')
      call history_add_field('CME', 'net_condensation_rate_due_to_microphysics', 'lev', 'avg', 'kg kg-1 s-1') ! qme.
      call history_add_field('CMEICE', 'rate_of_condensation_minus_evaporation_for_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CMELIQ', 'rate_of_condensation_minus_evaporation_for_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('ICE2PR', 'tendency_of_cloud_ice_mixing_ratio_wrt_to_moist_air_and_condensed_water_due_to_ice_to_snow_autoconversion', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('LIQ2PR', 'tendency_of_cloud_liquid_water_mixing_ratio_wrt_to_moist_air_and_condensed_water_due_to_liquid_to_rain_autoconversion', 'lev', 'avg', 'kg kg-1 s-1')

      call history_add_field('HPROGCLD', 'tendency_of_dry_air_enthalpy_at_constant_pressure_due_to_microphysics', 'lev', 'avg', 'J kg-1 s-1')
      call history_add_field('HEVAP', 'tendency_of_dry_air_enthalpy_at_constant_pressure_due_to_evaporation_of_precipitation', 'lev', 'avg', 'J kg-1 s-1')
      call history_add_field('HMELT', 'tendency_of_dry_air_enthalpy_at_constant_pressure_due_to_snow_melt', 'lev', 'avg', 'J kg-1 s-1')
      call history_add_field('HCME', 'tendency_of_dry_air_enthalpy_at_constant_pressure_due_to_condensation_minus_evaporation', 'lev', 'avg', 'J kg-1 s-1')
      call history_add_field('HFREEZ', 'tendency_of_dry_air_enthalpy_at_constant_pressure_due_to_freezing_of_precipitation', 'lev', 'avg', 'J kg-1 s-1')

      call history_add_field('PRODPREC', 'precipitation_production_due_to_microphysics', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('EVAPPREC', 'rate_of_evaporation_of_precipitation_due_to_microphysics', 'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('EVAPSNOW', 'rate_of_evaporation_of_falling_snow_due_to_microphysics', 'lev', 'avg', 'kg kg-1 s-1')

      ! ... for COSP/CFMIP
      call history_add_field('LS_FLXPRC', 'stratiform_rain_and_snow_flux_at_interface', 'ilev', 'avg', 'kg m-2 s-1')
      call history_add_field('LS_FLXSNW', 'stratiform_snow_flux_at_interface', 'ilev', 'avg', 'kg m-2 s-1')
      call history_add_field('PRACWO', 'accretion_of_cloud_liquid_water_by_rain', 'lev', 'avg', 's-1')
      call history_add_field('PSACWO', 'accretion_of_cloud_liquid_water_by_snow', 'lev', 'avg', 's-1')
      call history_add_field('PSACIO', 'accretion_of_cloud_ice_by_snow', 'lev', 'avg', 's-1')

      call history_add_field('CLDLIQSTR', 'stratiform_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')
      call history_add_field('CLDICESTR', 'stratiform_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')
      call history_add_field('CLDLIQCON', 'convective_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')
      call history_add_field('CLDICECON', 'convective_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')

      ! rk_stratiform_cloud_optical_properties_diagnostics
      call history_add_field('IWC', 'stratiform_cloud_ice_water_content', 'lev', 'avg', 'kg m-3')
      call history_add_field('LWC', 'stratiform_cloud_liquid_water_content', 'lev', 'avg', 'kg m-3')
      call history_add_field('ICIMR', 'in_cloud_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')
      call history_add_field('ICWMR', 'in_cloud_cloud_liquid_mixing_ratio_wrt_moist_air_and_condensed_water', 'lev', 'avg', 'kg kg-1')

      call history_add_field('REI', 'effective_radius_of_stratiform_cloud_ice_particle', 'lev', 'avg', 'um')
      call history_add_field('REL', 'effective_radius_of_stratiform_cloud_liquid_water_particle', 'lev', 'avg', 'um')


   end subroutine rk_stratiform_diagnostics_init

   !> \section arg_table_rk_stratiform_cloud_fraction_perturbation_diagnostics_run  Argument Table
   !! \htmlinclude rk_stratiform_cloud_fraction_perturbation_diagnostics_run.html
   subroutine rk_stratiform_cloud_fraction_perturbation_diagnostics_run( &
      cloud, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      real(kind_phys),    intent(in) :: cloud(:, :)      ! cloud_area_fraction [fraction]

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('AST', cloud)

   end subroutine rk_stratiform_cloud_fraction_perturbation_diagnostics_run

   !> \section arg_table_rk_stratiform_condensate_repartioning_diagnostics_run  Argument Table
   !! \htmlinclude rk_stratiform_condensate_repartioning_diagnostics_run.html
   subroutine rk_stratiform_condensate_repartioning_diagnostics_run( &
      fice, tend_cldice, tend_cldliq, repartht, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      real(kind_phys),    intent(in)    :: fice(:,:)      ! mass_fraction_of_ice_content_within_stratiform_cloud [fraction]
      real(kind_phys),    intent(in)    :: tend_cldice(:,:)  ! tendency_of_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1 s-1]
      real(kind_phys),    intent(in)    :: tend_cldliq(:,:)  ! tendency_of_cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1 s-1]
      real(kind_phys),    intent(in)    :: repartht(:,:)     ! [J kg-1 s-1]

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('FICE', fice)
      call history_out_field('REPARTICE', tend_cldice)
      call history_out_field('REPARTLIQ', tend_cldliq)
      call history_out_field('HREPART', repartht)

   end subroutine rk_stratiform_condensate_repartioning_diagnostics_run

   !> \section arg_table_rk_stratiform_prognostic_cloud_water_tendencies_diagnostics_run  Argument Table
   !! \htmlinclude rk_stratiform_prognostic_cloud_water_tendencies_diagnostics_run.html
   subroutine rk_stratiform_prognostic_cloud_water_tendencies_diagnostics_run( &
      ncol, pver, &
      cloud, concld, &
      cldliq, cldice, &
      fwaut, fsaut, fracw, fsacw, fsaci, &
      snow_pcw, cme, cmeice, cmeliq, ice2pr, liq2pr, &
      tend_s, &
      evapheat, meltheat, cmeheat, prfzheat, &
      prodprec, evapprec, evapsnow, &
      lsflxprc, lsflxsnw, &
      pracwo, psacwo, psacio, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input parameters
      integer,            intent(in)    :: ncol
      integer,            intent(in)    :: pver
      real(kind_phys),    intent(in)    :: cloud(:,:)  ! cloud_area_fraction [fraction]
      real(kind_phys),    intent(in)    :: concld(:,:) ! convective_cloud_area_fraction [fraction]
      real(kind_phys),    intent(in)    :: cldliq(:,:) ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
      real(kind_phys),    intent(in)    :: cldice(:,:) ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]

      real(kind_phys),    intent(in)    :: fwaut(:,:)
      real(kind_phys),    intent(in)    :: fsaut(:,:)
      real(kind_phys),    intent(in)    :: fracw(:,:)
      real(kind_phys),    intent(in)    :: fsacw(:,:)
      real(kind_phys),    intent(in)    :: fsaci(:,:)

      real(kind_phys),    intent(in)    :: snow_pcw(:)
      real(kind_phys),    intent(in)    :: cme(:,:)
      real(kind_phys),    intent(in)    :: cmeice(:,:)
      real(kind_phys),    intent(in)    :: cmeliq(:,:)
      real(kind_phys),    intent(in)    :: ice2pr(:,:)
      real(kind_phys),    intent(in)    :: liq2pr(:,:)

      real(kind_phys),    intent(in)    :: tend_s(:,:)
      real(kind_phys),    intent(in)    :: evapheat(:,:)
      real(kind_phys),    intent(in)    :: meltheat(:,:)
      real(kind_phys),    intent(in)    :: cmeheat(:,:)
      real(kind_phys),    intent(in)    :: prfzheat(:,:)

      real(kind_phys),    intent(in)    :: prodprec(:,:)
      real(kind_phys),    intent(in)    :: evapprec(:,:)
      real(kind_phys),    intent(in)    :: evapsnow(:,:)
      real(kind_phys),    intent(in)    :: lsflxprc(:,:)
      real(kind_phys),    intent(in)    :: lsflxsnw(:,:)
      real(kind_phys),    intent(in)    :: pracwo(:,:)
      real(kind_phys),    intent(in)    :: psacwo(:,:)
      real(kind_phys),    intent(in)    :: psacio(:,:)

      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      real(kind_phys) :: cldliqstr(ncol, pver) ! [kg kg-1]
      real(kind_phys) :: cldicestr(ncol, pver) ! [kg kg-1]
      real(kind_phys) :: cldliqcon(ncol, pver) ! [kg kg-1]
      real(kind_phys) :: cldicecon(ncol, pver) ! [kg kg-1]
      integer :: i, k

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('FWAUT', fwaut)
      call history_out_field('FSAUT', fsaut)
      call history_out_field('FRACW', fracw)
      call history_out_field('FSACW', fsacw)
      call history_out_field('FSACI', fsaci)

      call history_out_field('PCSNOW', snow_pcw)
      call history_out_field('CME', cme)
      call history_out_field('CMEICE', cmeice)
      call history_out_field('CMELIQ', cmeliq)
      call history_out_field('ICE2PR', ice2pr)
      call history_out_field('LIQ2PR', liq2pr)

      call history_out_field('HPROGCLD', tend_s)
      call history_out_field('HEVAP', evapheat)
      call history_out_field('HMELT', meltheat)
      call history_out_field('HCME', cmeheat)
      call history_out_field('HFREEZ', prfzheat)

      call history_out_field('PRODPREC', prodprec)
      call history_out_field('EVAPPREC', evapprec)
      call history_out_field('EVAPSNOW', evapsnow)

      call history_out_field('LS_FLXPRC', lsflxprc)
      call history_out_field('LS_FLXSNW', lsflxsnw)
      call history_out_field('PRACWO', pracwo)
      call history_out_field('PSACWO', psacwo)
      call history_out_field('PSACIO', psacio)

      ! Derived diagnostics -- mass mixing ratio for stratiform or convective cloud liquid / cloud ice
      cldliqstr(:,:) = 0._kind_phys
      cldicestr(:,:) = 0._kind_phys
      cldliqcon(:,:) = 0._kind_phys
      cldicecon(:,:) = 0._kind_phys
      do k = 1, pver
         do i = 1, ncol
            if(cloud(i,k) > 0._kind_phys) then
               ! convective mass mixing ratios
               cldliqcon(i,k) = cldliq(i,k)/cloud(i,k) * concld(i,k)
               cldicecon(i,k) = cldice(i,k)/cloud(i,k) * concld(i,k)

               ! stratiform (large-scale) mass mixing ratios
               cldliqstr(i,k) = cldliq(i,k)/cloud(i,k) * (cloud(i,k) - concld(i,k))
               cldicestr(i,k) = cldice(i,k)/cloud(i,k) * (cloud(i,k) - concld(i,k))
            endif
         enddo
      enddo

      call history_out_field('CLDLIQCON', cldliqcon)
      call history_out_field('CLDICECON', cldicecon)

      call history_out_field('CLDLIQSTR', cldliqstr)
      call history_out_field('CLDICESTR', cldicestr)

   end subroutine rk_stratiform_prognostic_cloud_water_tendencies_diagnostics_run

   !> \section arg_table_rk_stratiform_cloud_optical_properties_diagnostics_run  Argument Table
   !! \htmlinclude rk_stratiform_cloud_optical_properties_diagnostics_run.html
   subroutine rk_stratiform_cloud_optical_properties_diagnostics_run( &
      ncol, pver, &
      rair, &
      pmid, &
      t, &
      cldice, cldliq, &
      rhcloud, &
      rel, rei, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Input arguments
      integer,            intent(in)    :: ncol
      integer,            intent(in)    :: pver
      real(kind_phys),    intent(in)    :: rair
      real(kind_phys),    intent(in)    :: pmid(:,:)      ! air_pressure [Pa]
      real(kind_phys),    intent(in)    :: t(:,:)         ! air_temperature [K]
      real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
      real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
      real(kind_phys),    intent(in)    :: rhcloud(:,:)   ! cloud_area_fraction_from_relative_humidity_method [fraction]
      real(kind_phys),    intent(in)    :: rel(:,:)       ! effective_radius_of_stratiform_cloud_liquid_water_particle [um]
      real(kind_phys),    intent(in)    :: rei(:,:)       ! effective_radius_of_stratiform_cloud_ice_particle [um]

      ! Output arguments
      character(len=512), intent(out)   :: errmsg         ! error message
      integer,            intent(out)   :: errflg         ! error flag

      ! Temporaries for diagnostic output.
      real(kind_phys) :: iwc(ncol,pver)       ! stratiform_cloud_ice_water_content [kg m-3]
      real(kind_phys) :: lwc(ncol,pver)       ! stratiform_cloud_liquid_water_content [kg m-3]
      real(kind_phys) :: icimr(ncol,pver)     ! in_cloud_cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
      real(kind_phys) :: icwmr(ncol,pver)     ! in_cloud_cloud_liquid_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]

      integer :: i, k


      ! Prognostic cloud water diagnostics
      ! needs updated cloud fraction
      do k = 1, pver
         do i = 1, ncol
            iwc(i,k)   = cldice(i,k)*pmid(i,k)/(rair*t(i,k))
            lwc(i,k)   = cldliq(i,k)*pmid(i,k)/(rair*t(i,k))
            icimr(i,k) = cldice(i,k) / max(0.01_kind_phys, rhcloud(i,k))
            icwmr(i,k) = cldliq(i,k) / max(0.01_kind_phys, rhcloud(i,k))
         end do
      end do

      call history_out_field('IWC', iwc)
      call history_out_field('LWC', lwc)
      call history_out_field('ICIMR', icimr)
      call history_out_field('ICWMR', icwmr)

      ! Cloud optical properties
      call history_out_field('REL', rel)
      call history_out_field('REI', rei)

   end subroutine rk_stratiform_cloud_optical_properties_diagnostics_run


end module rk_stratiform_diagnostics
