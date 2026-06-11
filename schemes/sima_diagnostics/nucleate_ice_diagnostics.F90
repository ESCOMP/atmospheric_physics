! Diagnostics for nucleate_ice_ccpp scheme
module nucleate_ice_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: nucleate_ice_diagnostics_init
   public :: nucleate_ice_diagnostics_run

contains

   !> \section arg_table_nucleate_ice_diagnostics_init  Argument Table
   !! \htmlinclude nucleate_ice_diagnostics_init.html
   subroutine nucleate_ice_diagnostics_init(errmsg, errflg)
      use cam_history, only: history_add_field

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! Ice nucleation number tendencies by mechanism
      call history_add_field('NIHFTEN',  'tendency_of_ice_nucleation_number_concentration_from_homogeneous_freezing', &
           'lev', 'avg', 'm-3 s-1')
      call history_add_field('NIIMMTEN', 'tendency_of_ice_nucleation_number_concentration_from_immersion_freezing', &
           'lev', 'avg', 'm-3 s-1')
      call history_add_field('NIDEPTEN', 'tendency_of_ice_nucleation_number_concentration_from_deposition_nucleation', &
           'lev', 'avg', 'm-3 s-1')
      call history_add_field('NIMEYTEN', 'tendency_of_ice_nucleation_number_concentration_from_meyers_nucleation', &
           'lev', 'avg', 'm-3 s-1')

      ! Regime, subgrid, and tropopause diagnostics
      call history_add_field('NIREGM',    'ice_nucleation_regime_temperature_threshold',                    'lev', 'avg', 'degC')
      call history_add_field('NISUBGRID', 'subgrid_saturation_scaling_factor_diagnostic_for_ice_nucleation', 'lev', 'avg', '1')
      call history_add_field('NITROP_PD', 'vertical_layer_indicator_of_chemical_tropopause',                'lev', 'avg', '1')

      ! Pre-existing ice diagnostics
      call history_add_field('FHOM', 'fraction_of_cirrus_with_homogeneous_freezing',                      'lev', 'avg', 'fraction')
      call history_add_field('WICE', 'vertical_velocity_reduction_from_preexisting_ice_for_ice_nucleation', 'lev', 'avg', 'm s-1')
      call history_add_field('WEFF', 'effective_vertical_velocity_for_ice_nucleation',                     'lev', 'avg', 'm s-1')

      ! Aerosol input number concentration tendencies to ice nucleation
      call history_add_field('INnso4TEN',   'tendency_of_sulfate_number_concentration_input_to_ice_nucleation',     'lev', 'avg', 'm-3 s-1')
      call history_add_field('INnbcTEN',    'tendency_of_black_carbon_number_concentration_input_to_ice_nucleation', 'lev', 'avg', 'm-3 s-1')
      call history_add_field('INndustTEN',  'tendency_of_dust_number_concentration_input_to_ice_nucleation',        'lev', 'avg', 'm-3 s-1')
      call history_add_field('INondustTEN', 'tendency_of_dust_number_concentration_output_from_ice_nucleation',     'lev', 'avg', 'm-3 s-1')
      call history_add_field('INhetTEN',    'tendency_of_heterogeneous_ice_nuclei_number_concentration',            'lev', 'avg', 'm-3 s-1')
      call history_add_field('INhomTEN',    'tendency_of_homogeneous_ice_nuclei_number_concentration',              'lev', 'avg', 'm-3 s-1')

      ! Ice nucleation frequency diagnostics
      call history_add_field('INFrehom', 'frequency_of_homogeneous_ice_nucleation', 'lev', 'avg', '1')
      call history_add_field('INFreIN',  'frequency_of_ice_nucleation',             'lev', 'avg', '1')

   end subroutine nucleate_ice_diagnostics_init

   !> \section arg_table_nucleate_ice_diagnostics_run  Argument Table
   !! \htmlinclude nucleate_ice_diagnostics_run.html
   subroutine nucleate_ice_diagnostics_run( &
      nihf, niimm, nidep, nimey, &
      regm, subgrid_diag, trop_pd, &
      fhom, wice, weff, &
      INnso4, INnbc, INndust, INondust, &
      INhet, INhom, INFrehom, INFreIN, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      real(kind_phys),    intent(in) :: nihf(:,:)
      real(kind_phys),    intent(in) :: niimm(:,:)
      real(kind_phys),    intent(in) :: nidep(:,:)
      real(kind_phys),    intent(in) :: nimey(:,:)
      real(kind_phys),    intent(in) :: regm(:,:)
      real(kind_phys),    intent(in) :: subgrid_diag(:,:)
      real(kind_phys),    intent(in) :: trop_pd(:,:)
      real(kind_phys),    intent(in) :: fhom(:,:)
      real(kind_phys),    intent(in) :: wice(:,:)
      real(kind_phys),    intent(in) :: weff(:,:)
      real(kind_phys),    intent(in) :: INnso4(:,:)
      real(kind_phys),    intent(in) :: INnbc(:,:)
      real(kind_phys),    intent(in) :: INndust(:,:)
      real(kind_phys),    intent(in) :: INondust(:,:)
      real(kind_phys),    intent(in) :: INhet(:,:)
      real(kind_phys),    intent(in) :: INhom(:,:)
      real(kind_phys),    intent(in) :: INFrehom(:,:)
      real(kind_phys),    intent(in) :: INFreIN(:,:)
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call history_out_field('NIHFTEN',  nihf)
      call history_out_field('NIIMMTEN', niimm)
      call history_out_field('NIDEPTEN', nidep)
      call history_out_field('NIMEYTEN', nimey)

      call history_out_field('NIREGM',    regm)
      call history_out_field('NISUBGRID', subgrid_diag)
      call history_out_field('NITROP_PD', trop_pd)

      call history_out_field('FHOM', fhom)
      call history_out_field('WICE', wice)
      call history_out_field('WEFF', weff)

      call history_out_field('INnso4TEN',   INnso4)
      call history_out_field('INnbcTEN',    INnbc)
      call history_out_field('INndustTEN',  INndust)
      call history_out_field('INondustTEN', INondust)
      call history_out_field('INhetTEN',    INhet)
      call history_out_field('INhomTEN',    INhom)

      call history_out_field('INFrehom', INFrehom)
      call history_out_field('INFreIN',  INFreIN)

   end subroutine nucleate_ice_diagnostics_run

end module nucleate_ice_diagnostics
