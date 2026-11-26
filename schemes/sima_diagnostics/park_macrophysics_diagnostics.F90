module park_macrophysics_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private
   save

   public :: park_macrophysics_diagnostics_init
   public :: park_macrophysics_diagnostics_run

contains

   !> \section arg_table_park_macrophysics_diagnostics_init  Argument Table
   !! \htmlinclude park_macrophysics_diagnostics_init.html
   subroutine park_macrophysics_diagnostics_init(errmsg, errflg)
      use cam_history, only: history_add_field

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! Detrainment diagnostics
      call history_add_field('DPDLFLIQ', 'Detrained liquid water from deep convection', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('DPDLFICE', 'Detrained ice from deep convection', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('SHDLFLIQ', 'Detrained liquid water from shallow convection', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('SHDLFICE', 'Detrained ice from shallow convection', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('DPDLFT', 'T-tendency due to deep convective detrainment', &
                             'lev', 'avg', 'K s-1')
      call history_add_field('SHDLFT', 'T-tendency due to shallow convective detrainment', &
                             'lev', 'avg', 'K s-1')

      ! Macrophysics tendency diagnostics
      call history_add_field('MACPDT', 'Heating tendency - Park macrophysics', &
                             'lev', 'avg', 'W kg-1')
      call history_add_field('MACPDQ', 'Q tendency - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('MACPDLIQ', 'CLDLIQ tendency - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('MACPDICE', 'CLDICE tendency - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDVAPADJ', 'Q tendency - liq/ice adjustment - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDLIQADJ', 'CLDLIQ adjustment tendency - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDICEADJ', 'CLDICE adjustment tendency - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDLIQDET', 'Detrainment of conv cld liq into environment - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDICEDET', 'Detrainment of conv cld ice into environment - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDLIQLIM', 'CLDLIQ limiting tendency - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')
      call history_add_field('CLDICELIM', 'CLDICE limiting tendency - Park macrophysics', &
                             'lev', 'avg', 'kg kg-1 s-1')

      ! Cloud fraction and state diagnostics
      call history_add_field('AST', 'Stratus cloud fraction', &
                             'lev', 'avg', 'fraction')
      call history_add_field('LIQCLDF', 'Stratus liquid cloud fraction', &
                             'lev', 'avg', 'fraction')
      call history_add_field('ICECLDF', 'Stratus ICE cloud fraction', &
                             'lev', 'avg', 'fraction')
      call history_add_field('QLST', 'Physical in-stratus LWC', &
                             'lev', 'avg', 'kg kg-1')
      call history_add_field('QIST', 'Physical in-stratus IWC', &
                             'lev', 'avg', 'kg kg-1')
      call history_add_field('CMELIQ', 'Rate of cond-evap of liq within the cloud', &
                             'lev', 'avg', 'kg kg-1 s-1')

      ! Critical relative humidity diagnostics
      call history_add_field('RHMIN_LIQ', 'Default critical RH for liquid-stratus', &
                             'lev', 'avg', 'fraction')
      call history_add_field('RHMIN_ICE', 'Default critical RH for ice-stratus', &
                             'lev', 'avg', 'fraction')

   end subroutine park_macrophysics_diagnostics_init

   !> \section arg_table_park_macrophysics_diagnostics_run  Argument Table
   !! \htmlinclude park_macrophysics_diagnostics_run.html
   subroutine park_macrophysics_diagnostics_run( &
      dpdlfliq, dpdlfice, shdlfliq, shdlfice, dpdlft, shdlft, dlf, &
      tlat, qvlat, qcten, qiten, &
      qvadj, qladj, qiadj, dlf_ql, dlf_qi, qllim, qilim, &
      ast, alst, aist, qlst, qist, cmeliq, &
      rhmin_liq, rhmin_ice, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      ! Detrainment fields
      real(kind_phys), intent(in) :: dpdlfliq(:,:)
      real(kind_phys), intent(in) :: dpdlfice(:,:)
      real(kind_phys), intent(in) :: shdlfliq(:,:)
      real(kind_phys), intent(in) :: shdlfice(:,:)
      real(kind_phys), intent(in) :: dpdlft(:,:)
      real(kind_phys), intent(in) :: shdlft(:,:)
      real(kind_phys), intent(in) :: dlf(:,:)

      ! Macrophysics tendency fields
      real(kind_phys), intent(in) :: tlat(:,:)
      real(kind_phys), intent(in) :: qvlat(:,:)
      real(kind_phys), intent(in) :: qcten(:,:)
      real(kind_phys), intent(in) :: qiten(:,:)
      real(kind_phys), intent(in) :: qvadj(:,:)
      real(kind_phys), intent(in) :: qladj(:,:)
      real(kind_phys), intent(in) :: qiadj(:,:)
      real(kind_phys), intent(in) :: dlf_ql(:,:)
      real(kind_phys), intent(in) :: dlf_qi(:,:)
      real(kind_phys), intent(in) :: qllim(:,:)
      real(kind_phys), intent(in) :: qilim(:,:)

      ! Cloud fraction and state fields
      real(kind_phys), intent(in) :: ast(:,:)
      real(kind_phys), intent(in) :: alst(:,:)
      real(kind_phys), intent(in) :: aist(:,:)
      real(kind_phys), intent(in) :: qlst(:,:)
      real(kind_phys), intent(in) :: qist(:,:)
      real(kind_phys), intent(in) :: cmeliq(:,:)

      ! Critical relative humidity fields
      real(kind_phys), intent(in) :: rhmin_liq(:,:)
      real(kind_phys), intent(in) :: rhmin_ice(:,:)

      ! CCPP error handling
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! Detrainment diagnostics
      call history_out_field('DPDLFLIQ', dpdlfliq)
      call history_out_field('DPDLFICE', dpdlfice)
      call history_out_field('SHDLFLIQ', shdlfliq)
      call history_out_field('SHDLFICE', shdlfice)
      call history_out_field('DPDLFT', dpdlft)
      call history_out_field('SHDLFT', shdlft)

      ! Macrophysics tendency diagnostics
      call history_out_field('MACPDT', tlat)
      call history_out_field('MACPDQ', qvlat)
      call history_out_field('MACPDLIQ', qcten)
      call history_out_field('MACPDICE', qiten)
      call history_out_field('CLDVAPADJ', qvadj)
      call history_out_field('CLDLIQADJ', qladj)
      call history_out_field('CLDICEADJ', qiadj)
      call history_out_field('CLDLIQDET', dlf_ql)
      call history_out_field('CLDICEDET', dlf_qi)
      call history_out_field('CLDLIQLIM', qllim)
      call history_out_field('CLDICELIM', qilim)

      ! Cloud fraction and state diagnostics
      call history_out_field('AST', ast)
      call history_out_field('LIQCLDF', alst)
      call history_out_field('ICECLDF', aist)
      call history_out_field('QLST', qlst)
      call history_out_field('QIST', qist)
      call history_out_field('CMELIQ', cmeliq)

      ! Critical relative humidity diagnostics
      call history_out_field('RHMIN_LIQ', rhmin_liq)
      call history_out_field('RHMIN_ICE', rhmin_ice)

   end subroutine park_macrophysics_diagnostics_run

end module park_macrophysics_diagnostics
