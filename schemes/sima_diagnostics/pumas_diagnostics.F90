!REVIEWERS - This is a work in progress.  Please save reviews until a future PR

!CACNOTE - All outfld calls to history_out_field are converted
!         - All history_add_field calls are done
!               - Except the constituent history_add_field needs to be changed - see CACNOTE
!               - Need to replace "#" in units?
!         - need to line up calls
!!!!----------------------------------------------

module pumas_diagnostics

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: pumas_diagnostics_init ! init routine

CONTAINS

!> \section arg_table_pumas_diagnostics_init  Argument Table
!! \htmlinclude pumas_diagnostics_init.html
subroutine pumas_diagnostics_init(qprops, ncnst, errmsg, errflg)

   use cam_history,         only: history_add_field
   use cam_history_support, only: horiz_only

   type(ccpp_constituent_prop_ptr_t), intent(in)  :: qprops(:)
   integer,             intent(in) :: ncnst          ! number of tracers to transport

   character(len=512), intent(out) :: errmsg
   integer,            intent(out) :: errflg

   ! Local variables:

   errmsg = ''
   errflg = 0

   ! NOTE -subcyc comments below mean that the original CAM code had "sampled_on_subcycle=.true." in the old outfld call
   !        Will use this for guidance during future developmen

!CACNOTE - come back to this - did not convert get the constituent info yet
!--------------------
!   do m = 1, ncnst
!      call cnst_get_ind(cnst_names(m), mm)
!      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixsnow, ixgraupel /)) ) then
!  mass mixing ratios
!         call history_add_field(cnst_name(mm), cnst_longname(mm),                    'lev',      'avg', 'kg kg-1') !subcyc
!         call history_add_field(sflxnam(mm),   trim(cnst_name(mm))//' surface flux', horiz_only, 'avg', 'kg m-2 s-1') !subcyc
!      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain, ixnumsnow, ixnumgraupel /)) ) then
!         ! number concentrations
!         call history_add_field(cnst_name(mm), cnst_longname(mm),                    'lev',      'avg', 'kg-1') !subcyc
!         call history_add_field(sflxnam(mm),  trim(cnst_name(mm))//' surface flux',  horiz_only, 'avg', '1 m-2 s-1') !subcyc
!      else
!         call endrun( "micro_pumas_cam_init: &
!              &Could not call history_add_field for constituent with unknown units.")
!      endif
!   end do

!   call history_add_field(apcnst(ixcldliq), trim(cnst_name(ixcldliq))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
!   call history_add_field(apcnst(ixcldice), trim(cnst_name(ixcldice))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
!   call history_add_field(bpcnst(ixcldliq), trim(cnst_name(ixcldliq))//' before physics','lev', 'avg', 'kg kg-1') !subcyc
!   call history_add_field(bpcnst(ixcldice), trim(cnst_name(ixcldice))//' before physics','lev', 'avg', 'kg kg-1') !subcyc

!   call history_add_field(apcnst(ixrain), trim(cnst_name(ixrain))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
!   call history_add_field(apcnst(ixsnow), trim(cnst_name(ixsnow))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
!   call history_add_field(bpcnst(ixrain), trim(cnst_name(ixrain))//' before physics','lev', 'avg', 'kg kg-1') !subcyc
!   call history_add_field(bpcnst(ixsnow), trim(cnst_name(ixsnow))//' before physics','lev', 'avg', 'kg kg-1') !subcyc

!   if (micro_mg_version > 2) then
!      call history_add_field(apcnst(ixgraupel), trim(cnst_name(ixgraupel))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
!      call history_add_field(bpcnst(ixgraupel), trim(cnst_name(ixgraupel))//' before physics','lev', 'avg', 'kg kg-1') !subcyc
!   end if

!--------------------

   call history_add_field ('CME',      'Rate of cond-evap within the cloud',        'lev',         'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('PRODPREC', 'Rate of conversion of condensate to precip','lev',         'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('EVAPPREC', 'Rate of evaporation of falling precip',     'lev',         'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('EVAPSNOW', 'Rate of evaporation of falling snow',       'trop_cld_lev','avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('HPROGCLD', 'Heating from prognostic clouds',            'lev',         'avg', 'W kg-1'   ) !subcyc
   call history_add_field ('FICE',     'Fractional ice content within cloud',       'lev',         'avg', 'fraction') !subcyc
   call history_add_field ('CLDFSNOW', 'Cloud fraction adjusted for snow',          'lev',         'avg', '1' ) !subcyc
   call history_add_field ('ICWMRST',  'Prognostic in-stratus water mixing ratio',  'lev',         'avg', 'kg kg-1') !subcyc
   call history_add_field ('ICIMRST',  'Prognostic in-stratus ice mixing ratio',    'lev',         'avg', 'kg kg-1') !subcyc

   ! MG microphysics diagnostics
   call history_add_field ('QCSEVAP', 'Rate of evaporation of falling cloud water',  'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('QISEVAP', 'Rate of sublimation of falling cloud ice',    'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('QVRES',   'Rate of residual condensation term',          'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('CMEIOUT', 'Rate of deposition/sublimation of cloud ice', 'trop_cld_lev', 'avg', 'kg kg-1/s') !subcyc
   call history_add_field ('VTRMC',   'Mass-weighted cloud water fallspeed',         'trop_cld_lev', 'avg', 'm s-1') !subcyc
   call history_add_field ('VTRMI',   'Mass-weighted cloud ice fallspeed',           'trop_cld_lev', 'avg', 'm s-1') !subcyc

   call history_add_field ('QCSEDTEN','Cloud water mixing ratio tendency from sedimentation', 'trop_cld_lev', 'avg', &
                                      'kg kg-1 s-1') !subcyc
   call history_add_field ('QISEDTEN','Cloud ice mixing ratio tendency from sedimentation',   'trop_cld_lev', 'avg', &
                                      'kg kg-1 s-1') !subcyc

   call history_add_field ('PRAO',    'Accretion of cloud water by rain',   'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('PRCO',    'Autoconversion of cloud water',      'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MNUCCCO', 'Immersion freezing of cloud water',  'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MNUCCTO', 'Contact freezing of cloud water',    'lev', 'avg', 'kg kg-1 s-1') !subcyc

   call history_add_field ('MNUCCDO', 'Homogeneous and heterogeneous nucleation from vapor',  'trop_cld_lev', 'avg', &
                           'kg kg-1 s-1') !subcyc

   call history_add_field ('MNUCCDOhet', 'Heterogeneous nucleation from vapor',                 'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MSACWIO',    'Conversion of cloud water from rime-splintering',     'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('PSACWSO',    'Accretion of cloud water by snow',                    'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('BERGSO',     'Conversion of cloud water to snow from bergeron',     'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('BERGO',      'Conversion of cloud water to cloud ice from bergeron','lev', 'avg', 'kg kg-1 s-1') !subcyc

   call history_add_field ('MELTO',    'Melting of cloud ice',                       'lev',          'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MELTSTOT', 'Melting of snow',                            'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MNUDEPO',  'Deposition Nucleation',                      'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('HOMOO',    'Homogeneous freezing of cloud water',        'lev',          'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('QCRESO',   'Residual condensation term for cloud water', 'lev',          'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('PRCIO',    'Autoconversion of cloud ice to snow',        'lev',          'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('PRAIO',    'Accretion of cloud ice to snow',             'lev',          'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('QIRESO',   'Residual deposition term for cloud ice',     'lev',          'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MNUCCRO',  'Heterogeneous freezing of rain to snow',     'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MNUCCRIO', 'Heterogeneous freezing of rain to ice',      'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('PRACSO',   'Accretion of rain by snow',                  'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('VAPDEPSO', 'Vapor deposition onto snow',                 'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MELTSDT',  'Latent heating rate due to melting of snow', 'trop_cld_lev', 'avg', 'W kg-1') !subcyc

   call history_add_field ('FRZRDT',   'Latent heating rate due to homogeneous freezing of rain',       'trop_cld_lev', 'avg', &
                           'W kg-1') !subcyc
   call history_add_field ('QRSEDTEN', 'Rain mixing ratio tendency from sedimentation',                 'trop_cld_lev', 'avg', &
                           'kg kg-1 s-1') !subcyc
   call history_add_field ('QSSEDTEN', 'Snow mixing ratio tendency from sedimentation',                 'trop_cld_lev', 'avg', &
                           'kg kg-1 s-1')  !subcyc
   call history_add_field ('NNUCCCO',  'Number Tendency due to Immersion freezing of cloud water',      'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NNUCCTO',  'Number Tendency due to Contact freezing of cloud water',        'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NNUCCDO',  'Number Tendency due to Ice nucleation',                         'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NNUDEPO',  'Number Tendency due to Deposition Nucleation',                  'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NHOMO',    'Number Tendency due to Homogeneous freezing of cloud water',    'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NNUCCRO',  'Number Tendency due to heterogeneous freezing of rain to snow', 'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NNUCCRIO', 'Number Tendency due to Heterogeneous freezing of rain to ice',  'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NSACWIO',  'Number Tendency due to Ice Multiplication- Rime-splintering',   'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NPRAO',    'Number Tendency due to Accretion of cloud water by rain',       'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NPSACWSO', 'Number Tendency due to Accretion of cloud water by snow',       'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NPRAIO',   'Number Tendency due to Accretion of cloud ice to snow',         'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NPRACSO',  'Number Tendency due to Accretion of rain by snow',              'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NPRCO',    'Number Tendency due to Autoconversion of cloud water [to rain]','trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NPRCIO',   'Number Tendency due to Autoconversion of cloud ice to snow',    'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NCSEDTEN', 'Number Tendency due to cloud liquid sedimentation',             'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NISEDTEN', 'Number Tendency due to cloud ice sedimentation',                'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc

   call history_add_field ('NRSEDTEN', 'Number Tendency due to rain sedimentation',   'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
   call history_add_field ('NSSEDTEN', 'Number Tendency due to snow sedimentation',   'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
   call history_add_field ('NMELTO',   'Number Tendency due to Melting of cloud ice', 'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
   call history_add_field ('NMELTS',   'Number Tendency due to Melting of snow',      'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc

! CACNOTE
!   if (trim(micro_mg_warm_rain) == 'kk2000') then
      call history_add_field ('qctend_KK2000',   'cloud liquid mass tendency due to autoconversion accretion from KK2000', &
                              'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('nctend_KK2000',   'cloud number mass tendency due to autoconversion accretion from KK2000', &
                              'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
      call history_add_field ('qrtend_KK2000',   'rain mass tendency due to autoconversion accretion from KK2000',         &
                              'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('nrtend_KK2000',   'rain number tendency due to autoconversion accretion from KK2000',       &
                              'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
!   end if

! CACNOTE
!   if (trim(micro_mg_warm_rain) == 'sb2001') then
      call history_add_field ('qctend_SB2001',  'cloud liquid mass tendency due to autoconversion  accretion from SB2001',  &
                              'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('nctend_SB2001',  'cloud liquid number tendency due to autoconversion accretion from SB2001', &
                              'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
      call history_add_field ('qrtend_SB2001',  'rain mass tendency due to autoconversion accretion from SB2001',           &
                              'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('nrtend_SB2001',  'rain number tendency due to autoconversion accretion from SB2001',         &
                              'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
!   end if

   call history_add_field ('LAMC', 'Size distribution parameter lambda for liquid',    'trop_cld_lev', 'avg', 'unitless') !subcyc
   call history_add_field ('LAMR', 'Size distribution parameter lambda for rain',      'trop_cld_lev', 'avg', 'unitless') !subcyc
   call history_add_field ('PGAM', 'Size distribution parameter mu (pgam) for liquid', 'trop_cld_lev', 'avg', 'unitless') !subcyc
   call history_add_field ('N0R',  'Size distribution parameter n0 for rain',          'trop_cld_lev', 'avg', 'unitless') !subcyc

! CACNOTE
!   if (micro_mg_version > 2) then
         call history_add_field ('NMELTG',  'Number Tendency due to Melting of graupel',                  'trop_cld_lev', 'avg', &
                                 '# kg-1 s-1') !subcyc
         call history_add_field ('NGSEDTEN', 'Number Tendency due to graupel sedimentation',              'trop_cld_lev', 'avg', &
                                 '# kg-1 s-1') !subcyc
         call history_add_field ('PSACRO', 'Collisions between rain & snow (Graupel collecting snow)',    'lev',          'avg', &
                                 'kg kg-1 s-1') !subcyc
         call history_add_field ('PRACGO', 'Change in q collection rain by graupel',                      'lev',          'avg', &
                                 'kg kg-1 s-1') !subcyc
         call history_add_field ('PSACWGO',  'Change in q collection droplets by graupel',                'lev',          'avg', &
                                 'kg kg-1 s-1') !subcyc
         call history_add_field ('PGSACWO',  'Q conversion to graupel due to collection droplets by snow','lev',          'avg', &
                                 'kg kg-1 s-1') !subcyc
         call history_add_field ('PGRACSO',  'Q conversion to graupel due to collection rain by snow',    'lev',          'avg', &
                                 'kg kg-1 s-1') !subcyc

         call history_add_field ('PRDGO',     'Deposition of graupel',                                            'lev', &
                                 'avg', 'kg kg-1 s-1') !subcyc
         call history_add_field ('QMULTGO',  'Q change due to ice mult droplets/graupel',                         'lev', &
                                 'avg', 'kg kg-1 s-1') !subcyc
         call history_add_field ('QMULTRGO', 'Q change due to ice mult rain/graupel',                             'lev', &
                                 'avg', 'kg kg-1 s-1') !subcyc
         call history_add_field ('QGSEDTEN',  'Graupel/Hail mixing ratio tendency from sedimentation',            'trop_cld_lev', &
                                 'avg', 'kg kg-1 s-1') !subcyc
         call history_add_field ('NPRACGO',   'Change N collection rain by graupel',                              'lev', &
                                 'avg', '# kg-1 s-1') !subcyc
         call history_add_field ('NSCNGO',    'Change N conversion to graupel due to collection droplets by snow','lev', &
                                 'avg', '# kg-1 s-1' ) !subcyc
         call history_add_field ('NGRACSO',    'Change N conversion to graupel due to collection rain by snow',   'lev', &
                                 'avg', '# kg-1 s-1') !subcyc
         call history_add_field ('NMULTGO',    'Ice mult due to acc droplets by graupel',                         'lev', &
                                 'avg', '# kg-1 s-1') !subcyc
         call history_add_field ('NMULTRGO',  'Ice mult due to acc rain by graupel',                              'lev', &
                                 'avg', '# kg-1 s-1') !subcyc
         call history_add_field ('NPSACWGO',  'Change N collection droplets by graupel',                          'lev', &
                                 'avg', '# kg-1 s-1') !subcyc
         call history_add_field ('CLDFGRAU',  'Cloud fraction adjusted for graupel',                              'lev', &
                                 'avg', '1') !subcyc
         call history_add_field ('MELTGTOT',  'Melting of graupel',                                               'trop_cld_lev', &
                                 'avg', 'kg kg-1 s-1') !subcyc
!   end if


   call history_add_field ('RBFRAC',  'Fraction of sky covered by a potential rainbow', horiz_only, 'avg',  'Fraction')  !subcyc
   call history_add_field ('RBFREQ',  'Potential rainbow frequency',                    horiz_only, 'avg',  'Frequency') !subcyc
   call history_add_field( 'rbSZA', 'solar zenith angle',                               horiz_only, 'inst', 'degrees')   !subcyc

   ! History variables for CAM5 microphysics
   call history_add_field ('MPDT',      'Heating tendency - Morrison microphysics',             'lev', 'avg', 'W kg-1') !subcyc
   call history_add_field ('MPDQ',      'Q tendency - Morrison microphysics',                   'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDLIQ',    'CLDLIQ tendency - Morrison microphysics',              'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDICE',    'CLDICE tendency - Morrison microphysics',              'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDNLIQ',   'NUMLIQ tendency - Morrison microphysics',              'lev', 'avg', 'kg-1 s-1') !subcyc
   call history_add_field ('MPDNICE',   'NUMICE tendency - Morrison microphysics',              'lev', 'avg', 'kg-1 s-1') !subcyc
   call history_add_field ('MPDW2V',    'Water <--> Vapor tendency - Morrison microphysics',    'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDW2I',    'Water <--> Ice tendency - Morrison microphysics',      'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDW2P',    'Water <--> Precip tendency - Morrison microphysics',   'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDI2V',    'Ice <--> Vapor tendency - Morrison microphysics',      'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDI2W',    'Ice <--> Water tendency - Morrison microphysics',      'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDI2P',    'Ice <--> Precip tendency - Morrison microphysics',     'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('ICWNC',     'Prognostic in-cloud water number conc',                'lev', 'avg', 'm-3') !subcyc
   call history_add_field ('ICINC',     'Prognostic in-cloud ice number conc',                  'lev', 'avg', 'm-3') !subcyc
   call history_add_field ('EFFLIQ_IND','Prognostic droplet effective radius (indirect effect)','lev', 'avg','Micron') !subcyc

   call history_add_field ('CDNUMC',    'Vertically-integrated droplet concentration',                     horiz_only,  'avg', &
                           '1 m-2') !subcyc
   call history_add_field ('MPICLWPI',  'Vertically-integrated in-cloud Initial Liquid WP (Before Micro)', horiz_only,  'avg', &
                           'kg m-2') !subcyc
   call history_add_field ('MPICIWPI', 'Vertically-integrated in-cloud Initial Ice WP (Before Micro)',     horiz_only,  'avg', &
                           'kg m-2') !subcyc


   ! This is only if the coldpoint temperatures are being adjusted.
   ! NOTE: Some fields related to these and output later are added in tropopause.F90.
!CACNOTE
!   if (micro_mg_adjust_cpt) then
     call history_add_field ('TROPF_TADJ', 'Temperatures after cold point adjustment',      'lev',        'avg', 'K') !subcyc
     call history_add_field ('TROPF_RHADJ','Relative Hunidity after cold point adjustment', 'lev',        'avg', 'K') !subcyc
     call history_add_field ('TROPF_CDT',  'Cold point temperature adjustment',              horiz_only,  'avg', 'K') !subcyc
     call history_add_field ('TROPF_CDZ',  'Distance of coldpoint from coldest model level', horiz_only,  'avg', 'm') !subcyc
!   end if


   ! Averaging for cloud particle number and size
   call history_add_field ('AWNC', 'Average cloud water number conc',     'lev',  'avg', 'm-3') !subcyc
   call history_add_field ('AWNI',  'Average cloud ice number conc',      'lev',  'avg', 'm-3') !subcyc
   call history_add_field ('AREL',   'Average droplet effective radius',  'lev',  'avg', 'Micron') !subcyc
   call history_add_field ('AREI',   'Average ice effective radius',      'lev',  'avg', 'Micron') !subcyc
   ! Frequency arrays for above
   call history_add_field ('FREQL', 'Fractional occurrence of liquid',    'lev',  'avg', 'fraction') !subcyc
   call history_add_field ('FREQI', 'Fractional occurrence of ice',       'lev',  'avg', 'fraction') !subcyc


   ! Average cloud top particle size and number (liq, ice) and frequency
   call history_add_field ('ACTREL', 'Average Cloud Top droplet effective radius', horiz_only,   'avg', 'Micron') !subcyc
   call history_add_field ('ACTREI', 'Average Cloud Top ice effective radius',     horiz_only,   'avg', 'Micron') !subcyc
   call history_add_field ('ACTNL',  'Average Cloud Top droplet number',           horiz_only,   'avg', 'm-3') !subcyc
   call history_add_field ('ACTNI',  'Average Cloud Top ice number',               horiz_only,   'avg', 'm-3') !subcyc

   call history_add_field ('FCTL',   'Fractional occurrence of cloud top liquid',  horiz_only,   'avg', 'fraction') !subcyc
   call history_add_field ('FCTI',   'Fractional occurrence of cloud top ice',     horiz_only,   'avg', 'fraction') !subcyc

   ! New frequency arrays for mixed phase and supercooled liquid (only and mixed) for (a) Cloud Top and (b) everywhere..
   call history_add_field ('FREQM',   'Fractional occurrence of mixed phase',                 'lev',  'avg', 'fraction') !subcyc
   call history_add_field ('FREQSL',  'Fractional occurrence of only supercooled liquid',     'lev',  'avg', 'fraction') !subcyc
   call history_add_field ('FREQSLM', 'Fractional occurrence of super cooled liquid with ice','lev',  'avg', 'fraction') !subcyc

   call history_add_field ('FCTM',    'Fractional occurrence of cloud top mixed phase',                 horiz_only,   'avg', &
                           'fraction') !subcyc
   call history_add_field ('FCTSL',   'Fractional occurrence of cloud top only supercooled liquid',     horiz_only,   'avg', &
                           'fraction') !subcyc
   call history_add_field ('FCTSLM',  'Fractional occurrence of cloud top super cooled liquid with ice',horiz_only,   'avg', &
                           'fraction') !subcyc

   call history_add_field ('LS_FLXPRC',  'ls stratiform gbm interface rain+snow flux',  'ilev', 'avg', 'kg m-2 s-1') !subcyc
   call history_add_field ('LS_FLXSNW',  'ls stratiform gbm interface snow flux',       'ilev', 'avg', 'kg m-2 s-1') !subcyc

   call history_add_field ('REL',  'MG REL stratiform cloud effective radius liquid', 'lev',  'avg', 'micron') !subcyc
   call history_add_field ('REI',  'MG REI stratiform cloud effective radius ice',    'lev',  'avg', 'micron') !subcyc
   call history_add_field ('LS_REFFRAIN',  'ls stratiform rain effective radius',     'lev',  'avg', 'micron') !subcyc
   call history_add_field ('LS_REFFSNOW','ls stratiform snow effective radius',       'lev',  'avg', 'micron') !subcyc
   call history_add_field ('CV_REFFLIQ',  'convective cloud liq effective radius',    'lev',  'avg', 'micron') !subcyc
   call history_add_field ('CV_REFFICE',  'convective cloud ice effective radius',    'lev',  'avg', 'micron') !subcyc
   call history_add_field ('MG_SADICE',  'MG surface area density ice',               'lev',  'avg', 'cm2 cm-3') !subcyc
   call history_add_field ('MG_SADSNOW', 'MG surface area density snow',              'lev',  'avg', 'cm2 cm-3') !subcyc

   ! diagnostic precip
   call history_add_field ('QRAIN', 'Diagnostic grid-mean rain mixing ratio',  'lev',  'avg', 'kg kg-1') !subcyc
   call history_add_field ('QSNOW',  'Diagnostic grid-mean snow mixing ratio', 'lev',  'avg', 'kg kg-1') !subcyc
   call history_add_field ('NRAIN', 'Diagnostic grid-mean rain number conc',   'lev',  'avg', 'm-3') !subcyc
   call history_add_field ('NSNOW', 'Diagnostic grid-mean snow number conc',   'lev',  'avg', 'm-3') !subcyc

   ! size of precip
   call history_add_field ('RERCLD',  'Diagnostic effective radius of Liquid Cloud and Rain',    'lev',  'avg', 'm') !subcyc
   call history_add_field ('DSNOW', 'Diagnostic grid-mean snow diameter',                        'lev',  'avg', 'm') !subcyc

   ! diagnostic radar reflectivity, cloud-averaged
   call history_add_field ('REFL',  '94 GHz radar reflectivity',                   'lev',  'avg', 'DBz') !subcyc
   call history_add_field ('AREFL', 'Average 94 GHz radar reflectivity',           'lev',  'avg', 'DBz') !subcyc
   call history_add_field ('FREFL', 'Fractional occurrence of radar reflectivity', 'lev',  'avg', 'fraction') !subcyc

   call history_add_field ('CSRFL', '94 GHz radar reflectivity (CloudSat thresholds)',           'lev',  'avg', 'DBz') !subcyc
   call history_add_field ('ACSRFL', 'Average 94 GHz radar reflectivity (CloudSat thresholds)',  'lev',  'avg', 'DBz') !subcyc

   call history_add_field ('FCSRFL', 'Fractional occurrence of radar reflectivity (CloudSat thresholds)',  'lev',  'avg', &
                           'fraction') !subcyc

   call history_add_field ('AREFLZ', 'Average 94 GHz radar reflectivity',           'lev',  'avg', 'mm6 m-3') !subcyc

   ! 10cm (rain) radar reflectivity
   call history_add_field ('REFL10CM',  '10cm (Rain) radar reflectivity (Dbz)',     'lev',  'avg', 'DBz') !subcyc
   call history_add_field ('REFLZ10CM', '10cm (Rain) radar reflectivity (Z units)', 'lev',  'avg', 'mm6 m-3') !subcyc

   ! Aerosol information
   call history_add_field ('NCAL',  'Number Concentation Activated for Liquid',   'lev',  'avg', 'm-3') !subcyc
   call history_add_field ('NCAI',  'Number Concentation Activated for Ice',      'lev',  'avg', 'm-3') !subcyc

   ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
   call history_add_field ('AQRAIN',  'Average rain mixing ratio',        'lev',  'avg', 'kg kg-1') !subcyc
   call history_add_field ('AQSNOW',  'Average snow mixing ratio',        'lev',  'avg', 'kg kg-1') !subcyc
   call history_add_field ('ANRAIN',  'Average rain number conc',         'lev',  'avg', 'm-3') !subcyc
   call history_add_field ('ANSNOW',   'Average snow number conc',        'lev',  'avg', 'm-3') !subcyc
   call history_add_field ('ADRAIN',  'Average rain effective Diameter',  'lev',  'avg', 'm') !subcyc
   call history_add_field ('ADSNOW',   'Average snow effective Diameter', 'lev',  'avg', 'm') !subcyc
   call history_add_field ('FREQR',    'Fractional occurrence of rain',   'lev',  'avg', 'fraction') !subcyc
   call history_add_field ('FREQS',  'Fractional occurrence of snow',     'lev',  'avg', 'fraction') !subcyc

   ! precipitation efficiency & other diagnostic fields
   call history_add_field('PE'    ,  'Stratiform Precipitation Efficiency  (precip/cmeliq)',          horiz_only,   'avg', &
                          '1') !subcyc
   call history_add_field('APRL'  ,   'Average Stratiform Precip Rate over efficiency calculation',   horiz_only,   'avg', &
                          'm s-1') !subcyc

   call history_add_field('PEFRAC', 'Fraction of timesteps precip efficiency reported',    horiz_only, 'avg', '1') !subcyc
   call history_add_field('VPRCO' , 'Vertical average of autoconversion rate',             horiz_only, 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field('VPRAO' , 'Vertical average of accretion rate',                  horiz_only, 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field('RACAU' , 'Accretion/autoconversion ratio from vertical average',horiz_only, 'avg', 'kg kg-1 s-1') !subcyc

   call history_add_field('UMR','Mass-weighted rain  fallspeed', 'trop_cld_lev', 'avg',   'm s-1') !subcyc

! CACNOTE
!   if (micro_mg_version > 2) then
      call history_add_field('UMG',        'Mass-weighted graupel/hail  fallspeed',       'trop_cld_lev', 'avg', 'm s-1') !subcyc
      call history_add_field('FREQG',      'Fractional occurrence of Graupel',            'lev',          'avg', 'fraction') !subcyc
      call history_add_field('LS_REFFGRAU','ls stratiform graupel/hail effective radius', 'lev',          'avg', 'micron') !subcyc
      call history_add_field('AQGRAU',     'Average graupel/hail mixing ratio',           'lev',          'avg', 'kg kg-1') !subcyc
      call history_add_field('ANGRAU',     'Average graupel/hail number conc',            'lev',          'avg', 'm-3') !subcyc
!   end if


   ! qc limiter (only output in versions 1.5 and later)
   call history_add_field('QCRAT', 'Qc Limiter: Fraction of qc tendency applied', 'lev', 'avg', 'fraction') !subcyc

end subroutine pumas_diagnostics_init


end module pumas_diagnostics
