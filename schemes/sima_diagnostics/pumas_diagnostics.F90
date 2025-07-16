!CACNOTE - All outfld calls to history_out_field are converted
         - All history_add_field calls are done
               - Except the constituent history_add_field needs to be changed - see CACNOTE
               - Need to replace "#" in units?
         - need to line up calls
!!!----------------------------------------------

module pumas_diagnostics

   use ccpp_kinds, only:  kind_phys

   implicit none
   private
   save

   public :: pumas_diagnostics_init ! init routine
   public :: pumas_diagnostics_run  ! main routine

CONTAINS

!> \section arg_table_pumas_diagnostics_init  Argument Table
!! \htmlinclude pumas_diagnostics_init.html
subroutine pumas_diagnostics_init(errmsg, errflg)

   use cam_history,         only: history_add_field
   use cam_history_support, only: horiz_only

   character(len=512), intent(out) :: errmsg
   integer,            intent(out) :: errflg

   ! Local variables:

   errmsg = ''
   errflg = 0

   ! NOTE -subcyc comments below mean that the original CAM code had "sampled_on_subcycle=.true." in the old outfld call
   !        Will use this for guidance during future developmen

!CACNOTE - come back to this - did not convert get the constituent info yet

   do m = 1, ncnst
      call cnst_get_ind(cnst_names(m), mm)
      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixsnow, ixgraupel /)) ) then
         ! mass mixing ratios
         call history_add_field(cnst_name(mm), cnst_longname(mm),                    'lev',      'avg', 'kg kg-1') !subcyc
         call history_add_field(sflxnam(mm),   trim(cnst_name(mm))//' surface flux', horiz_only, 'avg', 'kg m-2 s-1') !subcyc
      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain, ixnumsnow, ixnumgraupel /)) ) then
         ! number concentrations
         call history_add_field(cnst_name(mm), cnst_longname(mm),                    'lev',      'avg', 'kg-1') !subcyc
         call history_add_field(sflxnam(mm),  trim(cnst_name(mm))//' surface flux',  horiz_only, 'avg', '1 m-2 s-1') !subcyc
      else
         call endrun( "micro_pumas_cam_init: &
              &Could not call history_add_field for constituent with unknown units.")
      endif
   end do

   call history_add_field(apcnst(ixcldliq), trim(cnst_name(ixcldliq))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
   call history_add_field(apcnst(ixcldice), trim(cnst_name(ixcldice))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
   call history_add_field(bpcnst(ixcldliq), trim(cnst_name(ixcldliq))//' before physics','lev', 'avg', 'kg kg-1') !subcyc
   call history_add_field(bpcnst(ixcldice), trim(cnst_name(ixcldice))//' before physics','lev', 'avg', 'kg kg-1') !subcyc

   call history_add_field(apcnst(ixrain), trim(cnst_name(ixrain))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
   call history_add_field(apcnst(ixsnow), trim(cnst_name(ixsnow))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
   call history_add_field(bpcnst(ixrain), trim(cnst_name(ixrain))//' before physics','lev', 'avg', 'kg kg-1') !subcyc
   call history_add_field(bpcnst(ixsnow), trim(cnst_name(ixsnow))//' before physics','lev', 'avg', 'kg kg-1') !subcyc

   if (micro_mg_version > 2) then
      call history_add_field(apcnst(ixgraupel), trim(cnst_name(ixgraupel))//' after physics', 'lev', 'avg', 'kg kg-1') !subcyc
      call history_add_field(bpcnst(ixgraupel), trim(cnst_name(ixgraupel))//' before physics','lev', 'avg', 'kg kg-1') !subcyc
   end if

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
                           'kg kg-1 s-1''  !subcyc
   call history_add_field ('NNUCCCO',  'Number Tendency due to Immersion freezing of cloud water',      'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NNUCCTO',  'Number Tendency due to Contact freezing of cloud water',        'trop_cld_lev', 'avg', &
                           '# kg-1 s-1') !subcyc
   call history_add_field ('NNUCCDO',  'Number Tendency due to Ice nucleation',                         'trop_cld_lev', 'avg', &
                           '# kg-1 s-1'') !subcyc
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

   if (trim(micro_mg_warm_rain) == 'kk2000') then
      call history_add_field ('qctend_KK2000',   'cloud liquid mass tendency due to autoconversion accretion from KK2000', &
                              'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('nctend_KK2000',   'cloud number mass tendency due to autoconversion accretion from KK2000', &
                              'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
      call history_add_field ('qrtend_KK2000',   'rain mass tendency due to autoconversion accretion from KK2000',         &
                              'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('nrtend_KK2000',   'rain number tendency due to autoconversion accretion from KK2000',       &
                              'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
   end if

   if (trim(micro_mg_warm_rain) == 'sb2001') then
      call history_add_field ('qctend_SB2001',  'cloud liquid mass tendency due to autoconversion  accretion from SB2001',  &
                              'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('nctend_SB2001',  'cloud liquid number tendency due to autoconversion accretion from SB2001', &
                              'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
      call history_add_field ('qrtend_SB2001',  'rain mass tendency due to autoconversion accretion from SB2001'),          &
                              'trop_cld_lev', 'avg', 'kg kg-1 s-1') !subcyc
      call history_add_field ('nrtend_SB2001',  'rain number tendency due to autoconversion accretion from SB2001',         &
                              'trop_cld_lev', 'avg', '# kg-1 s-1') !subcyc
   end if

   call history_add_field ('LAMC', 'Size distribution parameter lambda for liquid',    'trop_cld_lev', 'avg', 'unitless') !subcyc
   call history_add_field ('LAMR', 'Size distribution parameter lambda for rain',      'trop_cld_lev', 'avg', 'unitless') !subcyc
   call history_add_field ('PGAM', 'Size distribution parameter mu (pgam) for liquid', 'trop_cld_lev', 'avg', 'unitless') !subcyc
   call history_add_field ('N0R',  'Size distribution parameter n0 for rain',          'trop_cld_lev', 'avg', 'unitless') !subcyc

   if (micro_mg_version > 2) then
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
   end if


   call history_add_field ('RBFRAC',  'Fraction of sky covered by a potential rainbow', horiz_only, 'avg',  'Fraction')  !subcyc
   call history_add_field ('RBFREQ',  'Potential rainbow frequency',                    horiz_only, 'avg',  'Frequency') !subcyc
   call history_add_field( 'rbSZA', 'solar zenith angle',                               horiz_only, 'inst', 'degrees')   !subcyc

   ! History variables for CAM5 microphysics
   call history_add_field ('MPDT',      'Heating tendency - Morrison microphysics',             'lev', 'avg', 'W kg-1') !subcyc
   call history_add_field ('MPDQ',      'Q tendency - Morrison microphysics',                   'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDLIQ',    'CLDLIQ tendency - Morrison microphysics',              'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDICE',    'CLDICE tendency - Morrison microphysics',              'lev', 'avg', 'kg kg-1 s-1') !subcyc
   call history_add_field ('MPDNLIQ',   'NUMLIQ tendency - Morrison microphysics',              'lev', 'avg', 'kg-1 s-1') !subcyc
   call history_add_field ('MPDNICE',   'NUMICE tendency - Morrison microphysics'),             'lev', 'avg', 'kg-1 s-1') !subcyc
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
   if (micro_mg_adjust_cpt) then
     call history_add_field ('TROPF_TADJ', 'Temperatures after cold point adjustment',      'lev',        'avg', 'K') !subcyc
     call history_add_field ('TROPF_RHADJ','Relative Hunidity after cold point adjustment', 'lev',        'avg', 'K') !subcyc
     call history_add_field ('TROPF_CDT',  'Cold point temperature adjustment',              horiz_only,  'avg', 'K') !subcyc
     call history_add_field ('TROPF_CDZ',  'Distance of coldpoint from coldest model level', horiz_only,  'avg', 'm') !subcyc
   end if


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
   call history_add_field ('CV_REFFICE',   'convective cloud ice effective radius, ,  'lev'   'avg', 'micron') !subcyc
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

   if (micro_mg_version > 2) then
      call history_add_field('UMG',        'Mass-weighted graupel/hail  fallspeed',       'trop_cld_lev', 'avg', 'm s-1') !subcyc
      call history_add_field('FREQG',      'Fractional occurrence of Graupel',            'lev',          'avg', 'fraction') !subcyc
      call history_add_field('LS_REFFGRAU','ls stratiform graupel/hail effective radius', 'lev',          'avg', 'micron') !subcyc
      call history_add_field('AQGRAU',     'Average graupel/hail mixing ratio',           'lev',          'avg', 'kg kg-1') !subcyc
      call history_add_field('ANGRAU',     'Average graupel/hail number conc',            'lev',          'avg', 'm-3') !subcyc
   end if


   ! qc limiter (only output in versions 1.5 and later)
   call history_add_field('QCRAT', 'Qc Limiter: Fraction of qc tendency applied', 'lev', 'avg', 'fraction') !subcyc

end subroutine pumas_diagnostics_init

!===============================================================================

subroutine pumas_diagnostics_run(state, ptend, dtime, pbuf)


!> \section arg_table_pumas_diagnostics_run  Argument Table
!! \htmlinclude pumas_diagnostics_run.html
subroutine pumas_diagnostics_run(proc_rates, qcsinksum_rate1ord, naai, npccn, rndst, nacon, tlat, qvlat, qctend, qitend, &
                   nctend, nitend, qrtend, qstend, nrtend, nstend, qgtend, ngtend, effc, effc_fn, effi, sadice, sadsnow, &
                   prect, preci, nevapr, am_evp_st, prain, cmeout, deffi, pgamrad, lamcrad, qsout, dsout, qgout, ngout, &
                   dgout, lflx, iflx, gflx, rflx, sflx, qrout, reff_rain, reff_snow, reff_grau, nrout, nsout, refl, &
                   arefl, areflz, frefl, csrfl, acsrfl, fcsrfl, refl10cm,  reflz10cm, rercld, ncai, ncal, qrout2, qsout2, &
                   nrout2, nsout2, drout2, dsout2, qgout2, ngout2, dgout2, freqg, freqs, freqr, nfice, qcrat, proc_rates, &
                   errstring, tnd_qsnow, tnd_nsnow, re_ice, prer_evap, frzimm, frzcnt, frzdep,  errmsg,  errflg)

   use cam_history, only: history_in_field

   type (proc_rates_type), intent(inout)  :: proc_rates

   real(kind_phys), intent(in) :: qcsinksum_rate1ord(:,:) ! 1st order rate for direct cw to precip conversion
   real(kind_phys), intent(in) :: tlat(:,:)               ! latent heating rate       (W/kg)
   real(kind_phys), intent(in) :: qvlat(:,:)              ! microphysical tendency qv (1/s)
   real(kind_phys), intent(in) :: qctend(:,:)             ! microphysical tendency qc (1/s)
   real(kind_phys), intent(in) :: qitend(:,:)             ! microphysical tendency qi (1/s)
   real(kind_phys), intent(in) :: nctend(:,:)             ! microphysical tendency nc (1/(kg*s))
   real(kind_phys), intent(in) :: nitend(:,:)             ! microphysical tendency ni (1/(kg*s))

   real(kind_phys), intent(in) :: qrtend(:,:)             ! microphysical tendency qr (1/s)
   real(kind_phys), intent(in) :: qstend(:,:)             ! microphysical tendency qs (1/s)
   real(kind_phys), intent(in) :: nrtend(:,:)             ! microphysical tendency nr (1/(kg*s))
   real(kind_phys), intent(in) :: nstend(:,:)             ! microphysical tendency ns (1/(kg*s))
   real(kind_phys), intent(in) :: qgtend(:,:)             ! microphysical tendency qg (1/s)
   real(kind_phys), intent(in) :: ngtend(:,:)             ! microphysical tendency ng (1/(kg*s))

   real(kind_phys), intent(in) :: effc(:,:)               ! droplet effective radius (micron)
   real(kind_phys), intent(in) :: effc_fn(:,:)            ! droplet effective radius, assuming nc = 1.e8 kg-1
   real(kind_phys), intent(in) :: effi(:,:)               ! cloud ice effective radius (micron)
   real(kind_phys), intent(in) :: sadice(:,:)             ! cloud ice surface area density (cm2/cm3)
   real(kind_phys), intent(in) :: sadsnow(:,:)            ! cloud snow surface area density (cm2/cm3)
   real(kind_phys), intent(in) :: prect(:)                ! surface precip rate (m/s)
   real(kind_phys), intent(in) :: preci(:)                ! surface ice/snow precip rate (m/s)
   real(kind_phys), intent(in) :: nevapr(:,:)             ! evaporation rate of rain + snow (kg/kg/s)
   real(kind_phys), intent(in) :: am_evp_st(:,:)          ! stratiform evaporation area (frac)
   real(kind_phys), intent(in) :: prain(:,:)              ! production of rain + snow (kg/kg/s)
   real(kind_phys), intent(in) :: cmeout(:,:)             ! Rate of cond-evap of ice (kg/kg/s)
   real(kind_phys), intent(in) :: deffi(:,:)              ! ice effective diameter for optics (radiation) (micron)
   real(kind_phys), intent(in) :: pgamrad(:,:)            ! ice gamma parameter for optics (radiation) (no units)
   real(kind_phys), intent(in) :: lamcrad(:,:)            ! slope of droplet distribution for optics (radiation) (1/m)
   real(kind_phys), intent(in) :: qsout(:,:)              ! snow mixing ratio (kg/kg)
   real(kind_phys), intent(in) :: dsout(:,:)              ! snow diameter (m)
   real(kind_phys), intent(in) :: lflx(:,:)               ! grid-box average liquid condensate flux (kg m^-2 s^-1)
   real(kind_phys), intent(in) :: iflx(:,:)               ! grid-box average ice condensate flux (kg m^-2 s^-1)
   real(kind_phys), intent(in) :: rflx(:,:)               ! grid-box average rain flux (kg m^-2 s^-1)
   real(kind_phys), intent(in) :: sflx(:,:)               ! grid-box average snow flux (kg m^-2 s^-1)
   real(kind_phys), intent(in) :: gflx(:,:)               ! grid-box average graupel/hail flux (kg m^-2 s^-1)

   real(kind_phys), intent(in) :: qrout(:,:)              ! grid-box average rain mixing ratio (kg/kg)
   real(kind_phys), intent(in) :: reff_rain(:,:)          ! rain effective radius (micron)
   real(kind_phys), intent(in) :: reff_snow(:,:)          ! snow effective radius (micron)
   real(kind_phys), intent(in) :: reff_grau(:,:)          ! graupel effective radius (micron)

   real(kind_phys), intent(in) :: nrout(:,:)              ! rain number concentration (1/m3)
   real(kind_phys), intent(in) :: nsout(:,:)              ! snow number concentration (1/m3)
   real(kind_phys), intent(in) :: refl(:,:)               ! analytic radar reflectivity (94GHZ, cloud radar) (dBZ)
   real(kind_phys), intent(in) :: arefl(:,:)              ! average reflectivity will zero points inside valid range (dBZ)
   real(kind_phys), intent(in) :: areflz(:,:)             ! average reflectivity in z. (mm6 m-3)
   real(kind_phys), intent(in) :: frefl(:,:)              ! fractional occurrence of radar reflectivity
   real(kind_phys), intent(in) :: csrfl(:,:)              ! cloudsat reflectivity (dBZ)
   real(kind_phys), intent(in) :: acsrfl(:,:)             ! cloudsat average (dBZ)
   real(kind_phys), intent(in) :: fcsrfl(:,:)             ! cloudsat fractional occurrence of radar reflectivity
   real(kind_phys), intent(in) :: refl10cm(:,:)           ! 10cm (rain) analytic radar reflectivity (dBZ)
   real(kind_phys), intent(in) :: reflz10cm(:,:)          ! 10cm (rain) analytic radar reflectivity (mm6 m-3)
   real(kind_phys), intent(in) :: rercld(:,:)             ! effective radius calculation for rain + cloud
   real(kind_phys), intent(in) :: ncai(:,:)               ! input number conc of ice nuclei available (1/m3)
   real(kind_phys), intent(in) :: ncal(:,:)               ! input number conc of CCN (1/m3)
   real(kind_phys), intent(in) :: qrout2(:,:)             ! copy of qrin as used to compute drin2
   real(kind_phys), intent(in) :: qsout2(:,:)             ! copy of qsin as used to compute dsin2
   real(kind_phys), intent(in) :: nrout2(:,:)             ! copy of nrin as used to compute drin2
   real(kind_phys), intent(in) :: nsout2(:,:)             ! copy of nsin as used to compute dsin2
   real(kind_phys), intent(in) :: drout2(:,:)             ! mean rain particle diameter (m)
   real(kind_phys), intent(in) :: dsout2(:,:)             ! mean snow particle diameter (m)
   real(kind_phys), intent(in) :: freqs(:,:)              ! fractional occurrence of snow
   real(kind_phys), intent(in) :: freqr(:,:)              ! fractional occurrence of rain
   real(kind_phys), intent(in) :: nfice(:,:)              ! fraction of frozen water to total condensed water
   real(kind_phys), intent(in) :: qcrat(:,:)              ! limiter for qc process rates (1=no limit --> 0. no qc)
   real(kind_phys), intent(in) :: qgout(:,:)              ! graupel/hail mixing ratio (kg/kg)
   real(kind_phys), intent(in) :: dgout(:,:)              ! graupel/hail diameter (m)
   real(kind_phys), intent(in) :: ngout(:,:)              ! graupel/hail number concentration (1/m3)
   real(kind_phys), intent(in) :: qgout2(:,:)             ! copy of qgin as used to compute dgin2
   real(kind_phys), intent(in) :: ngout2(:,:)             ! copy of ngin as used to compute dgin2
   real(kind_phys), intent(in) :: dgout2(:,:)             ! mean graupel/hail particle diameter (m)
   real(kind_phys), intent(in) :: freqg(:,:)              ! fractional occurrence of graupel
   real(kind_phys), intent(in) :: prer_evap(:,:)          ! evaporation rate of rain (kg/kg/s)

   ! CCPP error handling variables
   character(len=512), intent(in) :: errmsg
   integer,            intent(in) :: errflg

   errmsg = ''
   errflg = 0


! KATES HACKATHON NOTES
!! The calls for history_out_field are in the pumas_diagnostcis.F90 file, but also with all of the other pumas code.
!! Those calls need to be hooked up to the parameters in the call list here.
!! There is an issue with subcolumn/grid scale I'm not sure how to address the grid averaging that is done after the
!! pumas_tend function call.


   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   psetcols = state%psetcols
   ngrdcol  = state%ngrdcol
   itim_old = pbuf_old_tim_idx()
   nlev = pver - top_lev + 1

   nan_array = nan

   ! Allocate the proc_rates DDT
   ! IMPORTANT NOTE -- elements in proc_rates are dimensioned to the nlev dimension while
   !     all the other arrays in this routine are dimensioned pver.  This is required because
   !     PUMAS only gets the top_lev:pver array subsection, and the proc_rates arrays
   !     need to be the same levels.
   call proc_rates%allocate(ncol, nlev, ncd, micro_mg_warm_rain, errstring)

   call handle_errmsg(errstring, subname="micro_pumas_cam_tend")


   call phys_getopts(use_subcol_microp_out=use_subcol_microp)

   ! Set the col_type flag to grid or subcolumn dependent on the value of use_subcol_microp
   call pbuf_col_type_index(use_subcol_microp, col_type=col_type)

   !-----------------------
   ! These physics buffer fields are read only and not set in this parameterization
   ! If these fields do not have subcolumn data, copy the grid to the subcolumn if subcolumns is turned on
   ! If subcolumns is not turned on, then these fields will be grid data

   call pbuf_get_field(pbuf, naai_hom_idx,    naai_hom,    col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, relvar_idx,      relvar,      col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan, col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq,      col_type=col_type, copy_if_needed=use_subcol_microp)

   call pbuf_get_field(pbuf, cld_idx,         cld,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, concld_idx,      concld,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, ast_idx,         ast,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)

   ! Get convective precip for rainbows
   if (prec_dp_idx > 0) then
      call pbuf_get_field(pbuf, prec_dp_idx, prec_dp, col_type=col_type, copy_if_needed=use_subcol_microp)
   else
      nullify(prec_dp)
   end if
   if (prec_sh_idx > 0) then
      call pbuf_get_field(pbuf, prec_sh_idx, prec_sh, col_type=col_type, copy_if_needed=use_subcol_microp)
   else
      nullify(prec_sh)
   end if

! Merge Precipitation rates (multi-process)
   if (associated(prec_dp) .and. associated(prec_sh)) then
      precc(:ncol) = prec_dp(:ncol)  + prec_sh(:ncol)
   else if (associated(prec_dp)) then
      precc(:ncol) = prec_dp(:ncol)
   else if (associated(prec_sh)) then
      precc(:ncol) = prec_sh(:ncol)
   else
      precc(:ncol) = 0._r8
   end if

   if (do_cldice) then
      ! If we ARE prognosing tendencies, then just point to an array of NaN fields to have
      ! something for PUMAS to use in call
      tnd_qsnow => nan_array
      tnd_nsnow => nan_array
      re_ice => nan_array
   end if

   if (.not. use_hetfrz_classnuc) then
      ! Needed to satisfy gnu compiler with optional argument - set to an array of Nan fields
      frzimm => nan_array
      frzcnt => nan_array
      frzdep => nan_array
   end if

   if (qsatfac_idx > 0) then
      call pbuf_get_field(pbuf, qsatfac_idx, qsatfac, col_type=col_type, copy_if_needed=use_subcol_microp)
   else
      allocate(qsatfac(ncol,pver),stat=ierr)
      call handle_allocate_error(ierr, 'micro_pumas_cam_tend', 'qsatfac')
      qsatfac = 1._r8
   end if

   ! initialize tendency variables
    preci  = 0._r8
    prect  = 0._r8


   !-----------------------
   ! These physics buffer fields are calculated and set in this parameterization
   ! If subcolumns is turned on, then these fields will be calculated on a subcolumn grid, otherwise they will be a normal grid

   call pbuf_get_field(pbuf, prec_str_idx,    prec_str,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_str_idx,    snow_str,    col_type=col_type)
   call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw,    col_type=col_type)
   call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed,    col_type=col_type)
   call pbuf_get_field(pbuf, dei_idx,         dei,         col_type=col_type)
   call pbuf_get_field(pbuf, mu_idx,          mu,          col_type=col_type)
   call pbuf_get_field(pbuf, lambdac_idx,     lambdac,     col_type=col_type)
   call pbuf_get_field(pbuf, des_idx,         des,         col_type=col_type)
   call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc,    col_type=col_type)
   call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw,    col_type=col_type)
   call pbuf_get_field(pbuf, ls_mrprc_idx,    mgmrprc,     col_type=col_type)
   call pbuf_get_field(pbuf, ls_mrsnw_idx,    mgmrsnw,     col_type=col_type)
   call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq,   col_type=col_type)
   call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice,   col_type=col_type)
   call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst,     col_type=col_type)
   call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst,     col_type=col_type)
   call pbuf_get_field(pbuf, icswp_idx,       icswp,       col_type=col_type)
   call pbuf_get_field(pbuf, rel_idx,         rel,         col_type=col_type)
   call pbuf_get_field(pbuf, rei_idx,         rei,         col_type=col_type)
   call pbuf_get_field(pbuf, wsedl_idx,       wsedl,       col_type=col_type)
   call pbuf_get_field(pbuf, qme_idx,         qme,         col_type=col_type)
   call pbuf_get_field(pbuf, bergso_idx,      bergstot,    col_type=col_type)

   ! Assign the pointer values to the non-pointer proc_rates element
   proc_rates%bergstot(:ncol,1:nlev) = bergstot(:ncol,top_lev:pver)

   if (degrau_idx > 0)   call pbuf_get_field(pbuf, degrau_idx,   degrau,   col_type=col_type)
   if (icgrauwp_idx > 0) call pbuf_get_field(pbuf, icgrauwp_idx, icgrauwp, col_type=col_type)
   if (cldfgrau_idx > 0) call pbuf_get_field(pbuf, cldfgrau_idx, cldfgrau, col_type=col_type)

   call pbuf_get_field(pbuf, cldo_idx,        cldo,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cldfsnow_idx,    cldfsnow, start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_t_idx,        CC_t,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qv_idx,       CC_qv,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_ql_idx,       CC_ql,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qi_idx,       CC_qi,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_nl_idx,       CC_nl,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_ni_idx,       CC_ni,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qlst_idx,     CC_qlst,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)

   if (rate1_cw2pr_st_idx > 0) then
      call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st, col_type=col_type)
   end if

   !-----------------------
   ! These are only on the grid regardless of whether subcolumns are turned on or not
   call pbuf_get_field(pbuf, ls_reffrain_idx, mgreffrain_grid)
   call pbuf_get_field(pbuf, ls_reffsnow_idx, mgreffsnow_grid)
   call pbuf_get_field(pbuf, acpr_idx,        acprecl_grid)
   call pbuf_get_field(pbuf, acgcme_idx,      acgcme_grid)
   call pbuf_get_field(pbuf, acnum_idx,       acnum_grid)
   call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq_grid)
   call pbuf_get_field(pbuf, ast_idx,         ast_grid, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, evprain_st_idx,  evprain_st_grid)
   call pbuf_get_field(pbuf, evpsnow_st_idx,  evpsnow_st_grid)
!!!!! CACNOTE?   call pbuf_get_field(pbuf, am_evp_st_idx,   am_evp_st_grid)

   !-----------------------------------------------------------------------
   !        ... Calculate cosine of zenith angle
   !            then cast back to angle (radians)
   !-----------------------------------------------------------------------

   zen_angle(:) = 0.0_r8
   rlats(:) = 0.0_r8
   rlons(:) = 0.0_r8
   calday = get_curr_calday()
   call get_rlat_all_p( lchnk, ncol, rlats )
   call get_rlon_all_p( lchnk, ncol, rlons )
   call zenith( calday, rlats, rlons, zen_angle, ncol )
   where (zen_angle(:) <= 1.0_r8 .and. zen_angle(:) >= -1.0_r8)
      zen_angle(:) = acos( zen_angle(:) )
   elsewhere
      zen_angle(:) = 0.0_r8
   end where

   sza(:) = zen_angle(:) * rad2deg
   call history_out_field( 'rbSZA',   sza)

   !-------------------------------------------------------------------------------------
   ! Microphysics assumes 'liquid stratus frac = ice stratus frac
   !                      = max( liquid stratus frac, ice stratus frac )'.
   alst_mic => ast
   aist_mic => ast

   ! Output initial in-cloud LWP (before microphysics)

   iclwpi = 0._r8
   iciwpi = 0._r8

   do i = 1, ncol
      do k = top_lev, pver
         iclwpi(i) = iclwpi(i) + &
              min(state%q(i,k,ixcldliq) / max(mincld,ast(i,k)),0.005_r8) &
              * state%pdel(i,k) / gravit
         iciwpi(i) = iciwpi(i) + &
              min(state%q(i,k,ixcldice) / max(mincld,ast(i,k)),0.005_r8) &
              * state%pdel(i,k) / gravit
      end do
   end do

   cldo(:ncol,top_lev:pver)=ast(:ncol,top_lev:pver)

   ! Initialize local state from input.
   call physics_state_copy(state, state_loc)

   ! Because of the of limited vertical resolution, there can be a signifcant
   ! warm bias at the cold point tropopause, which can create a wet bias in the
   ! stratosphere. For the microphysics only, update the cold point temperature, with
   ! an estimate of the coldest point between the model layers.
   if (micro_mg_adjust_cpt) then
      cp_rh(:ncol, :pver)  = 0._r8
      cp_dt(:ncol)         = 0._r8
      cp_dz(:ncol)         = 0._r8

      !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
      troplev(:) = 0
      cp_z(:) = 0._r8
      cp_t(:) = 0._r8
      !REMOVECAM_END
      call tropopause_find_cam(state_loc, troplev, primary=TROP_ALG_CPP, backup=TROP_ALG_NONE, &
                           tropZ=cp_z, tropT=cp_t)

      do i = 1, ncol

         ! Update statistics and output results.
         if (troplev(i) .ne. NOTFOUND) then
            cp_dt(i) = cp_t(i) - state_loc%t(i,troplev(i))
            cp_dz(i) = cp_z(i) - state_loc%zm(i,troplev(i))

            ! NOTE: This change in temperature is just for the microphysics
            ! and should not be added to any tendencies or used to update
            ! any states
            state_loc%t(i,troplev(i)) = state_loc%t(i,troplev(i)) + cp_dt(i)
         end if
      end do

      ! Output all of the statistics related to the cold point
      ! tropopause adjustment. Th cold point information itself is
      ! output in tropopause.F90.
      call history_out_field("TROPF_TADJ", state_loc%t)
      call history_out_field("TROPF_CDT",  cp_dt)
      call history_out_field("TROPF_CDZ",  cp_dz)
   end if

   ! Initialize ptend for output.
   lq = .false.
   lq(ixq) = .true.
   lq(ixcldliq) = .true.
   lq(ixcldice) = .true.
   lq(ixnumliq) = .true.
   lq(ixnumice) = .true.
   lq(ixrain) = .true.
   lq(ixsnow) = .true.
   lq(ixnumrain) = .true.
   lq(ixnumsnow) = .true.
   if (micro_mg_version > 2) then
      lq(ixgraupel) = .true.
      lq(ixnumgraupel) = .true.
   end if

   ! the name 'cldwat' triggers special tests on cldliq
   ! and cldice in physics_update
   call physics_ptend_init(ptend, psetcols, "cldwat", ls=.true., lq=lq)

   if (micro_mg_version > 2) then
      state_loc_graup(:ncol,:) = state_loc%q(:ncol,:,ixgraupel)
      state_loc_numgraup(:ncol,:) = state_loc%q(:ncol,:,ixnumgraupel)
   else
      state_loc_graup(:ncol,:) = 0._r8
      state_loc_numgraup(:ncol,:) = 0._r8
   end if

   ! Zero out diagnostic rainbow arrays
   rbfreq = 0._r8
   rbfrac = 0._r8

   ! Zero out values above top_lev before passing into _tend for some pbuf variables that are inputs
   naai(:ncol,:top_lev-1) = 0._r8
   npccn(:ncol,:top_lev-1) = 0._r8

   ! The null value for qsatfac is 1, not zero
   qsatfac(:ncol,:top_lev-1) = 1._r8

   ! Zero out values above top_lev for all output variables
   ! Note that elements in proc_rates do not have the extra levels as they are dimensioned to be nlev instead of pver
   tlat(:ncol,:top_lev-1)=0._r8
   qvlat(:ncol,:top_lev-1)=0._r8
   qcten(:ncol,:top_lev-1)=0._r8
   qiten(:ncol,:top_lev-1)=0._r8
   ncten(:ncol,:top_lev-1)=0._r8
   niten(:ncol,:top_lev-1)=0._r8
   qrten(:ncol,:top_lev-1)=0._r8
   qsten(:ncol,:top_lev-1)=0._r8
   nrten(:ncol,:top_lev-1)=0._r8
   nsten(:ncol,:top_lev-1)=0._r8
   qgten(:ncol,:top_lev-1)=0._r8
   ngten(:ncol,:top_lev-1)=0._r8
   rel(:ncol,:top_lev-1)=0._r8
   rel_fn_dum(:ncol,:top_lev-1)=0._r8
   rei(:ncol,:top_lev-1)=0._r8
   sadice(:ncol,:top_lev-1)=0._r8
   sadsnow(:ncol,:top_lev-1)=0._r8
   prect(:ncol)=0._r8
   preci(:ncol)=0._r8
   nevapr(:ncol,:top_lev-1)=0._r8
   am_evp_st(:ncol,:top_lev-1)=0._r8
   prain(:ncol,:top_lev-1)=0._r8
   cmeice(:ncol,:top_lev-1)=0._r8
   dei(:ncol,:top_lev-1)=0._r8
   mu(:ncol,:top_lev-1)=0._r8
   lambdac(:ncol,:top_lev-1)=0._r8
   qsout(:ncol,:top_lev-1)=0._r8
   des(:ncol,:top_lev-1)=0._r8
   qgout(:ncol,:top_lev-1)=0._r8
   ngout(:ncol,:top_lev-1)=0._r8
   dgout(:ncol,:top_lev-1)=0._r8
   cflx(:ncol,:top_lev-1)=0._r8
   iflx(:ncol,:top_lev-1)=0._r8
   gflx(:ncol,:top_lev-1)=0._r8
   rflx(:ncol,:top_lev-1)=0._r8
   sflx(:ncol,:top_lev-1)=0._r8
   qrout(:ncol,:top_lev-1)=0._r8
   reff_rain_dum(:ncol,:top_lev-1)=0._r8
   reff_snow_dum(:ncol,:top_lev-1)=0._r8
   reff_grau_dum(:ncol,:top_lev-1)=0._r8
   nrout(:ncol,:top_lev-1)=0._r8
   nsout(:ncol,:top_lev-1)=0._r8
   refl(:ncol,:top_lev-1)=0._r8
   arefl(:ncol,:top_lev-1)=0._r8
   areflz(:ncol,:top_lev-1)=0._r8
   frefl(:ncol,:top_lev-1)=0._r8
   csrfl(:ncol,:top_lev-1)=0._r8
   acsrfl(:ncol,:top_lev-1)=0._r8
   fcsrfl(:ncol,:top_lev-1)=0._r8
   refl10cm(:ncol,:top_lev-1)=-9999._r8
   reflz10cm(:ncol,:top_lev-1)=0._r8
   rercld(:ncol,:top_lev-1)=0._r8
   ncai(:ncol,:top_lev-1)=0._r8
   ncal(:ncol,:top_lev-1)=0._r8
   qrout2(:ncol,:top_lev-1)=0._r8
   qsout2(:ncol,:top_lev-1)=0._r8
   nrout2(:ncol,:top_lev-1)=0._r8
   nsout2(:ncol,:top_lev-1)=0._r8
   qgout2(:ncol,:top_lev-1)=0._r8
   ngout2(:ncol,:top_lev-1)=0._r8
   dgout2(:ncol,:top_lev-1)=0._r8
   freqg(:ncol,:top_lev-1)=0._r8
   freqs(:ncol,:top_lev-1)=0._r8
   freqr(:ncol,:top_lev-1)=0._r8
   nfice(:ncol,:top_lev-1)=0._r8
   qcrat(:ncol,:top_lev-1)=0._r8
   tnd_qsnow(:ncol,:top_lev-1)=0._r8
   tnd_nsnow(:ncol,:top_lev-1)=0._r8
   re_ice(:ncol,:top_lev-1)=0._r8
   prer_evap(:ncol,:top_lev-1)=0._r8
   frzimm(:ncol,:top_lev-1)=0._r8
   frzcnt(:ncol,:top_lev-1)=0._r8
   frzdep(:ncol,:top_lev-1)=0._r8

   do it = 1, num_steps

     call micro_pumas_tend( &
              ncol,         nlev,           dtime/num_steps,&
              state_loc%t(:ncol,top_lev:),              state_loc%q(:ncol,top_lev:,ixq),            &
              state_loc%q(:ncol,top_lev:,ixcldliq),     state_loc%q(:ncol,top_lev:,ixcldice),          &
              state_loc%q(:ncol,top_lev:,ixnumliq),     state_loc%q(:ncol,top_lev:,ixnumice),       &
              state_loc%q(:ncol,top_lev:,ixrain),       state_loc%q(:ncol,top_lev:,ixsnow),         &
              state_loc%q(:ncol,top_lev:,ixnumrain),    state_loc%q(:ncol,top_lev:,ixnumsnow),      &
              state_loc_graup(:ncol,top_lev:),    state_loc_numgraup(:ncol,top_lev:),     &
              relvar(:ncol,top_lev:),         accre_enhan(:ncol,top_lev:),     &
              state_loc%pmid(:ncol,top_lev:),                state_loc%pdel(:ncol,top_lev:),  state_loc%pint(:ncol,top_lev:), &
              ast(:ncol,top_lev:), alst_mic(:ncol,top_lev:), aist_mic(:ncol,top_lev:), qsatfac(:ncol,top_lev:), &
              rate1cld(:ncol,top_lev:),                         &
              naai(:ncol,top_lev:),            npccn(:ncol,top_lev:),           &
              rndst(:ncol,top_lev:,:),    nacon(:ncol,top_lev:,:),           &
              tlat(:ncol,top_lev:),            qvlat(:ncol,top_lev:),           &
              qcten(:ncol,top_lev:),          qiten(:ncol,top_lev:),          &
              ncten(:ncol,top_lev:),          niten(:ncol,top_lev:),          &
              qrten(:ncol,top_lev:),          qsten(:ncol,top_lev:),          &
              nrten(:ncol,top_lev:),          nsten(:ncol,top_lev:),          &
              qgten(:ncol,top_lev:),          ngten(:ncol,top_lev:),          &
              rel(:ncol,top_lev:),     rel_fn_dum(:ncol,top_lev:),     rei(:ncol,top_lev:),     &
              sadice(:ncol,top_lev:),          sadsnow(:ncol,top_lev:),         &
              prect(:ncol),           preci(:ncol),           &
              nevapr(:ncol,top_lev:),          am_evp_st(:ncol,top_lev:),       &
              prain(:ncol,top_lev:),                   &
              cmeice(:ncol,top_lev:),          dei(:ncol,top_lev:),             &
              mu(:ncol,top_lev:),              lambdac(:ncol,top_lev:),         &
              qsout(:ncol,top_lev:),           des(:ncol,top_lev:),             &
              qgout(:ncol,top_lev:),   ngout(:ncol,top_lev:),   dgout(:ncol,top_lev:),   &
              cflx(:ncol,top_lev:),    iflx(:ncol,top_lev:),                    &
              gflx(:ncol,top_lev:),                                    &
              rflx(:ncol,top_lev:),    sflx(:ncol,top_lev:),    qrout(:ncol,top_lev:),   &
              reff_rain_dum(:ncol,top_lev:),          reff_snow_dum(:ncol,top_lev:),   reff_grau_dum(:ncol,top_lev:),       &
              nrout(:ncol,top_lev:),           nsout(:ncol,top_lev:),           &
              refl(:ncol,top_lev:),    arefl(:ncol,top_lev:),   areflz(:ncol,top_lev:),  &
              frefl(:ncol,top_lev:),   csrfl(:ncol,top_lev:),   acsrfl(:ncol,top_lev:),  &
              fcsrfl(:ncol,top_lev:),   &
              refl10cm(:ncol,top_lev:), reflz10cm(:ncol,top_lev:),    rercld(:ncol,top_lev:),          &
              ncai(:ncol,top_lev:),            ncal(:ncol,top_lev:),            &
              qrout2(:ncol,top_lev:),          qsout2(:ncol,top_lev:),          &
              nrout2(:ncol,top_lev:),          nsout2(:ncol,top_lev:),          &
              drout_dum(:ncol,top_lev:),              dsout2_dum(:ncol,top_lev:),             &
              qgout2(:ncol,top_lev:), ngout2(:ncol,top_lev:), dgout2(:ncol,top_lev:), freqg(:ncol,top_lev:),   &
              freqs(:ncol,top_lev:),           freqr(:ncol,top_lev:),           &
              nfice(:ncol,top_lev:),           qcrat(:ncol,top_lev:),           &
              proc_rates,                                                       &
              errstring, &
              tnd_qsnow(:ncol,top_lev:),tnd_nsnow(:ncol,top_lev:),re_ice(:ncol,top_lev:),&
              prer_evap(:ncol,top_lev:),                                     &
              frzimm(:ncol,top_lev:),  frzcnt(:ncol,top_lev:),  frzdep(:ncol,top_lev:)   )

      call handle_errmsg(errstring, subname="micro_pumas_cam_tend")

      call physics_ptend_init(ptend_loc, psetcols, "micro_pumas", &
                              ls=.true., lq=lq)

      ! Set local tendency.
      ptend_loc%s(:ncol,top_lev:) = tlat(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixq) = qvlat(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixcldliq) = qcten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixcldice) = qiten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixnumliq) = ncten(:ncol,top_lev:)

      if (do_cldice) then
         ptend_loc%q(:ncol,top_lev:,ixnumice) = niten(:ncol,top_lev:)
      else
         ! In this case, the tendency should be all 0.
         if (any(niten(:ncol,:) /= 0._r8)) then
              call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
              " but micro_pumas_tend has ice number tendencies.")
         end if
         ptend_loc%q(:ncol,:,ixnumice) = 0._r8
      end if

      ptend_loc%q(:ncol,top_lev:,ixrain) = qrten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixsnow) = qsten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixnumrain) = nrten(:ncol,top_lev:)
      ptend_loc%q(:ncol,top_lev:,ixnumsnow) = nsten(:ncol,top_lev:)

      if (micro_mg_version > 2) then
         ptend_loc%q(:ncol,top_lev:,ixgraupel) = qgten(:ncol,top_lev:)
         ptend_loc%q(:ncol,top_lev:,ixnumgraupel) = ngten(:ncol,top_lev:)
      end if

      ! Sum into overall ptend
      call physics_ptend_sum(ptend_loc, ptend, ncol)

      ! Update local state
      call physics_update(state_loc, ptend_loc, dtime/num_steps)

      if (trim(micro_mg_warm_rain) == 'tau') then
         proc_rates%amk_c(:ncol,:,:) = proc_rates%amk_c(:ncol,:,:)/num_steps
         proc_rates%ank_c(:ncol,:,:) = proc_rates%ank_c(:ncol,:,:)/num_steps
         proc_rates%amk_r(:ncol,:,:) = proc_rates%amk_r(:ncol,:,:)/num_steps
         proc_rates%ank_r(:ncol,:,:) = proc_rates%ank_r(:ncol,:,:)/num_steps
         proc_rates%amk(:ncol,:,:) = proc_rates%amk(:ncol,:,:)/num_steps
         proc_rates%ank(:ncol,:,:) = proc_rates%ank(:ncol,:,:)/num_steps
         proc_rates%amk_out(:ncol,:,:) = proc_rates%amk_out(:ncol,:,:)/num_steps
      end if

   end do

   ! Divide ptend by substeps.
   call physics_ptend_scale(ptend, 1._r8/num_steps, ncol)

   ! Check to make sure that the microphysics code is respecting the flags that control
   ! whether MG should be prognosing cloud ice and cloud liquid or not.
   if (.not. do_cldice) then
      if (any(ptend%q(:ncol,top_lev:pver,ixcldice) /= 0.0_r8)) &
           call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
           " but micro_pumas_tend has ice mass tendencies.")
      if (any(ptend%q(:ncol,top_lev:pver,ixnumice) /= 0.0_r8)) &
           call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
           " but micro_pumas_tend has ice number tendencies.")
   end if
   if (.not. do_cldliq) then
      if (any(ptend%q(:ncol,top_lev:pver,ixcldliq) /= 0.0_r8)) &
           call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
           " but micro_pumas_tend has liquid mass tendencies.")
      if (any(ptend%q(:ncol,top_lev:pver,ixnumliq) /= 0.0_r8)) &
           call endrun("micro_pumas_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
           " but micro_pumas_tend has liquid number tendencies.")
   end if

   mnuccdohet = 0._r8
   do k=top_lev,pver
      do i=1,ncol
         if (naai(i,k) > 0._r8) then
            mnuccdohet(i,k) = proc_rates%mnuccdtot(i,k-top_lev+1) - (naai_hom(i,k)/naai(i,k))*proc_rates%mnuccdtot(i,k-top_lev+1)
         end if
      end do
   end do

   mgflxprc(:ncol,top_lev:pverp) = rflx(:ncol,top_lev:pverp) + sflx(:ncol,top_lev:pverp)
   mgflxsnw(:ncol,top_lev:pverp) = sflx(:ncol,top_lev:pverp)

   !add condensate fluxes for MG2 (ice and snow already added for MG1)
   if (micro_mg_version >= 2) then
      mgflxprc(:ncol,top_lev:pverp) = mgflxprc(:ncol,top_lev:pverp)+ iflx(:ncol,top_lev:pverp) + cflx(:ncol,top_lev:pverp)
      mgflxsnw(:ncol,top_lev:pverp) = mgflxsnw(:ncol,top_lev:pverp) + iflx(:ncol,top_lev:pverp)
   end if

   !add graupel fluxes for MG3 to snow flux
   if (micro_mg_version >= 3) then
      mgflxprc(:ncol,top_lev:pverp) = mgflxprc(:ncol,top_lev:pverp)+gflx(:ncol,top_lev:pverp)
      mgflxsnw(:ncol,top_lev:pverp) = mgflxsnw(:ncol,top_lev:pverp)+gflx(:ncol,top_lev:pverp)
   end if

   mgmrprc(:ncol,top_lev:pver) = qrout(:ncol,top_lev:pver) + qsout(:ncol,top_lev:pver)
   mgmrsnw(:ncol,top_lev:pver) = qsout(:ncol,top_lev:pver)

   !! calculate effective radius of convective liquid and ice using dcon and deicon (not used by code, not useful for COSP)
   !! hard-coded as average of hard-coded values used for deep/shallow convective detrainment (near line 1502/1505)
   cvreffliq(:ncol,top_lev:pver) = 9.0_r8
   cvreffice(:ncol,top_lev:pver) = 37.0_r8

   ! Reassign rate1 if modal aerosols
   if (rate1_cw2pr_st_idx > 0) then
      rate1ord_cw2pr_st(:ncol,top_lev:pver) = rate1cld(:ncol,top_lev:pver)
   end if

   ! Sedimentation velocity for liquid stratus cloud droplet
   wsedl(:ncol,top_lev:pver) = proc_rates%vtrmc(:ncol,1:nlev)

   ! Microphysical tendencies for use in the macrophysics at the next time step
   CC_T(:ncol,top_lev:pver)    = tlat(:ncol,top_lev:pver)/cpair
   CC_qv(:ncol,top_lev:pver)   = qvlat(:ncol,top_lev:pver)
   CC_ql(:ncol,top_lev:pver)   = qcten(:ncol,top_lev:pver)
   CC_qi(:ncol,top_lev:pver)   = qiten(:ncol,top_lev:pver)
   CC_nl(:ncol,top_lev:pver)   = ncten(:ncol,top_lev:pver)
   CC_ni(:ncol,top_lev:pver)   = niten(:ncol,top_lev:pver)
   CC_qlst(:ncol,top_lev:pver) = qcten(:ncol,top_lev:pver)/max(0.01_r8,alst_mic(:ncol,top_lev:pver))

   ! Net micro_pumas_cam condensation rate
   qme(:ncol,:top_lev-1) = 0._r8
   qme(:ncol,top_lev:pver) = cmeliq(:ncol,top_lev:pver) + proc_rates%cmeitot(:ncol,1:nlev)

   ! For precip, accumulate only total precip in prec_pcw and snow_pcw variables.
   ! Other precip output variables are set to 0
   ! Do not subscript by ncol here, because in physpkg we divide the whole
   ! array and need to avoid an FPE due to uninitialized data.
   prec_pcw = prect
   snow_pcw = preci
   prec_sed = 0._r8
   snow_sed = 0._r8
   prec_str = prec_pcw + prec_sed
   snow_str = snow_pcw + snow_sed

   icecldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)
   liqcldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)

   ! ------------------------------------------------------------ !
   ! Compute in cloud ice and liquid mixing ratios                !
   ! Note that 'iclwp, iciwp' are used for radiation computation. !
   ! ------------------------------------------------------------ !

   icinc = 0._r8
   icwnc = 0._r8
   iciwpst = 0._r8
   iclwpst = 0._r8
   icswp = 0._r8
   cldfsnow = 0._r8
   if (micro_mg_version > 2) then
      icgrauwp = 0._r8
      cldfgrau = 0._r8
   end if

   do k = top_lev, pver
      do i = 1, ncol
         ! Limits for in-cloud mixing ratios consistent with MG microphysics
         ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
         icimrst(i,k)   = min( state_loc%q(i,k,ixcldice) / max(mincld,icecldf(i,k)),0.005_r8 )
         icwmrst(i,k)   = min( state_loc%q(i,k,ixcldliq) / max(mincld,liqcldf(i,k)),0.005_r8 )
         icinc(i,k)     = state_loc%q(i,k,ixnumice) / max(mincld,icecldf(i,k)) * &
              state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
         icwnc(i,k)     = state_loc%q(i,k,ixnumliq) / max(mincld,liqcldf(i,k)) * &
              state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
         ! Calculate micro_pumas_cam cloud water paths in each layer
         ! Note: uses stratiform cloud fraction!
         iciwpst(i,k)   = min(state_loc%q(i,k,ixcldice)/max(mincld,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit
         iclwpst(i,k)   = min(state_loc%q(i,k,ixcldliq)/max(mincld,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit

         ! ------------------------------ !
         ! Adjust cloud fraction for snow !
         ! ------------------------------ !
         cldfsnow(i,k) = cld(i,k)
         ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
         if( ( cldfsnow(i,k) .gt. 1.e-4_r8 ) .and. &
            ( concld(i,k)   .lt. 1.e-4_r8 ) .and. &
            ( state_loc%q(i,k,ixcldliq) .lt. 1.e-10_r8 ) ) then
            cldfsnow(i,k) = 0._r8
         end if
         ! If no cloud and snow, then set to 0.25
         if( ( cldfsnow(i,k) .le. 1.e-4_r8 ) .and. ( qsout(i,k) .gt. 1.e-6_r8 ) ) then
            cldfsnow(i,k) = 0.25_r8
         end if
         ! Calculate in-cloud snow water path
         icswp(i,k) = qsout(i,k) / max( mincld, cldfsnow(i,k) ) * state_loc%pdel(i,k) / gravit

         ! --------------------------------- !
         ! Adjust cloud fraction for graupel !
         ! --------------------------------- !
       if (micro_mg_version > 2) then
          cldfgrau(i,k) = cld(i,k)
         ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
          if( ( cldfgrau(i,k) .gt. 1.e-4_r8 ) .and. &
              ( concld(i,k)   .lt. 1.e-4_r8 ) .and. &
              ( state_loc%q(i,k,ixcldliq) .lt. 1.e-10_r8 ) ) then
              cldfgrau(i,k) = 0._r8
           end if
         ! If no cloud and graupel, then set to 0.25
           if( ( cldfgrau(i,k) .le. 1.e-4_r8 ) .and. ( qgout(i,k) .gt. 1.e-9_r8 ) ) then
              cldfgrau(i,k) = 0.25_r8
           end if

         ! Calculate in-cloud snow water path
           icgrauwp(i,k) = qgout(i,k) / max( 1.e-2_r8, cldfgrau(i,k) ) * state_loc%pdel(i,k) / gravit
        end if

      end do
   end do

   ! Calculate cloud fraction for prognostic precip sizes.
   ! Cloud fraction for purposes of precipitation is maximum cloud
   ! fraction out of all the layers that the precipitation may be
   ! falling down from.
   cldmax(:ncol,top_lev:) = max(mincld, ast(:ncol,top_lev:))
   do k = top_lev+1, pver
      where (state_loc%q(:ncol,k-1,ixrain) >= qsmall .or. &
           state_loc%q(:ncol,k-1,ixsnow) >= qsmall)
         cldmax(:ncol,k) = max(cldmax(:ncol,k-1), cldmax(:ncol,k))
      end where
   end do

   !Copy pbuf field from proc_rates back to pbuf pointer
   bergstot(:ncol,top_lev:) = proc_rates%bergstot(:ncol,1:nlev)
   bergstot(:ncol,1:top_lev-1) = 0._r8

   ! ------------------------------------------------------ !
   ! ------------------------------------------------------ !
   ! All code from here to the end is on grid columns only  !
   ! ------------------------------------------------------ !
   ! ------------------------------------------------------ !

   ! Average the fields which are needed later in this paramterization to be on the grid
   if (use_subcol_microp) then
      call subcol_field_avg(prec_str,  ngrdcol, lchnk, prec_str_grid)
      call subcol_field_avg(iclwpst,   ngrdcol, lchnk, iclwpst_grid)
      call subcol_field_avg(cvreffliq, ngrdcol, lchnk, cvreffliq_grid)
      call subcol_field_avg(cvreffice, ngrdcol, lchnk, cvreffice_grid)
      call subcol_field_avg(mgflxprc,  ngrdcol, lchnk, mgflxprc_grid)
      call subcol_field_avg(mgflxsnw,  ngrdcol, lchnk, mgflxsnw_grid)
      call subcol_field_avg(qme,       ngrdcol, lchnk, qme_grid)
      call subcol_field_avg(nevapr,    ngrdcol, lchnk, nevapr_grid)
      call subcol_field_avg(prain,     ngrdcol, lchnk, prain_grid)
      call subcol_field_avg(proc_rates%evapsnow,  ngrdcol, lchnk, evpsnow_st_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%bergstot,    ngrdcol, lchnk, bergso_grid(:,top_lev:))

      call subcol_field_avg(am_evp_st, ngrdcol, lchnk, am_evp_st_grid)

      ! Average fields which are not in pbuf
      call subcol_field_avg(qrout,     ngrdcol, lchnk, qrout_grid)
      call subcol_field_avg(qsout,     ngrdcol, lchnk, qsout_grid)
      call subcol_field_avg(nsout,     ngrdcol, lchnk, nsout_grid)
      call subcol_field_avg(nrout,     ngrdcol, lchnk, nrout_grid)
      call subcol_field_avg(cld,       ngrdcol, lchnk, cld_grid)
      call subcol_field_avg(proc_rates%qcrestot,    ngrdcol, lchnk, qcreso_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%melttot,     ngrdcol, lchnk, melto_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%mnuccctot,   ngrdcol, lchnk, mnuccco_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%mnuccttot,   ngrdcol, lchnk, mnuccto_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%bergtot,     ngrdcol, lchnk, bergo_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%homotot,     ngrdcol, lchnk, homoo_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%msacwitot,   ngrdcol, lchnk, msacwio_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%psacwstot,   ngrdcol, lchnk, psacwso_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%cmeitot,   ngrdcol, lchnk, cmeiout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qirestot,    ngrdcol, lchnk, qireso_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%prcitot,     ngrdcol, lchnk, prcio_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%praitot,     ngrdcol, lchnk, praio_grid(:,top_lev:))
      call subcol_field_avg(icwmrst,   ngrdcol, lchnk, icwmrst_grid)
      call subcol_field_avg(icimrst,   ngrdcol, lchnk, icimrst_grid)
      call subcol_field_avg(liqcldf,   ngrdcol, lchnk, liqcldf_grid)
      call subcol_field_avg(icecldf,   ngrdcol, lchnk, icecldf_grid)
      call subcol_field_avg(icwnc,     ngrdcol, lchnk, icwnc_grid)
      call subcol_field_avg(icinc,     ngrdcol, lchnk, icinc_grid)
      call subcol_field_avg(state_loc%pdel,            ngrdcol, lchnk, pdel_grid)
      call subcol_field_avg(proc_rates%pratot,      ngrdcol, lchnk, prao_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%prctot,      ngrdcol, lchnk, prco_grid(:,top_lev:))

      call subcol_field_avg(state_loc%q(:,:,ixnumliq), ngrdcol, lchnk, nc_grid(:,top_lev:))
      call subcol_field_avg(state_loc%q(:,:,ixnumice), ngrdcol, lchnk, ni_grid(:,top_lev:))

      call subcol_field_avg(proc_rates%qcsedten,  ngrdcol, lchnk, qcsedtenout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qisedten,  ngrdcol, lchnk, qisedtenout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%vtrmc,     ngrdcol, lchnk, vtrmcout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%vtrmi,     ngrdcol, lchnk, vtrmiout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qcsevap,  ngrdcol, lchnk, qcsevapout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qisevap,  ngrdcol, lchnk, qisevapout_grid(:,top_lev:))

      call subcol_field_avg(cldmax,    ngrdcol, lchnk, cldmax_grid)

      call subcol_field_avg(state_loc%q(:,:,ixrain),    ngrdcol, lchnk, qr_grid)
      call subcol_field_avg(state_loc%q(:,:,ixnumrain), ngrdcol, lchnk, nr_grid)
      call subcol_field_avg(state_loc%q(:,:,ixsnow),    ngrdcol, lchnk, qs_grid)
      call subcol_field_avg(state_loc%q(:,:,ixnumsnow), ngrdcol, lchnk, ns_grid)
      call subcol_field_avg(proc_rates%qrsedten,  ngrdcol, lchnk, qrsedtenout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%qssedten,  ngrdcol, lchnk, qssedtenout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%umr,       ngrdcol, lchnk, umrout_grid(:,top_lev:))
      call subcol_field_avg(proc_rates%ums,       ngrdcol, lchnk, umsout_grid(:,top_lev:))

      if (micro_mg_version > 2) then
            call subcol_field_avg(state_loc%q(:,:,ixgraupel),    ngrdcol, lchnk, qg_grid)
            call subcol_field_avg(state_loc%q(:,:,ixnumgraupel), ngrdcol, lchnk, ng_grid)
            call subcol_field_avg(proc_rates%psacrtot,       ngrdcol, lchnk, psacro_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%pracgtot,       ngrdcol, lchnk, pracgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%psacwgtot,      ngrdcol, lchnk, psacwgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%pgsacwtot,      ngrdcol, lchnk, pgsacwo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%pgracstot,      ngrdcol, lchnk, pgracso_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%prdgtot,        ngrdcol, lchnk, prdgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%qmultgtot,      ngrdcol, lchnk, qmultgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%qmultrgtot,     ngrdcol, lchnk, qmultrgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%npracgtot,      ngrdcol, lchnk, npracgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%nscngtot,       ngrdcol, lchnk, nscngo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%ngracstot,      ngrdcol, lchnk, ngracso_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%nmultgtot,      ngrdcol, lchnk, nmultgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%nmultrgtot,     ngrdcol, lchnk, nmultrgo_grid(:,top_lev:))
            call subcol_field_avg(proc_rates%npsacwgtot,     ngrdcol, lchnk, npsacwgo_grid(:,top_lev:))
      end if

   else
      qcreso_grid(:ncol,:top_lev-1)     = 0._r8
      melto_grid(:ncol,:top_lev-1)      = 0._r8
      mnuccco_grid(:ncol,:top_lev-1)    = 0._r8
      mnuccto_grid(:ncol,:top_lev-1)    = 0._r8
      bergo_grid(:ncol,:top_lev-1)      = 0._r8
      homoo_grid(:ncol,:top_lev-1)      = 0._r8
      msacwio_grid(:ncol,:top_lev-1)    = 0._r8
      psacwso_grid(:ncol,:top_lev-1)    = 0._r8
      cmeiout_grid(:ncol,:top_lev-1)    = 0._r8
      qireso_grid(:ncol,:top_lev-1)     = 0._r8
      prcio_grid(:ncol,:top_lev-1)      = 0._r8
      praio_grid(:ncol,:top_lev-1)      = 0._r8
      prao_grid(:ncol,:top_lev-1)       = 0._r8
      prco_grid(:ncol,:top_lev-1)       = 0._r8
      qcsedtenout_grid(:ncol,:top_lev-1) = 0._r8
      qisedtenout_grid(:ncol,:top_lev-1) = 0._r8
      vtrmcout_grid(:ncol,:top_lev-1)    = 0._r8
      vtrmiout_grid(:ncol,:top_lev-1)    = 0._r8
      qcsevapout_grid(:ncol,:top_lev-1) = 0._r8
      qisevapout_grid(:ncol,:top_lev-1) = 0._r8
      qrsedtenout_grid(:ncol,:top_lev-1) = 0._r8
      qssedtenout_grid(:ncol,:top_lev-1) = 0._r8
      umrout_grid(:ncol,:top_lev-1) = 0._r8
      umsout_grid(:ncol,:top_lev-1) = 0._r8
      psacro_grid(:ncol,:top_lev-1) = 0._r8
      pracgo_grid(:ncol,:top_lev-1) = 0._r8
      psacwgo_grid(:ncol,:top_lev-1) = 0._r8
      pgsacwo_grid(:ncol,:top_lev-1) = 0._r8
      pgracso_grid(:ncol,:top_lev-1) = 0._r8
      prdgo_grid(:ncol,:top_lev-1) = 0._r8
      qmultgo_grid(:ncol,:top_lev-1) = 0._r8
      qmultrgo_grid(:ncol,:top_lev-1) = 0._r8
      npracgo_grid(:ncol,:top_lev-1) = 0._r8
      nscngo_grid(:ncol,:top_lev-1) = 0._r8
      ngracso_grid(:ncol,:top_lev-1) = 0._r8
      nmultgo_grid(:ncol,:top_lev-1) = 0._r8
      nmultrgo_grid(:ncol,:top_lev-1) = 0._r8
      npsacwgo_grid(:ncol,:top_lev-1) = 0._r8
      bergso_grid(:ncol,:top_lev-1) = 0._r8

      ! These pbuf fields need to be assigned.  There is no corresponding subcol_field_avg
      ! as they are reset before being used, so it would be a needless calculation
      lambdac_grid    => lambdac
      mu_grid         => mu
      rel_grid        => rel
      rei_grid        => rei
      sadice_grid     => sadice
      sadsnow_grid    => sadsnow
      dei_grid        => dei
      des_grid        => des
      degrau_grid     => degrau

      ! fields already on grids, so just assign
      prec_str_grid   => prec_str
      iclwpst_grid    => iclwpst
      cvreffliq_grid  => cvreffliq
      cvreffice_grid  => cvreffice
      mgflxprc_grid   => mgflxprc
      mgflxsnw_grid   => mgflxsnw
      qme_grid        => qme
      nevapr_grid     => nevapr
      prain_grid      => prain

      bergso_grid(:ncol,top_lev:)    =  proc_rates%bergstot
      am_evp_st_grid  = am_evp_st

      evpsnow_st_grid(:ncol,top_lev:) = proc_rates%evapsnow
      qrout_grid      = qrout
      qsout_grid      = qsout
      nsout_grid      = nsout
      nrout_grid      = nrout
      cld_grid        = cld
      qcreso_grid(:ncol,top_lev:)     = proc_rates%qcrestot
      melto_grid(:ncol,top_lev:)      = proc_rates%melttot
      mnuccco_grid(:ncol,top_lev:)    = proc_rates%mnuccctot
      mnuccto_grid(:ncol,top_lev:)    = proc_rates%mnuccttot
      bergo_grid(:ncol,top_lev:)      = proc_rates%bergtot
      homoo_grid(:ncol,top_lev:)      = proc_rates%homotot
      msacwio_grid(:ncol,top_lev:)    = proc_rates%msacwitot
      psacwso_grid(:ncol,top_lev:)    = proc_rates%psacwstot
      cmeiout_grid(:ncol,top_lev:)    = proc_rates%cmeitot
      qireso_grid(:ncol,top_lev:)     = proc_rates%qirestot
      prcio_grid(:ncol,top_lev:)      = proc_rates%prcitot
      praio_grid(:ncol,top_lev:)      = proc_rates%praitot
      icwmrst_grid    = icwmrst
      icimrst_grid    = icimrst
      liqcldf_grid    = liqcldf
      icecldf_grid    = icecldf
      icwnc_grid      = icwnc
      icinc_grid      = icinc
      pdel_grid       = state_loc%pdel
      prao_grid(:ncol,top_lev:)       = proc_rates%pratot
      prco_grid(:ncol,top_lev:)       = proc_rates%prctot

      nc_grid = state_loc%q(:,:,ixnumliq)
      ni_grid = state_loc%q(:,:,ixnumice)

      qcsedtenout_grid(:ncol,top_lev:) = proc_rates%qcsedten
      qisedtenout_grid(:ncol,top_lev:) = proc_rates%qisedten
      vtrmcout_grid(:ncol,top_lev:)    = proc_rates%vtrmc
      vtrmiout_grid(:ncol,top_lev:)    = proc_rates%vtrmi
      qcsevapout_grid(:ncol,top_lev:) = proc_rates%qcsevap
      qisevapout_grid(:ncol,top_lev:) = proc_rates%qisevap

      cldmax_grid = cldmax

      qr_grid = state_loc%q(:,:,ixrain)
      nr_grid = state_loc%q(:,:,ixnumrain)
      qs_grid = state_loc%q(:,:,ixsnow)
      ns_grid = state_loc%q(:,:,ixnumsnow)
      qrsedtenout_grid(:ncol,top_lev:) = proc_rates%qrsedten
      qssedtenout_grid(:ncol,top_lev:) = proc_rates%qssedten
      umrout_grid(:ncol,top_lev:) = proc_rates%umr
      umsout_grid(:ncol,top_lev:) = proc_rates%ums

! Zero out terms for budgets if not mg3....
      psacwgo_grid = 0._r8
      pgsacwo_grid = 0._r8
      qmultgo_grid = 0._r8

      if (micro_mg_version > 2) then
            qg_grid = state_loc%q(:,:,ixgraupel)
            ng_grid = state_loc%q(:,:,ixnumgraupel)
            psacro_grid(:ncol,top_lev:) =     proc_rates%psacrtot
            pracgo_grid(:ncol,top_lev:) =     proc_rates%pracgtot
            psacwgo_grid(:ncol,top_lev:) =    proc_rates%psacwgtot
            pgsacwo_grid(:ncol,top_lev:) =    proc_rates%pgsacwtot
            pgracso_grid(:ncol,top_lev:) =    proc_rates%pgracstot
            prdgo_grid(:ncol,top_lev:) =      proc_rates%prdgtot
            qmultgo_grid(:ncol,top_lev:) =    proc_rates%qmultgtot
            qmultrgo_grid(:ncol,top_lev:) =   proc_rates%qmultrgtot
            npracgo_grid(:ncol,top_lev:) =   proc_rates%npracgtot
            nscngo_grid(:ncol,top_lev:) =   proc_rates%nscngtot
            ngracso_grid(:ncol,top_lev:) =   proc_rates%ngracstot
            nmultgo_grid(:ncol,top_lev:) =   proc_rates%nmultgtot
            nmultrgo_grid(:ncol,top_lev:) =   proc_rates%nmultrgtot
            npsacwgo_grid(:ncol,top_lev:) =   proc_rates%npsacwgtot
      end if


   end if

   ! If on subcolumns, average the rest of the pbuf fields which were modified on subcolumns but are not used further in
   ! this parameterization  (no need to assign in the non-subcolumn case -- the else step)
   if (use_subcol_microp) then
      call subcol_field_avg(snow_str,    ngrdcol, lchnk, snow_str_grid)
      call subcol_field_avg(prec_pcw,    ngrdcol, lchnk, prec_pcw_grid)
      call subcol_field_avg(snow_pcw,    ngrdcol, lchnk, snow_pcw_grid)
      call subcol_field_avg(prec_sed,    ngrdcol, lchnk, prec_sed_grid)
      call subcol_field_avg(snow_sed,    ngrdcol, lchnk, snow_sed_grid)
      call subcol_field_avg(cldo,        ngrdcol, lchnk, cldo_grid)
      call subcol_field_avg(mgmrprc,     ngrdcol, lchnk, mgmrprc_grid)
      call subcol_field_avg(mgmrsnw,     ngrdcol, lchnk, mgmrsnw_grid)
      call subcol_field_avg(wsedl,       ngrdcol, lchnk, wsedl_grid)
      call subcol_field_avg(cc_t,        ngrdcol, lchnk, cc_t_grid)
      call subcol_field_avg(cc_qv,       ngrdcol, lchnk, cc_qv_grid)
      call subcol_field_avg(cc_ql,       ngrdcol, lchnk, cc_ql_grid)
      call subcol_field_avg(cc_qi,       ngrdcol, lchnk, cc_qi_grid)
      call subcol_field_avg(cc_nl,       ngrdcol, lchnk, cc_nl_grid)
      call subcol_field_avg(cc_ni,       ngrdcol, lchnk, cc_ni_grid)
      call subcol_field_avg(cc_qlst,     ngrdcol, lchnk, cc_qlst_grid)
      call subcol_field_avg(iciwpst,     ngrdcol, lchnk, iciwpst_grid)
      call subcol_field_avg(icswp,       ngrdcol, lchnk, icswp_grid)
      call subcol_field_avg(cldfsnow,    ngrdcol, lchnk, cldfsnow_grid)

      if (micro_mg_version > 2) then
         call subcol_field_avg(icgrauwp,    ngrdcol, lchnk, icgrauwp_grid)
         call subcol_field_avg(cldfgrau,    ngrdcol, lchnk, cldfsnow_grid)
      end if

      if (rate1_cw2pr_st_idx > 0) then
         call subcol_field_avg(rate1ord_cw2pr_st,    ngrdcol, lchnk, rate1ord_cw2pr_st_grid)
      end if

   end if

   ! ------------------------------------- !
   ! Size distribution calculation         !
   ! ------------------------------------- !

   ! Calculate rho (on subcolumns if turned on) for size distribution
   ! parameter calculations and average it if needed
   !
   ! State instead of state_loc to preserve answers for MG1 (and in any
   ! case, it is unlikely to make much difference).
   rho(:ncol,top_lev:) = state%pmid(:ncol,top_lev:) / &
        (rair*state%t(:ncol,top_lev:))
   if (use_subcol_microp) then
      call subcol_field_avg(rho, ngrdcol, lchnk, rho_grid)
   else
      rho_grid = rho
   end if

   ! Effective radius for cloud liquid, fixed number.
   mu_grid = 0._r8
   lambdac_grid = 0._r8
   rel_fn_grid = 10._r8

   ncic_grid = 1.e8_r8

   do k = top_lev, pver
      !$acc data copyin  (mg_liq_props,icwmrst_grid(:ngrdcol,k),rho_grid(:ngrdcol,k)) &
      !$acc      copy    (ncic_grid(:ngrdcol,k)) &
      !$acc      copyout (mu_grid(:ngrdcol,k),lambdac_grid(:ngrdcol,k))
      call size_dist_param_liq(mg_liq_props, icwmrst_grid(:ngrdcol,k), &
                               ncic_grid(:ngrdcol,k), rho_grid(:ngrdcol,k), &
                               mu_grid(:ngrdcol,k), lambdac_grid(:ngrdcol,k), ngrdcol)
      !$acc end data
   end do

   where (icwmrst_grid(:ngrdcol,top_lev:) > qsmall)
      rel_fn_grid(:ngrdcol,top_lev:) = &
           (mu_grid(:ngrdcol,top_lev:) + 3._r8)/ &
           lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
   end where

   ! Effective radius for cloud liquid, and size parameters
   ! mu_grid and lambdac_grid.
   mu_grid = 0._r8
   lambdac_grid = 0._r8
   rel_grid = 10._r8

   ! Calculate ncic on the grid
   ncic_grid(:ngrdcol,top_lev:) = nc_grid(:ngrdcol,top_lev:) / &
        max(mincld,liqcldf_grid(:ngrdcol,top_lev:))

   do k = top_lev, pver
      !$acc data copyin  (mg_liq_props,icwmrst_grid(:ngrdcol,k), rho_grid(:ngrdcol,k)) &
      !$acc      copy    (ncic_grid(:ngrdcol,k)) &
      !$acc      copyout (mu_grid(:ngrdcol,k),lambdac_grid(:ngrdcol,k))
      call size_dist_param_liq(mg_liq_props, icwmrst_grid(:ngrdcol,k), &
           ncic_grid(:ngrdcol,k), rho_grid(:ngrdcol,k), &
           mu_grid(:ngrdcol,k), lambdac_grid(:ngrdcol,k), ngrdcol)
      !$acc end data
   end do

   where (icwmrst_grid(:ngrdcol,top_lev:) >= qsmall)
      rel_grid(:ngrdcol,top_lev:) = &
           (mu_grid(:ngrdcol,top_lev:) + 3._r8) / &
           lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
   elsewhere
      ! Deal with the fact that size_dist_param_liq sets mu_grid to -100
      ! wherever there is no cloud.
      mu_grid(:ngrdcol,top_lev:) = 0._r8
   end where

   ! Rain/snow effective diameter.
   drout2_grid = 0._r8
   reff_rain_grid = 0._r8
   des_grid = 0._r8
   dsout2_grid = 0._r8
   reff_snow_grid = 0._r8
   reff_grau_grid = 0._r8

   ! Prognostic precipitation

   where (qr_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
      drout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
           qr_grid(:ngrdcol,top_lev:), &
           nr_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
           rho_grid(:ngrdcol,top_lev:), rhow)

      reff_rain_grid(:ngrdcol,top_lev:) = drout2_grid(:ngrdcol,top_lev:) * &
           shapeparam * micron2meter
   end where

   where (qs_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
      dsout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
           qs_grid(:ngrdcol,top_lev:), &
           ns_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
           rho_grid(:ngrdcol,top_lev:), rhosn)

      des_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) *&
           3._r8 * rhosn/rhows

      reff_snow_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) * &
           shapeparam * micron2meter
   end where


! Graupel/Hail size distribution Placeholder
   if (micro_mg_version > 2) then
      degrau_grid = 0._r8
      where (qg_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
         dgout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
              qg_grid(:ngrdcol,top_lev:), &
              ng_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
              rho_grid(:ngrdcol,top_lev:), rhog)

         reff_grau_grid(:ngrdcol,top_lev:) = dgout2_grid(:ngrdcol,top_lev:) * &
              1.5_r8 * 1.e6_r8
         degrau_grid(:ngrdcol,top_lev:) = dgout2_grid(:ngrdcol,top_lev:) *&
              3._r8 * rhog/rhows
      end where
   end if

   ! Effective radius and diameter for cloud ice.
   rei_grid = 25._r8

   niic_grid(:ngrdcol,top_lev:) = ni_grid(:ngrdcol,top_lev:) / &
        max(mincld,icecldf_grid(:ngrdcol,top_lev:))

   do k = top_lev, pver
      !$acc data copyin  (mg_ice_props, icimrst_grid(:ngrdcol,k)) &
      !$acc      copy    (niic_grid(:ngrdcol,k)) &
      !$acc      copyout (rei_grid(:ngrdcol,k))
      call size_dist_param_basic(mg_ice_props,icimrst_grid(:ngrdcol,k), &
                                 niic_grid(:ngrdcol,k),rei_grid(:ngrdcol,k),ngrdcol)
      !$acc end data
   end do

   where (icimrst_grid(:ngrdcol,top_lev:) >= qsmall)
      rei_grid(:ngrdcol,top_lev:) = 1.5_r8/rei_grid(:ngrdcol,top_lev:) &
           * 1.e6_r8
   elsewhere
      rei_grid(:ngrdcol,top_lev:) = 25._r8
   end where

   dei_grid = rei_grid * rhoi/rhows * 2._r8

   ! Limiters for low cloud fraction.
   do k = top_lev, pver
      do i = 1, ngrdcol
         ! Convert snow effective diameter to microns
         des_grid(i,k) = des_grid(i,k) * 1.e6_r8
         if ( ast_grid(i,k) < 1.e-4_r8 ) then
            mu_grid(i,k) = mucon
            lambdac_grid(i,k) = (mucon + 1._r8)/dcon
            dei_grid(i,k) = deicon
         end if
      end do
   end do

   mgreffrain_grid(:ngrdcol,top_lev:pver) = reff_rain_grid(:ngrdcol,top_lev:pver)
   mgreffsnow_grid(:ngrdcol,top_lev:pver) = reff_snow_grid(:ngrdcol,top_lev:pver)

   ! ------------------------------------- !
   ! Precipitation efficiency Calculation  !
   ! ------------------------------------- !

   !-----------------------------------------------------------------------
   ! Liquid water path

   ! Compute liquid water paths, and column condensation
   tgliqwp_grid(:ngrdcol) = 0._r8
   tgcmeliq_grid(:ngrdcol) = 0._r8
   do k = top_lev, pver
      do i = 1, ngrdcol
         tgliqwp_grid(i)  = tgliqwp_grid(i) + iclwpst_grid(i,k)*cld_grid(i,k)

         if (cmeliq_grid(i,k) > 1.e-12_r8) then
            !convert cmeliq to right units:  kgh2o/kgair/s  *  kgair/m2  / kgh2o/m3  = m/s
            tgcmeliq_grid(i) = tgcmeliq_grid(i) + cmeliq_grid(i,k) * &
                 (pdel_grid(i,k) / gravit) / rhoh2o
         end if
      end do
   end do

   ! note: 1e-6 kgho2/kgair/s * 1000. pa / (9.81 m/s2) / 1000 kgh2o/m3 = 1e-7 m/s
   ! this is 1ppmv of h2o in 10hpa
   ! alternatively: 0.1 mm/day * 1.e-4 m/mm * 1/86400 day/s = 1.e-9

   !-----------------------------------------------------------------------
   ! precipitation efficiency calculation  (accumulate cme and precip)

   minlwp = 0.01_r8        !minimum lwp threshold (kg/m3)

   ! zero out precip efficiency and total averaged precip
   pe_grid(:ngrdcol)     = 0._r8
   tpr_grid(:ngrdcol)    = 0._r8
   pefrac_grid(:ngrdcol) = 0._r8

   ! accumulate precip and condensation
   do i = 1, ngrdcol

      acgcme_grid(i)  = acgcme_grid(i) + tgcmeliq_grid(i)
      acprecl_grid(i) = acprecl_grid(i) + prec_str_grid(i)
      acnum_grid(i)   = acnum_grid(i) + 1

      ! if LWP is zero, then 'end of cloud': calculate precip efficiency
      if (tgliqwp_grid(i) < minlwp) then
         if (acprecl_grid(i) > 5.e-8_r8) then
            tpr_grid(i) = max(acprecl_grid(i)/acnum_grid(i), 1.e-15_r8)
            if (acgcme_grid(i) > 1.e-10_r8) then
               pe_grid(i) = min(max(acprecl_grid(i)/acgcme_grid(i), 1.e-15_r8), 1.e5_r8)
               pefrac_grid(i) = 1._r8
            end if
         end if

         ! reset counters
!        if (pe_grid(i) /= 0._r8 .and. (pe_grid(i) < 1.e-8_r8 .or. pe_grid(i) > 1.e3_r8)) then
!           write (iulog,*) 'PE_grid:ANOMALY  pe_grid, acprecl_grid, acgcme_grid, tpr_grid, acnum_grid ', &
!                           pe_grid(i),acprecl_grid(i), acgcme_grid(i), tpr_grid(i), acnum_grid(i)
!        endif

         acprecl_grid(i) = 0._r8
         acgcme_grid(i)  = 0._r8
         acnum_grid(i)   = 0
      end if               ! end LWP zero conditional

      ! if never find any rain....(after 10^3 timesteps...)
      if (acnum_grid(i) > 1000) then
         acnum_grid(i)   = 0
         acprecl_grid(i) = 0._r8
         acgcme_grid(i)  = 0._r8
      end if

   end do

   !-----------------------------------------------------------------------
   ! vertical average of non-zero accretion, autoconversion and ratio.
   ! vars: vprco_grid(i),vprao_grid(i),racau_grid(i),cnt_grid

   vprao_grid = 0._r8
   cnt_grid = 0
   do k = top_lev, pver
      vprao_grid(:ngrdcol) = vprao_grid(:ngrdcol) + prao_grid(:ngrdcol,k)
      where (prao_grid(:ngrdcol,k) /= 0._r8) cnt_grid(:ngrdcol) = cnt_grid(:ngrdcol) + 1
   end do

   where (cnt_grid > 0) vprao_grid = vprao_grid/cnt_grid

   vprco_grid = 0._r8
   cnt_grid = 0
   do k = top_lev, pver
      vprco_grid(:ngrdcol) = vprco_grid(:ngrdcol) + prco_grid(:ngrdcol,k)
      where (prco_grid(:ngrdcol,k) /= 0._r8) cnt_grid(:ngrdcol) = cnt_grid(:ngrdcol) + 1
   end do

   where (cnt_grid > 0)
      vprco_grid = vprco_grid/cnt_grid
      racau_grid = vprao_grid/vprco_grid
   elsewhere
      racau_grid = 0._r8
   end where

   racau_grid = min(racau_grid, 1.e10_r8)

!-----------------------------------------------------------------------
! Diagnostic Rainbow Calculation. Seriously.
!-----------------------------------------------------------------------

   do i = 1, ngrdcol

      top_idx = pver
      convmx = 0._r8
      frlow  = 0._r8
      cldmx  = 0._r8
      cldtot  = maxval(ast(i,top_lev:))

! Find levels in surface layer
      do k = top_lev, pver
         if (state%pmid(i,k) > rb_pmin) then
            top_idx = min(k,top_idx)
         end if
      end do

!For all fractional precip calculated below, use maximum in surface layer.
!For convective precip, base on convective cloud area
      convmx = maxval(concld(i,top_idx:))
!For stratiform precip, base on precip fraction
      cldmx= maxval(freqr(i,top_idx:))
! Combine and use maximum of strat or conv fraction
      frlow= max(cldmx,convmx)

!max precip
      rmax=maxval(qrout_grid(i,top_idx:))

!     Stratiform precip mixing ratio OR some convective precip
!      (rval = true  if any sig precip)

      rval = ((precc(i) > rb_rcmin) .or. (rmax > rb_rmin))

!Now can find conditions for a rainbow:
! Maximum cloud cover (CLDTOT) < 0.5
! 48 < SZA < 90
! freqr (below rb_pmin) > 0.25
! Some rain (liquid > 1.e-6 kg/kg, convective precip > 1.e-7 m/s

      if ((cldtot < 0.5_r8) .and. (sza(i) > 48._r8) .and. (sza(i) < 90._r8) .and. rval) then

!Rainbow 'probability' (area) derived from solid angle theory
!as the fraction of the hemisphere for a spherical cap with angle phi=sza-48.
! This is only valid between 48 < sza < 90 (controlled for above).

          rbfrac(i) =  max(0._r8,(1._r8-COS((sza(i)-48._r8)*deg2rad))/2._r8) * frlow
          rbfreq(i) =  1.0_r8
      end if

   end do                    ! end column loop for rainbows

   call history_out_field('RBFRAC',   rbfrac) ! subcols
   call history_out_field('RBFREQ',   rbfreq) ! subcols


   ! --------------------- !
   ! History Output Fields !
   ! --------------------- !

   ! Column droplet concentration
   cdnumc_grid(:ngrdcol) = sum(nc_grid(:ngrdcol,top_lev:pver) * &
        pdel_grid(:ngrdcol,top_lev:pver)/gravit, dim=2)

   ! Averaging for new output fields
   efcout_grid      = 0._r8
   efiout_grid      = 0._r8
   ncout_grid       = 0._r8
   niout_grid       = 0._r8
   freql_grid       = 0._r8
   freqi_grid       = 0._r8
   icwmrst_grid_out = 0._r8
   icimrst_grid_out = 0._r8
   freqm_grid       = 0._r8
   freqsl_grid      = 0._r8
   freqslm_grid     = 0._r8

   do k = top_lev, pver
      do i = 1, ngrdcol
         if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 5.e-5_r8 ) then
            efcout_grid(i,k) = rel_grid(i,k) * liqcldf_grid(i,k)
            ncout_grid(i,k)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
            freql_grid(i,k)  = liqcldf_grid(i,k)
            icwmrst_grid_out(i,k) = icwmrst_grid(i,k)
         end if
         if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-6_r8 ) then
            efiout_grid(i,k) = rei_grid(i,k) * icecldf_grid(i,k)
            niout_grid(i,k)  = icinc_grid(i,k) * icecldf_grid(i,k)
            freqi_grid(i,k)  = icecldf_grid(i,k)
            icimrst_grid_out(i,k) = icimrst_grid(i,k)
         end if

         ! Supercooled liquid
         if (freql_grid(i,k) > 0.01_r8 .and. freqi_grid(i,k) > 0.01_r8 ) then
            freqm_grid(i,k)=min(liqcldf_grid(i,k),icecldf_grid(i,k))
         end if
         if (freql_grid(i,k) > 0.01_r8 .and. freqi_grid(i,k) < 0.01_r8 .and. state_loc%t(i,k) < tmelt ) then
            freqsl_grid(i,k)=liqcldf_grid(i,k)
         end if
         if (freql_grid(i,k) > 0.01_r8 .and. freqi_grid(i,k) > 0.01_r8 .and. state_loc%t(i,k) < tmelt ) then
            freqslm_grid(i,k)=liqcldf_grid(i,k)
         end if

      end do
   end do

   ! Cloud top effective radius and number.
   fcti_grid  = 0._r8
   fctl_grid  = 0._r8
   ctrel_grid = 0._r8
   ctrei_grid = 0._r8
   ctnl_grid  = 0._r8
   ctni_grid  = 0._r8
   fctm_grid  = 0._r8
   fctsl_grid = 0._r8
   fctslm_grid= 0._r8

   do i = 1, ngrdcol
      do k = top_lev, pver
         if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 1.e-7_r8 ) then
            ctrel_grid(i) = rel_grid(i,k) * liqcldf_grid(i,k)
            ctnl_grid(i)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
            fctl_grid(i)  = liqcldf_grid(i,k)

            ! Cloud Top Mixed phase, supercooled liquid only and supercooled liquid mixed
            if (freqi_grid(i,k) > 0.01_r8) then
               fctm_grid(i)=min(liqcldf_grid(i,k),icecldf_grid(i,k))
            end if
            if (freqi_grid(i,k) < 0.01_r8 .and. state_loc%t(i,k) < tmelt ) then
               fctsl_grid(i)=liqcldf_grid(i,k)
            end if
            if (freqi_grid(i,k) > 0.01_r8 .and. state_loc%t(i,k) < tmelt ) then
               fctslm_grid(i)=liqcldf_grid(i,k)
            end if

            exit
         end if

         if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-7_r8 ) then
            ctrei_grid(i) = rei_grid(i,k) * icecldf_grid(i,k)
            ctni_grid(i)  = icinc_grid(i,k) * icecldf_grid(i,k)
            fcti_grid(i)  = icecldf_grid(i,k)
            exit
         end if
      end do
   end do

   ! Evaporation of stratiform precipitation fields for UNICON
   evprain_st_grid(:ngrdcol,:pver) = nevapr_grid(:ngrdcol,:pver) - evpsnow_st_grid(:ngrdcol,:pver)
   do k = top_lev, pver
      do i = 1, ngrdcol
         evprain_st_grid(i,k) = max(evprain_st_grid(i,k), 0._r8)
         evpsnow_st_grid(i,k) = max(evpsnow_st_grid(i,k), 0._r8)
      end do
   end do

   ! Assign the values to the pbuf pointers if they exist in pbuf
   if (qrain_idx > 0)  qrout_grid_ptr = qrout_grid
   if (qsnow_idx > 0)  qsout_grid_ptr = qsout_grid
   if (nrain_idx > 0)  nrout_grid_ptr = nrout_grid
   if (nsnow_idx > 0)  nsout_grid_ptr = nsout_grid
   if (qcsedten_idx > 0) qcsedtenout_grid_ptr = qcsedtenout_grid
   if (qrsedten_idx > 0) qrsedtenout_grid_ptr = qrsedtenout_grid
   if (qisedten_idx > 0) qisedtenout_grid_ptr = qisedtenout_grid
   if (qssedten_idx > 0) qssedtenout_grid_ptr = qssedtenout_grid
   if (vtrmc_idx > 0)    vtrmcout_grid_ptr    = vtrmcout_grid
   if (umr_idx > 0)      umrout_grid_ptr      = umrout_grid
   if (vtrmi_idx > 0)    vtrmiout_grid_ptr    = vtrmiout_grid
   if (ums_idx > 0)      umsout_grid_ptr      = umsout_grid
   if (qcsevap_idx > 0 ) qcsevapout_grid_ptr  = qcsevapout_grid
   if (qisevap_idx > 0 ) qisevapout_grid_ptr  = qisevapout_grid

   ! --------------------------------------------- !
   ! General outfield calls for microphysics       !
   ! --------------------------------------------- !

   ! Output a handle of variables which are calculated on the fly

   ftem_grid = 0._r8

   ftem_grid(:ngrdcol,top_lev:pver) =  qcreso_grid(:ngrdcol,top_lev:pver)
   call history_out_field( 'MPDW2V', ftem_grid)

   ftem_grid(:ngrdcol,top_lev:pver) =  melto_grid(:ngrdcol,top_lev:pver) - mnuccco_grid(:ngrdcol,top_lev:pver)&
        - mnuccto_grid(:ngrdcol,top_lev:pver) -  bergo_grid(:ngrdcol,top_lev:pver) - homoo_grid(:ngrdcol,top_lev:pver)&
        - msacwio_grid(:ngrdcol,top_lev:pver)
   call history_out_field( 'MPDW2I', ftem_grid)

   if (micro_mg_version > 2) then
      ftem_grid(:ngrdcol,top_lev:pver) = -prao_grid(:ngrdcol,top_lev:pver) - prco_grid(:ngrdcol,top_lev:pver)&
          - psacwso_grid(:ngrdcol,top_lev:pver) - bergso_grid(:ngrdcol,top_lev:pver)&
          - psacwgo_grid(:ngrdcol,top_lev:pver) - pgsacwo_grid(:ngrdcol,top_lev:pver)
   else
      ftem_grid(:ngrdcol,top_lev:pver) = -prao_grid(:ngrdcol,top_lev:pver) - prco_grid(:ngrdcol,top_lev:pver)&
          - psacwso_grid(:ngrdcol,top_lev:pver) - bergso_grid(:ngrdcol,top_lev:pver)
   endif

   call history_out_field( 'MPDW2P', ftem_grid)

   ftem_grid(:ngrdcol,top_lev:pver) =  cmeiout_grid(:ngrdcol,top_lev:pver) + qireso_grid(:ngrdcol,top_lev:pver)
   call history_out_field( 'MPDI2V', ftem_grid)

   if (micro_mg_version > 2) then
      ftem_grid(:ngrdcol,top_lev:pver) = -melto_grid(:ngrdcol,top_lev:pver) + mnuccco_grid(:ngrdcol,top_lev:pver) &
          + mnuccto_grid(:ngrdcol,top_lev:pver) +  bergo_grid(:ngrdcol,top_lev:pver) + homoo_grid(:ngrdcol,top_lev:pver)&
          + msacwio_grid(:ngrdcol,top_lev:pver)&
          - qmultgo_grid(:ngrdcol,top_lev:pver)
   else
      ftem_grid(:ngrdcol,top_lev:pver) = -melto_grid(:ngrdcol,top_lev:pver) + mnuccco_grid(:ngrdcol,top_lev:pver) &
          + mnuccto_grid(:ngrdcol,top_lev:pver) +  bergo_grid(:ngrdcol,top_lev:pver) + homoo_grid(:ngrdcol,top_lev:pver)&
          + msacwio_grid(:ngrdcol,top_lev:pver)
   endif

   call history_out_field( 'MPDI2W', ftem_grid)

   ftem_grid(:ngrdcol,top_lev:pver) = -prcio_grid(:ngrdcol,top_lev:pver) - praio_grid(:ngrdcol,top_lev:pver)
   call history_out_field( 'MPDI2P', ftem_grid)

   ! Output fields which have not been averaged already, averaging if use_subcol_microp is true
   if (trim(micro_mg_warm_rain) == 'tau' .or. trim(micro_mg_warm_rain) == 'emulated') then
      call history_out_field('scale_qc',    proc_rates%scale_qc) ! subcols
      call history_out_field('scale_nc',    proc_rates%scale_nc) ! subcols
      call history_out_field('scale_qr',    proc_rates%scale_qr) ! subcols
      call history_out_field('scale_nr',    proc_rates%scale_nr) ! subcols
      call history_out_field('amk_c',       proc_rates%amk_c) ! subcols
      call history_out_field('ank_c',       proc_rates%ank_c) ! subcols
      call history_out_field('amk_r',       proc_rates%amk_r) ! subcols
      call history_out_field('ank_r',       proc_rates%ank_r) ! subcols
      call history_out_field('amk',         proc_rates%amk) ! subcols
      call history_out_field('ank',         proc_rates%ank) ! subcols
      call history_out_field('amk_out',     proc_rates%amk_out) ! subcols
      call history_out_field('ank_out',     proc_rates%ank_out) ! subcols
      call history_out_field('QC_TAU_out',  proc_rates%qc_out_TAU) ! subcols
      call history_out_field('NC_TAU_out',  proc_rates%nc_out_TAU) ! subcols
      call history_out_field('QR_TAU_out',  proc_rates%qr_out_TAU) ! subcols
      call history_out_field('NR_TAU_out',  proc_rates%nr_out_TAU) ! subcols
      call history_out_field('qctend_TAU',  proc_rates%qctend_TAU) ! subcols
      call history_out_field('nctend_TAU',  proc_rates%nctend_TAU) ! subcols
      call history_out_field('qrtend_TAU',  proc_rates%qrtend_TAU) ! subcols
      call history_out_field('nrtend_TAU',  proc_rates%nrtend_TAU) ! subcols
      call history_out_field('gmnnn_lmnnn_TAU',  proc_rates%gmnnn_lmnnn_TAU) ! subcols
      call history_out_field('ML_fixer',     proc_rates%ML_fixer) ! subcols
      call history_out_field('qc_fixer',     proc_rates%qc_fixer) ! subcols
      call history_out_field('nc_fixer',     proc_rates%nc_fixer) ! subcols
      call history_out_field('qr_fixer',     proc_rates%qr_fixer) ! subcols
      call history_out_field('nr_fixer',     proc_rates%nr_fixer) ! subcols
      call history_out_field('QC_TAU_in',   proc_rates%qc_in_TAU) ! subcols
      call history_out_field('NC_TAU_in',   proc_rates%nc_in_TAU) ! subcols
      call history_out_field('QR_TAU_in',   proc_rates%qr_in_TAU) ! subcols
      call history_out_field('NR_TAU_in',   proc_rates%nr_in_TAU) ! subcols
   end if

   if (trim(micro_mg_warm_rain) == 'sb2001') then
      call history_out_field('qctend_SB2001',  proc_rates%qctend_SB2001) ! subcols
      call history_out_field('nctend_SB2001',  proc_rates%nctend_SB2001) ! subcols
      call history_out_field('qrtend_SB2001',  proc_rates%qrtend_SB2001) ! subcols
      call history_out_field('nrtend_SB2001',  proc_rates%nrtend_SB2001) ! subcols
   end if
   if (trim(micro_mg_warm_rain) == 'kk2000') then
      call history_out_field('qctend_KK2000',  proc_rates%qctend_KK2000) ! subcols
      call history_out_field('nctend_KK2000',  proc_rates%nctend_KK2000) ! subcols
      call history_out_field('qrtend_KK2000',  proc_rates%qrtend_KK2000) ! subcols
      call history_out_field('nrtend_KK2000',  proc_rates%nrtend_KK2000) ! subcols
   end if
   call history_out_field('LAMC',  proc_rates%lamc_out) ! subcols
   call history_out_field('LAMR',  proc_rates%lamr_out) ! subcols
   call history_out_field('PGAM',  proc_rates%pgam_out) ! subcols
   call history_out_field('N0R',  proc_rates%n0r_out) ! subcols

   call history_out_field('MPICLWPI',    iclwpi) ! subcols
   call history_out_field('MPICIWPI',    iciwpi) ! subcols
   call history_out_field('REFL',        refl) ! subcols
   call history_out_field('AREFL',       arefl) ! subcols
   call history_out_field('AREFLZ',      areflz) ! subcols
   call history_out_field('FREFL',       frefl) ! subcols
   call history_out_field('CSRFL',       csrfl) ! subcols
   call history_out_field('ACSRFL',      acsrfl) ! subcols
   call history_out_field('FCSRFL',      fcsrfl) ! subcols
   call history_out_field('REFL10CM',    refl10cm) ! subcols
   call history_out_field('REFLZ10CM',   reflz10cm) ! subcols
   call history_out_field('RERCLD',      rercld) ! subcols
   call history_out_field('NCAL',        ncal) ! subcols
   call history_out_field('NCAI',        ncai) ! subcols
   call history_out_field('AQRAIN',      qrout2) ! subcols
   call history_out_field('AQSNOW',      qsout2) ! subcols
   call history_out_field('ANRAIN',      nrout2) ! subcols
   call history_out_field('ANSNOW',      nsout2) ! subcols
   call history_out_field('FREQR',       freqr) ! subcols
   call history_out_field('FREQS',       freqs) ! subcols
   call history_out_field('MPDT',        tlat) ! subcols
   call history_out_field('MPDQ',        qvlat) ! subcols
   call history_out_field('MPDLIQ',      qcten) ! subcols
   call history_out_field('MPDICE',      qiten) ! subcols
   call history_out_field('MPDNLIQ',     ncten) ! subcols
   call history_out_field('MPDNICE',     niten) ! subcols
   call history_out_field('EVAPSNOW',    proc_rates%evapsnow) ! subcols
   call history_out_field('QCSEVAP',     proc_rates%qcsevap) ! subcols
   call history_out_field('QISEVAP',     proc_rates%qisevap) ! subcols
   call history_out_field('QVRES',       proc_rates%qvres) ! subcols
   call history_out_field('VTRMC',       proc_rates%vtrmc) ! subcols
   call history_out_field('VTRMI',       proc_rates%vtrmi) ! subcols
   call history_out_field('QCSEDTEN',    proc_rates%qcsedten) ! subcols
   call history_out_field('QISEDTEN',    proc_rates%qisedten) ! subcols
   call history_out_field('QRSEDTEN',    proc_rates%qrsedten) ! subcols
   call history_out_field('QSSEDTEN',    proc_rates%qssedten) ! subcols
   call history_out_field('MNUCCRIO',    proc_rates%mnuccritot) ! subcols
   call history_out_field('MNUDEPO',     proc_rates%mnudeptot) ! subcols
   call history_out_field('MELTSTOT',    proc_rates%meltstot) ! subcols
   call history_out_field('MNUCCDO',     proc_rates%mnuccdtot) ! subcols
   call history_out_field('MNUCCDOhet',  mnuccdohet) ! subcols
   call history_out_field('MNUCCRO',     proc_rates%mnuccrtot) ! subcols
   call history_out_field('PRACSO',      proc_rates%pracstot ) ! subcols
   call history_out_field('VAPDEPSO',    proc_rates%vapdepstot) ! subcols
   call history_out_field('MELTSDT',     proc_rates%meltsdttot) ! subcols
   call history_out_field('FRZRDT',      proc_rates%frzrdttot ) ! subcols
   call history_out_field('FICE',        nfice) ! subcols
   call history_out_field('CLDFSNOW',    cldfsnow) ! subcols
   call history_out_field ('NNUCCCO',  proc_rates%nnuccctot  ) ! subcols
   call history_out_field ('NNUCCTO',  proc_rates%nnuccttot  ) ! subcols
   call history_out_field ('NNUCCDO',  proc_rates%nnuccdtot  ) ! subcols
   call history_out_field ('NNUDEPO',  proc_rates%nnudeptot  ) ! subcols
   call history_out_field ('NHOMO',    proc_rates%nhomotot   ) ! subcols
   call history_out_field ('NNUCCRO',  proc_rates%nnuccrtot  ) ! subcols
   call history_out_field ('NNUCCRIO', proc_rates%nnuccritot ) ! subcols
   call history_out_field ('NSACWIO',  proc_rates%nsacwitot  ) ! subcols
   call history_out_field ('NPRAO',    proc_rates%npratot    ) ! subcols
   call history_out_field ('NPSACWSO', proc_rates%npsacwstot ) ! subcols
   call history_out_field ('NPRAIO',   proc_rates%npraitot   ) ! subcols
   call history_out_field ('NPRACSO',  proc_rates%npracstot  ) ! subcols
   call history_out_field ('NPRCO',    proc_rates%nprctot    ) ! subcols
   call history_out_field ('NPRCIO',   proc_rates%nprcitot   ) ! subcols
   call history_out_field ('NCSEDTEN', proc_rates%ncsedten ) ! subcols
   call history_out_field ('NISEDTEN', proc_rates%nisedten ) ! subcols
   call history_out_field ('NRSEDTEN', proc_rates%nrsedten ) ! subcols
   call history_out_field ('NSSEDTEN', proc_rates%nssedten ) ! subcols
   call history_out_field ('NMELTO',   proc_rates%nmelttot   ) ! subcols
   call history_out_field ('NMELTS',   proc_rates%nmeltstot  ) ! subcols

   call history_out_field('UMR',      proc_rates%umr) ! subcols
   call history_out_field('UMS',      proc_rates%ums) ! subcols

   call history_out_field('QCRAT',    qcrat) ! subcols

   if (micro_mg_version > 2) then
      call history_out_field('UMG',        proc_rates%umg) ! subcols
      call history_out_field('QGSEDTEN',   proc_rates%qgsedten) ! subcols
      call history_out_field('FREQG',       freqg) ! subcols
      call history_out_field('AQGRAU',      qgout2) ! subcols
      call history_out_field('ANGRAU',      ngout2) ! subcols
      call history_out_field('CLDFGRAU',    cldfgrau) ! subcols
      call history_out_field('MELTGTOT',    proc_rates%meltgtot) ! subcols
      call history_out_field('NMELTG',      proc_rates%nmeltgtot) ! subcols
      call history_out_field('NGSEDTEN',    proc_rates%ngsedten ) ! subcols

   end if

   ! Example subcolumn history_out_field call
   if (use_subcol_microp) then
      call history_out_field('FICE_SCOL',   nfice) ! subcols
      call history_out_field('MPDLIQ_SCOL', ptend%q(:,:,ixcldliq)) ! subcols
      call history_out_field('MPDICE_SCOL', qiten) ! subcols
   end if

   ! Output fields which are already on the grid
   call history_out_field('QRAIN',       qrout_grid)
   call history_out_field('QSNOW',       qsout_grid)
   call history_out_field('NRAIN',       nrout_grid)
   call history_out_field('NSNOW',       nsout_grid)
   call history_out_field('CV_REFFLIQ',  cvreffliq_grid)
   call history_out_field('CV_REFFICE',  cvreffice_grid)
   call history_out_field('LS_FLXPRC',   mgflxprc_grid)
   call history_out_field('LS_FLXSNW',   mgflxsnw_grid)
   call history_out_field('CME',         qme_grid)
   call history_out_field('PRODPREC',    prain_grid)
   call history_out_field('EVAPPREC',    nevapr_grid)
   call history_out_field('QCRESO',      qcreso_grid)
   call history_out_field('LS_REFFRAIN', mgreffrain_grid)
   call history_out_field('LS_REFFSNOW', mgreffsnow_grid)
   call history_out_field('DSNOW',       des_grid)
   call history_out_field('ADRAIN',      drout2_grid)
   call history_out_field('ADSNOW',      dsout2_grid)
   call history_out_field('PE',          pe_grid)
   call history_out_field('PEFRAC',      pefrac_grid)
   call history_out_field('APRL',        tpr_grid)
   call history_out_field('VPRAO',       vprao_grid)
   call history_out_field('VPRCO',       vprco_grid)
   call history_out_field('RACAU',       racau_grid)
   call history_out_field('AREL',        efcout_grid)
   call history_out_field('AREI',        efiout_grid)
   call history_out_field('AWNC' ,       ncout_grid)
   call history_out_field('AWNI' ,       niout_grid)
   call history_out_field('FREQL',       freql_grid)
   call history_out_field('FREQI',       freqi_grid)
   call history_out_field('ACTREL',      ctrel_grid)
   call history_out_field('ACTREI',      ctrei_grid)
   call history_out_field('ACTNL',       ctnl_grid)
   call history_out_field('ACTNI',       ctni_grid)
   call history_out_field('FCTL',        fctl_grid)
   call history_out_field('FCTI',        fcti_grid)
   call history_out_field('ICINC',       icinc_grid)
   call history_out_field('ICWNC',       icwnc_grid)
   call history_out_field('EFFLIQ_IND',  rel_fn_grid)
   call history_out_field('CDNUMC',      cdnumc_grid)
   call history_out_field('REL',         rel_grid)
   call history_out_field('REI',         rei_grid)
   call history_out_field('MG_SADICE',   sadice_grid)
   call history_out_field('MG_SADSNOW',  sadsnow_grid)
   call history_out_field('ICIMRST',     icimrst_grid_out)
   call history_out_field('ICWMRST',     icwmrst_grid_out)
   call history_out_field('CMEIOUT',     cmeiout_grid)
   call history_out_field('PRAO',        prao_grid)
   call history_out_field('PRCO',        prco_grid)
   call history_out_field('MNUCCCO',     mnuccco_grid)
   call history_out_field('MNUCCTO',     mnuccto_grid)
   call history_out_field('MSACWIO',     msacwio_grid)
   call history_out_field('PSACWSO',     psacwso_grid)
   call history_out_field('BERGSO',      bergso_grid)
   call history_out_field('BERGO',       bergo_grid)
   call history_out_field('MELTO',       melto_grid)
   call history_out_field('HOMOO',       homoo_grid)
   call history_out_field('PRCIO',       prcio_grid)
   call history_out_field('PRAIO',       praio_grid)
   call history_out_field('QIRESO',      qireso_grid)
   call history_out_field('FREQM',       freqm_grid)
   call history_out_field('FREQSL',      freqsl_grid)
   call history_out_field('FREQSLM',     freqslm_grid)
   call history_out_field('FCTM',        fctm_grid)
   call history_out_field('FCTSL',       fctsl_grid)
   call history_out_field('FCTSLM',      fctslm_grid)

   if (micro_mg_version > 2) then
      call history_out_field('PRACGO',      pracgo_grid)
      call history_out_field('PSACRO',      psacro_grid)
      call history_out_field('PSACWGO',     psacwgo_grid)
      call history_out_field('PGSACWO',     pgsacwo_grid)
      call history_out_field('PGRACSO',     pgracso_grid)
      call history_out_field('PRDGO',       prdgo_grid)
      call history_out_field('QMULTGO',     qmultgo_grid)
      call history_out_field('QMULTRGO',    qmultrgo_grid)
      call history_out_field('LS_REFFGRAU', reff_grau_grid)
      call history_out_field ('NPRACGO',    npracgo_grid)
      call history_out_field ('NSCNGO',     nscngo_grid)
      call history_out_field ('NGRACSO',    ngracso_grid)
      call history_out_field ('NMULTGO',    nmultgo_grid)
      call history_out_field ('NMULTRGO',   nmultrgo_grid)
      call history_out_field ('NPSACWGO',   npsacwgo_grid)
   end if

   if (micro_mg_adjust_cpt) then
      cp_rh(:ncol, :pver)  = 0._r8

      do i = 1, ncol

         ! Calculate the RH including any T change that we make.
         do k = top_lev, pver
           call qsat(state_loc%t(i,k), state_loc%pmid(i,k), es, qs)
           cp_rh(i,k) = state_loc%q(i, k, ixq) / qs * 100._r8
         end do
      end do

      call history_out_field("TROPF_RHADJ", cp_rh)
   end if

   ! deallocate the temporary pbuf grid variable which was allocated if subcolumns are not used
   if (.not. use_subcol_microp) then
      deallocate(bergso_grid)
   end if

   ! deallocate the proc_rates DDT
   call proc_rates%deallocate(micro_mg_warm_rain)

   ! ptend_loc is deallocated in physics_update above
   call physics_state_dealloc(state_loc)

   if (qsatfac_idx <= 0) then
      deallocate(qsatfac)
   end if

end subroutine pumas_diagnostics_run

end module pumas_diagnostics
