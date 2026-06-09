module clubb

  use ccpp_kinds, only: kind_phys

  implicit none

  private

  save

  ! Subroutines to make public
  public :: clubb_init, clubb2_run, clubb3_run, stats_zero

  contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_init(pcols, pver, pverp, begchunk, endchunk, masterproc, &
                        mpicom, mstrid, mpi_character, max_fieldname_len, &
                        iulog, subcol_scheme, pcnst, cnst_ndropmixed, lq, &
                        stats_metadata, l_input_fields, clubb_config_flags, &
                        clubb_l_do_expldiff_rtm_thlm, l_implemented, sclr_idx, &
                        clubb_C2rtthl, clubb_C8, clubb_c11, clubb_c11b, clubb_c14, &
                        clubb_C_wp3_pr_turb, clubb_c_K10, clubb_mult_coef, &
                        clubb_Skw_denom_coef, clubb_C2rt, clubb_C2thl, clubb_beta, &
                        clubb_c6rt, clubb_c6rtb, clubb_c6rtc, clubb_c6thl, clubb_c6thlb, &
                        clubb_c6thlc, clubb_wpxp_L_thresh, clubb_C7, clubb_C7b, &
                        clubb_gamma_coef, clubb_c_K10h, clubb_lambda0_stability_coef, &
                        clubb_lmin_coef, clubb_C8b, clubb_skw_max_mag, clubb_C1, clubb_C1b, &
                        clubb_gamma_coefb, clubb_up2_sfc_coef, clubb_C4, clubb_C_uu_shr, &
                        clubb_C_uu_buoy, clubb_c_K1, clubb_c_K2, clubb_nu2, clubb_c_K8, &
                        clubb_c_K9, clubb_nu9, clubb_C_wp2_splat, clubb_C_invrs_tau_bkgnd, &
                        clubb_C_invrs_tau_sfc, clubb_C_invrs_tau_shear, clubb_C_invrs_tau_N2, &
                        clubb_C_invrs_tau_N2_wp2, clubb_C_invrs_tau_N2_xp2, clubb_C_invrs_tau_N2_wpxp, &
                        clubb_C_invrs_tau_N2_clear_wp3, clubb_bv_efold, clubb_wpxp_Ri_exp, clubb_z_displace, &
                        edsclr_dim, sclr_dim, hydromet_dim, &
                        nzm_clubb, nzt_clubb, hm_metadata, clubb_params_single_col, &
                        clubb_vars_zt, clubb_vars_zm, clubb_vars_sfc, clubb_vars_rad_zt, clubb_vars_rad_zm, & 
                        stats_zt, stats_zm, stats_sfc, stats_rad_zt, stats_rad_zm, &
                        pdf_params_chnk, pdf_params_zm_chnk, pdf_implicit_coefs_terms_chnk, &
                        out_zt, out_zm, out_sfc, out_radzt, out_radzm, &
                        error_message, clubb_init_errcode )

!-------------------------------------------------------------------------------
! Description:
!   Initialize UWM CLUBB.
! Author: Cheryl Craig March 2011
! Modifications: Pete Bogenschutz 2011 March and onward
! Modifications: K Thayer-Calder 2013 July and onward
! Origin: Based heavily on UWM clubb_init.F90
! References:
!   None
!-------------------------------------------------------------------------------

    ! These are needed to set parameters
    use clubb_api_module, only: &
         core_rknd, em_min, &
         ilambda0_stability_coef, ic_K10, ic_K10h, iC7, iC7b, iC8, iC8b, iC11, iC11b, iC4, iC_uu_shr, iC_uu_buoy, &
         iC1, iC1b, iC6rt, iC6rtb, iC6rtc, iC6thl, iC6thlb, iC6thlc, iup2_sfc_coef, iwpxp_L_thresh, &
         iC14, iC_wp3_pr_turb, igamma_coef, igamma_coefb, imult_coef, ilmin_coef, &
         iSkw_denom_coef, ibeta, iskw_max_mag, &
         iC_invrs_tau_bkgnd,iC_invrs_tau_sfc,iC_invrs_tau_shear,iC_invrs_tau_N2,iC_invrs_tau_N2_wp2, &
         iC_invrs_tau_N2_xp2,iC_invrs_tau_N2_wpxp,iC_invrs_tau_N2_clear_wp3, &
         iC2rt, iC2thl, iC2rtthl, ic_K1, ic_K2, inu2, ic_K8, ic_K9, inu9, iC_wp2_splat, ibv_efold, &
         iwpxp_Ri_exp, iz_displace, &
         params_list

    use clubb_api_module, only: &
         nvarmax_zm, &
         nvarmax_zt, &
         nvarmax_rad_zt, &
         nvarmax_rad_zm, &
         nvarmax_sfc, &
         var_length, &
         print_clubb_config_flags_api, &
         check_clubb_settings_api, &
         init_pdf_params_api, &
         time_precision, &
         stats_metadata_type, &
         clubb_config_flags_type, &
         sclr_idx_type, &
         hm_metadata_type, &
         stats, &
         pdf_parameter, &
         implicit_coefs_terms, &
         set_clubb_debug_level_api, &
         clubb_fatal_error, &     ! Error code value to indicate a fatal error
         err_info_type, &
         init_default_err_info_api, &
         cleanup_err_info_api, &
         nparams, &
         init_clubb_params_api, &
         saturation_gfdl,   & ! Constant for the GFDL approximation of saturation
         saturation_flatau    ! Constant for Flatau approximations of saturation


    implicit none

    !  Input Variables

    integer, intent(in) :: iulog, pcnst, max_fieldname_len
    integer, intent(in) :: mpicom, mpi_character, mstrid
    integer, intent(inout) :: edsclr_dim, nzt_clubb, nzm_clubb
    integer, intent(in) :: sclr_dim, hydromet_dim
    logical, intent(in) :: masterproc, l_implemented
    logical, intent(in) :: cnst_ndropmixed(:)
    logical, intent(inout) :: lq(:)
    character(len=16), intent(in) :: subcol_scheme

    real(kind=time_precision) :: dum1, dum2, dum3

    type(err_info_type) :: &
      err_info          ! err_info struct used in CLUBB containing err_code and err_header

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type(clubb_config_flags_type), intent(inout) :: &
      clubb_config_flags

    type (sclr_idx_type), intent(inout) :: &
      sclr_idx

    type (hm_metadata_type), intent(inout) :: &
      hm_metadata

    real(kind=kind_phys), intent(inout) :: &
      clubb_params_single_col(:,:)    ! Adjustable CLUBB parameters (C1, C2 ...)

    integer, intent(in) :: pcols, pver, pverp, begchunk, endchunk

    ! Variables that contains all the statistics
    type (stats), intent(inout) :: &
      stats_zt(:),      & ! stats_zt grid
      stats_zm(:),      & ! stats_zm grid
      stats_rad_zt(:),  & ! stats_rad_zt grid
      stats_rad_zm(:),  & ! stats_rad_zm grid
      stats_sfc(:)        ! stats_sfc

    type(pdf_parameter), allocatable, intent(inout) :: &
      pdf_params_chnk(:)                ! PDF parameters (thermo. levs.) [units vary]
  
    type(pdf_parameter), allocatable, intent(inout) :: &
      pdf_params_zm_chnk(:)             ! PDF parameters on momentum levs. [units vary]
 
    type(implicit_coefs_terms), allocatable, intent(inout) :: &
      pdf_implicit_coefs_terms_chnk(:)  ! PDF impl. coefs. & expl. terms      [units vary]

    character(len=var_length), intent(in) :: clubb_vars_zt(:)      ! Variables on the thermodynamic levels
    character(len=var_length), intent(in) :: clubb_vars_zm(:)      ! Variables on the momentum levels
    character(len=var_length), intent(in) :: clubb_vars_rad_zt(:)  ! Variables on the radiation levels
    character(len=var_length), intent(in) :: clubb_vars_rad_zm(:)  ! Variables on the radiation levels
    character(len=var_length), intent(in) :: clubb_vars_sfc(:)     ! Variables at the model surface
 
    real(kind=kind_phys), allocatable, intent(inout) :: &
      out_zt(:,:,:), out_zm(:,:,:), out_radzt(:,:,:), out_radzm(:,:,:), out_sfc(:,:,:)

    logical, intent(in) :: clubb_l_do_expldiff_rtm_thlm

    integer :: j, m

    logical, intent(in) :: l_input_fields

    integer :: ierr = 0
 
    integer, intent(out) :: clubb_init_errcode
    character(len=200), intent(inout) :: error_message 

    real(kind=kind_phys), intent(inout) :: clubb_C2rtthl
    real(kind=kind_phys), intent(inout) :: clubb_C8
    real(kind=kind_phys), intent(inout) :: clubb_c11
    real(kind=kind_phys), intent(inout) :: clubb_c11b
    real(kind=kind_phys), intent(inout) :: clubb_c14
    real(kind=kind_phys), intent(inout) :: clubb_C_wp3_pr_turb
    real(kind=kind_phys), intent(inout) :: clubb_c_K10
    real(kind=kind_phys), intent(inout) :: clubb_mult_coef
    real(kind=kind_phys), intent(inout) :: clubb_Skw_denom_coef
    real(kind=kind_phys), intent(inout) :: clubb_C2rt
    real(kind=kind_phys), intent(inout) :: clubb_C2thl
    real(kind=kind_phys), intent(inout) :: clubb_beta
    real(kind=kind_phys), intent(inout) :: clubb_c6rt
    real(kind=kind_phys), intent(inout) :: clubb_c6rtb
    real(kind=kind_phys), intent(inout) :: clubb_c6rtc
    real(kind=kind_phys), intent(inout) :: clubb_c6thl
    real(kind=kind_phys), intent(inout) :: clubb_c6thlb
    real(kind=kind_phys), intent(inout) :: clubb_c6thlc
    real(kind=kind_phys), intent(inout) :: clubb_wpxp_L_thresh
    real(kind=kind_phys), intent(inout) :: clubb_C7
    real(kind=kind_phys), intent(inout) :: clubb_C7b
    real(kind=kind_phys), intent(inout) :: clubb_gamma_coef
    real(kind=kind_phys), intent(inout) :: clubb_c_K10h
    real(kind=kind_phys), intent(inout) :: clubb_lambda0_stability_coef
    real(kind=kind_phys), intent(inout) :: clubb_lmin_coef
    real(kind=kind_phys), intent(inout) :: clubb_C8b
    real(kind=kind_phys), intent(inout) :: clubb_skw_max_mag
    real(kind=kind_phys), intent(inout) :: clubb_C1
    real(kind=kind_phys), intent(inout) :: clubb_C1b
    real(kind=kind_phys), intent(inout) :: clubb_gamma_coefb
    real(kind=kind_phys), intent(inout) :: clubb_up2_sfc_coef
    real(kind=kind_phys), intent(inout) :: clubb_C4
    real(kind=kind_phys), intent(inout) :: clubb_C_uu_shr
    real(kind=kind_phys), intent(inout) :: clubb_C_uu_buoy
    real(kind=kind_phys), intent(inout) :: clubb_c_K1
    real(kind=kind_phys), intent(inout) :: clubb_c_K2
    real(kind=kind_phys), intent(inout) :: clubb_nu2
    real(kind=kind_phys), intent(inout) :: clubb_c_K8
    real(kind=kind_phys), intent(inout) :: clubb_c_K9
    real(kind=kind_phys), intent(inout) :: clubb_nu9
    real(kind=kind_phys), intent(inout) :: clubb_C_wp2_splat
    real(kind=kind_phys), intent(inout) :: clubb_C_invrs_tau_bkgnd
    real(kind=kind_phys), intent(inout) :: clubb_C_invrs_tau_sfc
    real(kind=kind_phys), intent(inout) :: clubb_C_invrs_tau_shear
    real(kind=kind_phys), intent(inout) :: clubb_C_invrs_tau_N2
    real(kind=kind_phys), intent(inout) :: clubb_C_invrs_tau_N2_wp2
    real(kind=kind_phys), intent(inout) :: clubb_C_invrs_tau_N2_xp2
    real(kind=kind_phys), intent(inout) :: clubb_C_invrs_tau_N2_wpxp
    real(kind=kind_phys), intent(inout) :: clubb_C_invrs_tau_N2_clear_wp3
    real(kind=kind_phys), intent(inout) :: clubb_bv_efold
    real(kind=kind_phys), intent(inout) :: clubb_wpxp_Ri_exp
    real(kind=kind_phys), intent(inout) :: clubb_z_displace

    !----- Begin Code -----

    clubb_init_errcode = 0 

    if (core_rknd /= kind_phys) then
      error_message = 'clubb_ini_cam:  CLUBB library core_rknd must match CAM kind_phys and it does not'
      clubb_init_errcode = 1
      return
    end if

    ! Allocate PDF parameters across columns and chunks
    allocate( &
       pdf_params_chnk(begchunk:endchunk),   &
       pdf_params_zm_chnk(begchunk:endchunk), &
       pdf_implicit_coefs_terms_chnk(begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
      error_message = ' clubb_ini_cam: failed to allocate pdf_params'
      clubb_init_errcode = 1
      return
    end if

    ! ----------------------------------------------------------------- !
    ! Determine how many constituents CLUBB will transport.  Note that
    ! CLUBB does not transport aerosol consituents.  Therefore, need to
    ! determine how many aerosols constituents there are and subtract that
    ! off of pcnst (the total consituents)
    ! ----------------------------------------------------------------- !

    !  Select variables to apply tendencies back to CAM

    ! Initialize all consituents to true to start
    lq(1:pcnst) = .true.
    edsclr_dim  = pcnst

    do m = 1, pcnst
       if (cnst_ndropmixed(m)) then
          lq(m)=.false.
          !  Droplet number is transported in dropmixnuc, therefore we
          !  do NOT want CLUBB to apply transport tendencies to avoid double
          !  counting.  Else, we apply tendencies.
          edsclr_dim = edsclr_dim-1
       endif
    enddo

    ! ----------------------------------------------------------------- !
    ! Set the debug level.  Level 2 has additional computational expense since
    ! it checks the array variables in CLUBB for invalid values.
    ! ----------------------------------------------------------------- !
    call set_clubb_debug_level_api( 0 )

    ! ----------------------------------------------------------------- !
    ! use pbuf_get_fld_idx to get existing physics buffer fields from other
    ! physics packages (e.g. tke)
    ! ----------------------------------------------------------------- !

    !  Defaults
    stats_metadata%l_stats_samp = .false.
    stats_metadata%l_grads = .false.

    !  Overwrite defaults if needed
    if (stats_metadata%l_stats) stats_metadata%l_stats_samp = .true.

    ! Scalars aren't in use, set all indices to -1
    sclr_idx%iisclr_rt  = -1
    sclr_idx%iisclr_thl = -1
    sclr_idx%iisclr_CO2 = -1
    sclr_idx%iiedsclr_rt  = -1
    sclr_idx%iiedsclr_thl = -1
    sclr_idx%iiedsclr_CO2 = -1

    ! ----------------------------------------------------------------- !
    ! Define number of tracers for CLUBB to diffuse
    ! ----------------------------------------------------------------- !

    if (clubb_l_do_expldiff_rtm_thlm) then
      ! add 2 since we want to diffuse temperature and moisture explicitly as well
      edsclr_dim = edsclr_dim + 2
    endif

    ! ----------------------------------------------------------------- !
    ! Setup CLUBB core
    ! ----------------------------------------------------------------- !

     call init_clubb_params_api( 1, -99, "", &
                                     clubb_params_single_col )

    clubb_params_single_col(1,iC2rtthl)                       = clubb_C2rtthl
    clubb_params_single_col(1,iC8)                            = clubb_C8
    clubb_params_single_col(1,iC11)                           = clubb_c11
    clubb_params_single_col(1,iC11b)                          = clubb_c11b
    clubb_params_single_col(1,iC14)                           = clubb_c14
    clubb_params_single_col(1,iC_wp3_pr_turb)                 = clubb_C_wp3_pr_turb
    clubb_params_single_col(1,ic_K10)                         = clubb_c_K10
    clubb_params_single_col(1,imult_coef)                     = clubb_mult_coef
    clubb_params_single_col(1,iSkw_denom_coef)                = clubb_Skw_denom_coef
    clubb_params_single_col(1,iC2rt)                          = clubb_C2rt
    clubb_params_single_col(1,iC2thl)                         = clubb_C2thl
    clubb_params_single_col(1,ibeta)                          = clubb_beta
    clubb_params_single_col(1,iC6rt)                          = clubb_c6rt
    clubb_params_single_col(1,iC6rtb)                         = clubb_c6rtb
    clubb_params_single_col(1,iC6rtc)                         = clubb_c6rtc
    clubb_params_single_col(1,iC6thl)                         = clubb_c6thl
    clubb_params_single_col(1,iC6thlb)                        = clubb_c6thlb
    clubb_params_single_col(1,iC6thlc)                        = clubb_c6thlc
    clubb_params_single_col(1,iwpxp_L_thresh)                 = clubb_wpxp_L_thresh
    clubb_params_single_col(1,iC7)                            = clubb_C7
    clubb_params_single_col(1,iC7b)                           = clubb_C7b
    clubb_params_single_col(1,igamma_coef)                    = clubb_gamma_coef
    clubb_params_single_col(1,ic_K10h)                        = clubb_c_K10h
    clubb_params_single_col(1,ilambda0_stability_coef)        = clubb_lambda0_stability_coef
    clubb_params_single_col(1,ilmin_coef)                     = clubb_lmin_coef
    clubb_params_single_col(1,iC8b)                           = clubb_C8b
    clubb_params_single_col(1,iskw_max_mag)                   = clubb_skw_max_mag
    clubb_params_single_col(1,iC1)                            = clubb_C1
    clubb_params_single_col(1,iC1b)                           = clubb_C1b
    clubb_params_single_col(1,igamma_coefb)                   = clubb_gamma_coefb
    clubb_params_single_col(1,iup2_sfc_coef)                  = clubb_up2_sfc_coef
    clubb_params_single_col(1,iC4)                            = clubb_C4
    clubb_params_single_col(1,iC_uu_shr)                      = clubb_C_uu_shr
    clubb_params_single_col(1,iC_uu_buoy)                     = clubb_C_uu_buoy
    clubb_params_single_col(1,ic_K1)                          = clubb_c_K1
    clubb_params_single_col(1,ic_K2)                          = clubb_c_K2
    clubb_params_single_col(1,inu2)                           = clubb_nu2
    clubb_params_single_col(1,ic_K8)                          = clubb_c_K8
    clubb_params_single_col(1,ic_K9)                          = clubb_c_K9
    clubb_params_single_col(1,inu9)                           = clubb_nu9
    clubb_params_single_col(1,iC_wp2_splat)                   = clubb_C_wp2_splat
    clubb_params_single_col(1,iC_invrs_tau_bkgnd)             = clubb_C_invrs_tau_bkgnd
    clubb_params_single_col(1,iC_invrs_tau_sfc)               = clubb_C_invrs_tau_sfc
    clubb_params_single_col(1,iC_invrs_tau_shear)             = clubb_C_invrs_tau_shear
    clubb_params_single_col(1,iC_invrs_tau_N2)                = clubb_C_invrs_tau_N2
    clubb_params_single_col(1,iC_invrs_tau_N2_wp2)            = clubb_C_invrs_tau_N2_wp2
    clubb_params_single_col(1,iC_invrs_tau_N2_xp2)            = clubb_C_invrs_tau_N2_xp2
    clubb_params_single_col(1,iC_invrs_tau_N2_wpxp)           = clubb_C_invrs_tau_N2_wpxp
    clubb_params_single_col(1,iC_invrs_tau_N2_clear_wp3)      = clubb_C_invrs_tau_N2_clear_wp3
    clubb_params_single_col(1,ibv_efold)                      = clubb_bv_efold
    clubb_params_single_col(1,iwpxp_Ri_exp)                   = clubb_wpxp_Ri_exp
    clubb_params_single_col(1,iz_displace)                    = clubb_z_displace

    ! Override clubb default
    if ( trim(subcol_scheme) == 'SILHS' ) then
      clubb_config_flags%saturation_formula = saturation_flatau
    else
      clubb_config_flags%saturation_formula = saturation_gfdl     ! Goff & Gratch (1946) approximation for SVP
    end if

    !  Set up CLUBB core.  Note that some of these inputs are overwritten
    !  when clubb_tend_cam is called.  The reason is that heights can change
    !  at each time step, which is why dummy arrays are read in here for heights
    !  as they are immediately overwrote.
    !! Initialize err_info with default values since info is not available here
    call init_default_err_info_api(1, err_info)
!$OMP PARALLEL
    call check_clubb_settings_api( 1, clubb_params_single_col,  & ! Intent(in)
                                   l_implemented,               & ! Intent(in)
                                   l_input_fields,              & ! Intent(in)
                                   clubb_config_flags,          & ! intent(in)
                                   err_info )                     ! Intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) ) then
       error_message = 'clubb_ini_cam: FATAL ERROR CALLING CHECK_CLUBB_SETTINGS_API'
       clubb_init_errcode = 1
       return
    end if
!$OMP END PARALLEL

    ! Cleanup err_info since it is not needed anymore
    call cleanup_err_info_api(err_info)

    ! Print the list of CLUBB parameters
    if ( masterproc ) then ! masterproc and iulog need to be passed in.
       do j = 1, nparams, 1
          write(iulog,*) params_list(j), " = ", clubb_params_single_col(1,j)
       enddo
    endif

    ! Print configurable CLUBB flags
    if ( masterproc ) then
       write(iulog,'(a,i0,a)') " CLUBB configurable flags "
       call print_clubb_config_flags_api( iulog, clubb_config_flags ) ! Intent(in)
    end if

    if ( trim(subcol_scheme) /= 'SILHS' ) then
       ! hm_metadata is set up by calling init_pdf_hydromet_arrays_api in subcol_init_SILHS.
       ! So if we are not using silhs, we allocate the parts of hm_metadata that need allocating
       ! in order to making intel debug tests happy.
       allocate( hm_metadata%hydromet_list(1), stat=ierr)
       if( ierr /= 0 ) then
         error_message = 'clubb_ini_cam: Unable to allocate hm_metadata%hydromet_list'
         clubb_init_errcode = 1
         return
       end if
       allocate( hm_metadata%l_mix_rat_hm(1), stat=ierr)
       if( ierr /= 0 ) then
         error_message = 'clubb_ini_cam: Unable to allocate hm_metadata%l_mix_rat_hm'
         clubb_init_errcode = 1
         return
       end if
    end if

    !  Initialize statistics, below are dummy variables
    dum1 = 300._kind_phys
    dum2 = 1200._kind_phys
    dum3 = 300._kind_phys

    if (stats_metadata%l_stats) then

      call stats_init_clubb( pcols, masterproc, mpicom, mstrid, mpi_character, &
                             .true., dum1, dum2, max_fieldname_len, &
                             nzm_clubb, nzt_clubb, nzm_clubb, dum3, &
                             edsclr_dim, sclr_dim, hydromet_dim, &
                             hm_metadata, stats_metadata, &
                             clubb_vars_zt, clubb_vars_zm, clubb_vars_sfc, &
                             clubb_vars_rad_zt, clubb_vars_rad_zm, &
                             stats_zt(:), stats_zm(:), stats_sfc(:), &
                             stats_rad_zt(:), stats_rad_zm(:), &
                             error_message, clubb_init_errcode)

      if (clubb_init_errcode /= 0) return 

       allocate(out_zt(pcols,pver,stats_zt(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) then 
         error_message = 'clubb_ini_cam: Unable to allocate out_zt'
         clubb_init_errcode = 1
         return
       end if
       allocate(out_zm(pcols,pverp,stats_zm(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) then 
         error_message = 'clubb_ini_cam: Unable to allocate out_zm'
         clubb_init_errcode = 1
         return
       end if
       allocate(out_sfc(pcols,1,stats_sfc(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) then
         error_message = 'clubb_ini_cam: Unable to allocate out_sfc'
         clubb_init_errcode = 1
         return
       end if

       if ( stats_metadata%l_output_rad_files ) then
          allocate(out_radzt(pcols,pver,stats_rad_zt(1)%num_output_fields), stat=ierr)
          if( ierr /= 0 ) then
            error_message = 'clubb_ini_cam: Unable to allocate out_radzt'
            clubb_init_errcode = 1
            return
          end if
          allocate(out_radzm(pcols,pverp,stats_rad_zm(1)%num_output_fields), stat=ierr)
          if( ierr /= 0 ) then
            error_message = 'clubb_ini_cam: Unable to allocate out_radzm'
            clubb_init_errcode = 1
            return
          end if
       end if

    endif

    end subroutine clubb_init


  subroutine clubb2_run(ncol, pver, meltpt_temp, latice, rga, &
                        ixcldliq, ixcldice, ixnumliq, ixnumice, &
                        dlf, dlf_liq_out, dlf_ice_out, &
                        clubb_detliq_rad, clubb_detice_rad, clubb_detphase_lowtemp, &
                        pdel, pdeldry, &
                        s, t, q, det_s, det_ice )

    ! Incoming variables
    integer, intent(in) :: ncol, pver, ixcldliq, ixcldice, ixnumliq, ixnumice
    real(kind_phys), intent(in) :: clubb_detliq_rad, clubb_detice_rad, clubb_detphase_lowtemp
    real(kind_phys), intent(in) :: meltpt_temp, latice, rga
    real(kind_phys), intent(in) :: dlf(:,:)
    real(kind_phys), intent(in) :: t(:,:), pdel(:,:), pdeldry(:,:)
    real(kind_phys), intent(inout) :: q(:,:,:)
    real(kind_phys), intent(inout) :: s(:,:)
    real(kind_phys), intent(inout) :: det_s(:), det_ice(:)
    real(kind_phys), intent(out) :: dlf_liq_out(:,:), dlf_ice_out(:,:)

    ! Local variables
    integer :: i, k
    real(kind_phys) :: dlf2, dum1, dl_rad, di_rad, dt_low

    ! ------------------------------------------------------------ !
    ! The rest of the code deals with diagnosing variables         !
    ! for microphysics/radiation computation and macrophysics      !
    ! ------------------------------------------------------------ !

    ! --------------------------------------------------------------------------------- !
    !  COMPUTE THE ICE CLOUD DETRAINMENT                                                !
    !  Detrainment of convective condensate into the environment or stratiform cloud    !
    ! --------------------------------------------------------------------------------- !

    !  Initialize the shallow convective detrainment rate, will always be zero
    dlf2 = 0.0_kind_phys
    dlf_liq_out(:,:) = 0.0_kind_phys
    dlf_ice_out(:,:) = 0.0_kind_phys

    dl_rad = clubb_detliq_rad
    di_rad = clubb_detice_rad
    dt_low = clubb_detphase_lowtemp

    do k = 1, pver
      do i = 1, ncol

        if( t(i,k) > meltpt_temp ) then
          dum1 = 0.0_kind_phys
        elseif ( t(i,k) < dt_low ) then
          dum1 = 1.0_kind_phys
        else
          dum1 = ( meltpt_temp - t(i,k) ) / ( meltpt_temp - dt_low )
        endif

        q(i,k,ixcldliq) = dlf(i,k) * ( 1._kind_phys - dum1 )
        q(i,k,ixcldice) = dlf(i,k) * dum1
        q(i,k,ixnumliq) = 3._kind_phys * ( max(0._kind_phys, ( dlf(i,k) - dlf2 )) * ( 1._kind_phys - dum1 ) ) &
                                   / (4._kind_phys*3.14_kind_phys*dl_rad**3*997._kind_phys) + & ! Deep    Convection
                                   3._kind_phys * (                         dlf2    * ( 1._kind_phys - dum1 ) ) &
                                   / (4._kind_phys*3.14_kind_phys*10.e-6_kind_phys**3*997._kind_phys)     ! Shallow Convection
        q(i,k,ixnumice) = 3._kind_phys * ( max(0._kind_phys, ( dlf(i,k) - dlf2 )) *  dum1 ) &
                                   / (4._kind_phys*3.14_kind_phys*di_rad**3*500._kind_phys) + & ! Deep    Convection
                                   3._kind_phys * (                         dlf2    *  dum1 ) &
                                   / (4._kind_phys*3.14_kind_phys*50.e-6_kind_phys**3*500._kind_phys)     ! Shallow Convection
        s(i,k)          = dlf(i,k) * dum1 * latice

        dlf_liq_out(i,k) = dlf(i,k) * ( 1._kind_phys - dum1 )
        dlf_ice_out(i,k) = dlf(i,k) * dum1

        ! convert moist dlf tendencies to dry
        q(i,k,ixcldliq) = q(i,k,ixcldliq)*pdel(i,k)/pdeldry(i,k)
        q(i,k,ixcldice) = q(i,k,ixcldice)*pdel(i,k)/pdeldry(i,k)

        ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
        !   track of the integrals of ice and static energy that is effected from conversion to ice
        !   so that the energy checker doesn't complain.
        det_s(i)                  = det_s(i)   + s(i,k)          * pdel(i,k)    * rga
        det_ice(i)                = det_ice(i) - q(i,k,ixcldice) * pdeldry(i,k) * rga
      enddo
    enddo

    det_ice(:ncol) = det_ice(:ncol) / 1000._kind_phys  ! divide by density of water


  end subroutine clubb2_run


  subroutine clubb3_run(ncol, pver, pverp, pcnst, top_lev, zvir, rair, cpair, gravit, karman, &
                        ixq, ixcldice, ixcldliq, ixnumice, calday, tropp_days, tropLev, &
                        rhminis_const, rhmaxis_const, rhmini_const, rhmaxi_const, &
                        single_column, scm_cambfb_mode, scm_clubb_iop_name, subcol_scheme, &
                        dp1, dp2, cmfmc, cmfmc_sh_pbuf, dp_icwmr_pbuf, concld_pbuf, &
                        aist_pbuf, qsatfac_pbuf, ast_pbuf, qist_pbuf, cld_pbuf, &
                        pblh_pbuf, deepcu_pbuf, shalcu_pbuf, &
                        lq, cnst_type, &
                        alst_pbuf, qlst_pbuf, rcm, cloud_frac, exner, &
                        t, state_q, ptend_q, pmid, landfrac, snowhland, pdel, pdeldry, &
                        wsx, wsy, shf, cflx, zm, zi, u, v, lat, pint, phis, &
                        errmsg, errflg )

!    use tropopause_find,       only: tropopause_findChemTrop
    use holtslag_boville_diff, only: hb_pbl_dependent_coefficients_run
!    use cldfrc2m,              only: aist_vector
    use atmos_phys_pbl_utils,  only: calc_friction_velocity, calc_obukhov_length, calc_ideal_gas_rrho, &
                                     calc_kinematic_heat_flux, calc_kinematic_water_vapor_flux, &
                                     calc_kinematic_buoyancy_flux

    ! Incoming variables
    integer, intent(in) :: ncol, pver, pverp, pcnst, top_lev, ixq, ixcldice, ixcldliq, ixnumice
    logical, intent(in) :: single_column, scm_cambfb_mode
    logical, intent(in) :: lq(:)
    character(len=3), intent(in) :: cnst_type(:)
    real(kind_phys), intent(in) :: dp1, dp2, zvir, rair, cpair, gravit, karman, calday
    real(kind_phys), intent(in) :: rhminis_const, rhmaxis_const, rhmini_const, rhmaxi_const
    real(kind_phys), intent(in) :: lat(:), phis(:)
    real(kind_phys), intent(in) :: pint(:,:)
    real(kind_phys), intent(in) :: pmid(:,:)
    real(kind_phys), intent(in) :: landfrac(:), snowhland(:)
    real(kind_phys), intent(in) :: pdel(:,:), pdeldry(:,:), rcm(:,:), cloud_frac(:,:), &
                                   t(:,:), exner(:,:), zm(:,:), zi(:,:), u(:,:), v(:,:), cmfmc(:,:)
    real(kind_phys), intent(in) :: state_q(:,:,:)
    real(kind_phys), intent(in) :: wsx(:), wsy(:), shf(:)
    real(kind_phys), intent(in) :: cflx(:,:)
    real(kind_phys), intent(inout) :: ptend_q(:,:,:)
    real(kind_phys), intent(inout) :: pblh_pbuf(:)
    real(kind_phys), intent(inout) :: alst_pbuf(:,:), qlst_pbuf(:,:), deepcu_pbuf(:,:), shalcu_pbuf(:,:), &
                                      cmfmc_sh_pbuf(:,:), dp_icwmr_pbuf(:,:), concld_pbuf(:,:), &
                                      aist_pbuf(:,:), qsatfac_pbuf(:,:), ast_pbuf(:,:), qist_pbuf(:,:), &
                                      cld_pbuf(:,:)
    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    !real(kind_phys), intent(in)         :: tropp_p_loc(:,:)
    real(kind_phys), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12
    character(len=20), intent(in) :: scm_clubb_iop_name
    character(len=16), intent(in) :: subcol_scheme
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Local variables
    integer :: i, k, ixind, k_clubb, rhmini, rhmaxi
    real(kind_phys) :: frac_limit, ic_limit
    real(kind_phys) :: rrho(ncol), ustar2(ncol), kinheat(ncol), kinwat(ncol), kbfs(ncol), obklen(ncol), &
                       dummy2(ncol), dummy3(ncol)
    real(kind_phys) :: th(ncol,pver), thv(ncol,pver)
    integer, intent(in) :: tropLev(:)

    ! ---------------------------------------------------------

    errmsg = ''
    errflg = 0

    ! ptend_all now has all accumulated tendencies.  Convert the tendencies for the
    ! wet constituents to wet air basis.
    do ixind = 1, pcnst
      if (lq(ixind) .and. cnst_type(ixind) == 'wet') then
        do k = 1, pver
          do i = 1, ncol
            ptend_q(i,k,ixind) = ptend_q(i,k,ixind)*pdeldry(i,k)/pdel(i,k)
          end do
        end do
      end if
    end do

    ! --------------------------------------------------------------------------------- !
    !  Diagnose some quantities that are computed in macrop_tend here.                  !
    !  These are inputs required for the microphysics calculation.                      !
    !                                                                                   !
    !  FIRST PART COMPUTES THE STRATIFORM CLOUD FRACTION FROM CLUBB CLOUD FRACTION      !
    ! --------------------------------------------------------------------------------- !

    !  initialize variables
    alst_pbuf(:,:) = 0.0_kind_phys
    qlst_pbuf(:,:) = 0.0_kind_phys

    do k = top_lev, pver
      do i = 1, ncol
        k_clubb     = k + 1 - top_lev
        alst_pbuf(i,k)    = cloud_frac(i,k_clubb)
        qlst_pbuf(i,k)    = rcm(i,k_clubb) / max( 0.01_kind_phys, alst_pbuf(i,k) )  ! Incloud stratus condensate mixing ratio
      enddo
    enddo

    ! --------------------------------------------------------------------------------- !
    !  THIS PART COMPUTES CONVECTIVE AND DEEP CONVECTIVE CLOUD FRACTION                 !
    ! --------------------------------------------------------------------------------- !

    frac_limit = 0.01_kind_phys
    ic_limit   = 1.e-12_kind_phys
    deepcu_pbuf(:,:) = 0.0_kind_phys
    shalcu_pbuf(:,:) = 0.0_kind_phys

    do k = 1, pver-1
      do i = 1, ncol
        !  diagnose the deep convective cloud fraction, as done in macrophysics based on the
        !  deep convective mass flux, read in from pbuf.  Since shallow convection is never
        !  called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
        !  fraction is purely from deep convection scheme.
        deepcu_pbuf(i,k) = max(0.0_kind_phys,min(dp1*log(1.0_kind_phys+dp2*(cmfmc(i,k+1)-cmfmc_sh_pbuf(i,k+1))),0.6_kind_phys))

        if (deepcu_pbuf(i,k) <= frac_limit .or. dp_icwmr_pbuf(i,k) < ic_limit) then
          deepcu_pbuf(i,k) = 0._kind_phys
        endif

        !  using the deep convective cloud fraction, and CLUBB cloud fraction (variable
        !  "cloud_frac"), compute the convective cloud fraction.  This follows the formulation
        !  found in macrophysics code.  Assumes that convective cloud is all nonstratiform cloud
        !  from CLUBB plus the deep convective cloud fraction
        ! NOTE: concld_pbuf used to be calculated in the commented-out version below, but since we
        ! set alst_pbuf=cloud_frac_pbuf, this simplifies to only using deepcu_pbuf.
        ! This is potentially a bug, but there's not really a "right" way to combine the different
        ! cloud factions, so it has been left to only use deepcu_pbuf for now
        !concld_pbuf(i,k) = min(cloud_frac_pbuf(i,k)-alst_pbuf(i,k)+deepcu_pbuf(i,k),0.80_kind_phys)
        concld_pbuf(i,k) = min(deepcu_pbuf(i,k),0.80_kind_phys)
      enddo
    enddo

    if (single_column .and. .not. scm_cambfb_mode) then
      if (trim(scm_clubb_iop_name)  ==  'ATEX_48hr'       .or. &
          trim(scm_clubb_iop_name)  ==  'BOMEX_5day'      .or. &
          trim(scm_clubb_iop_name)  ==  'DYCOMSrf01_4day' .or. &
          trim(scm_clubb_iop_name)  ==  'DYCOMSrf02_06hr' .or. &
          trim(scm_clubb_iop_name)  ==  'RICO_3day'       .or. &
          trim(scm_clubb_iop_name)  ==  'ARM_CC') then

         deepcu_pbuf(:,:) = 0.0_kind_phys
         concld_pbuf(:,:) = 0.0_kind_phys

      endif
    endif

    ! --------------------------------------------------------------------------------- !
    !  COMPUTE THE ICE CLOUD FRACTION PORTION                                           !
    !  use the aist_vector function to compute the ice cloud fraction                   !
    ! --------------------------------------------------------------------------------- !

    !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
!    troplev(:) = 0
    !REMOVECAM_END
!    call tropopause_findChemTrop( state, troplev )
!    call tropopause_findChemTrop(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
!                                 calday, tropp_p_loc, tropp_days, &
!                                 tropLev, & 
!                                 errmsg, errflg)

    if (errflg /= 0) return

    aist_pbuf(:,:top_lev-1) = 0._kind_phys
    qsatfac_pbuf(:, :) = 0._kind_phys ! Zero out entire profile in case qsatfac is left undefined in aist_vector below

!    do k = top_lev, pver
!
!      ! For Type II PSC and for thin cirrus, the clouds can be thin, but
!      ! extensive and they should start forming when the gridbox mean saturation
!      ! reaches 1.0.
!      !
!      ! For now, use the tropopause diagnostic to determine where the Type II
!      ! PSC should be, but in the future wold like a better metric that can also
!      ! identify the level for thin cirrus. Include the tropopause level so that
!      ! the cold point tropopause will use the stratospheric values.
!      where (k <= troplev)
!        rhmini = rhminis_const
!        rhmaxi = rhmaxis_const
!      elsewhere
!        rhmini = rhmini_const
!        rhmaxi = rhmaxi_const
!      end where
!
!      if ( trim(subcol_scheme) == 'SILHS' ) then
!        call aist_vector(state_q(:,k,ixq),t(:,k),pmid(:,k),state_q(:,k,ixcldice), &
!             state_q(:,k,ixnumice), landfrac(:),snowhland(:),aist_pbuf(:,k),ncol )
!      else
!        call aist_vector(state_q(:,k,ixq),t(:,k),pmid(:,k),state_q(:,k,ixcldice), &
!              state_q(:,k,ixnumice), landfrac(:),snowhland(:),aist_pbuf(:,k),ncol,&
!              qsatfac_out=qsatfac_pbuf(:,k), rhmini_in=rhmini, rhmaxi_in=rhmaxi)
!      endif
!    enddo

    ! --------------------------------------------------------------------------------- !
    !  THIS PART COMPUTES THE LIQUID STRATUS FRACTION                                   !
    !                                                                                   !
    !  For now leave the computation of ice stratus fraction from macrop_driver intact  !
    !  because CLUBB does nothing with ice.  Here I simply overwrite the liquid stratus !
    !  fraction that was coded in macrop_driver                                         !
    ! --------------------------------------------------------------------------------- !

    do k = 1, pver
      do i = 1, ncol

        !  Recompute net stratus fraction using maximum over-lapping assumption, as done
        !  in macrophysics code, using alst computed above and aist read in from physics buffer
        ast_pbuf(i,k) = max(alst_pbuf(i,k),aist_pbuf(i,k))
        qist_pbuf(i,k) = state_q(i,k,ixcldice)/max(0.01_kind_phys,aist_pbuf(i,k))

        !  Probably need to add deepcu cloud fraction to the cloud fraction array, else would just
        !  be outputting the shallow convective cloud fraction
        cld_pbuf(i,k) = min(ast_pbuf(i,k)+deepcu_pbuf(i,k),1.0_kind_phys)

      enddo
    enddo

    ! --------------------------------------------------------------------------------- !
    !  DIAGNOSE THE PBL DEPTH                                                           !
    !  this is needed for aerosol code                                                  !
    ! --------------------------------------------------------------------------------- !
    do k = 1, pver
      do i = 1, ncol
         !subroutine pblind expects "Stull" definition of Exner
         th(i,k) = t(i,k)*exner(i,k)
         !thv should have condensate loading to be consistent with earlier def's in this module
         thv(i,k) = th(i,k)*(1.0_kind_phys+zvir*state_q(i,k,ixq) - state_q(i,k,ixcldliq))
      enddo
    enddo

    ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)
    rrho   (1:ncol) = calc_ideal_gas_rrho(rair, t(1:ncol,pver), pmid(1:ncol,pver))
    ustar2 (1:ncol) = calc_friction_velocity(wsx(1:ncol), wsy(1:ncol), rrho(1:ncol))
    ! use correct qflux from coupler
    kinheat(1:ncol) = calc_kinematic_heat_flux(shf(1:ncol), rrho(1:ncol), cpair)
    kinwat (1:ncol) = calc_kinematic_water_vapor_flux(cflx(1:ncol,1), rrho(1:ncol))
    kbfs   (1:ncol) = calc_kinematic_buoyancy_flux(kinheat(1:ncol), zvir, th(1:ncol,pver), kinwat(1:ncol))
    obklen (1:ncol) = calc_obukhov_length(thv(1:ncol,pver), ustar2(1:ncol), gravit, karman, kbfs(1:ncol))

    where (kbfs(:ncol)  ==  -0.0_kind_phys) kbfs(:ncol) = 0.0_kind_phys

    ! Compute PBL depth according to Holtslag-Boville Scheme -- only pblh is needed here
    ! and other outputs are discarded
    !REMOVECAM - no longer need this when CAM is retired and pcols no longer exists
    pblh_pbuf(:) = 0._kind_phys
    dummy2(:) = 0._kind_phys
    dummy3(:) = 0._kind_phys
    !REMOVECAM_END
    call hb_pbl_dependent_coefficients_run( &
      ncol      = ncol,                                      &
      pver      = pver,                                      &
      pverp     = pverp,                                     &
      gravit    = gravit,                                    &
      z         = zm(:ncol,:pver),                    &
      zi        = zi(:ncol,:pverp),                   &
      u         = u(:ncol,:pver),                     &
      v         = v(:ncol,:pver),                     &
      cldn      = cld_pbuf(:ncol,:pver),                   &
      ! Inputs from CLUBB (not HB coefficients)
      thv       = thv(:ncol,:pver),                          &
      ustar     = ustar2(:ncol),                             &
      kbfs      = kbfs(:ncol),                               &
      obklen    = obklen(:ncol),                             &
      ! Output variables
      pblh      = pblh_pbuf(:ncol),                               &
      wstar     = dummy2(:ncol),                             &
      bge       = dummy3(:ncol),                             &
      errmsg    = errmsg,                                    &
      errflg    = errflg)

    if (errflg /= 0) return
    ! --------------------------------------------------------------------------------- !
    !                              END CLOUD FRACTION DIAGNOSIS                         !
    ! --------------------------------------------------------------------------------- !

  end subroutine clubb3_run


#ifdef CLUBB_SGS

  subroutine stats_init_clubb( pcols, masterproc, mpicom, mstrid, mpi_character, &
                               l_stats_in, stats_tsamp_in, stats_tout_in, max_fieldname_len, &
                               nnzp, nnrad_zt,nnrad_zm, delt, &
                               edsclr_dim, sclr_dim, hydromet_dim, &
                               hm_metadata, stats_metadata, &
                               clubb_vars_zt, clubb_vars_zm, clubb_vars_sfc, &
                               clubb_vars_rad_zt, clubb_vars_rad_zm, &
                               stats_zt, stats_zm, stats_sfc, &
                               stats_rad_zt, stats_rad_zm, &
                               error_message, clubb_init_errcode )
    !
    ! Description: Initializes the statistics saving functionality of
    !   the CLUBB model.  This is for purpose of CAM-CLUBB interface.  Here
    !   the traditional stats_init of CLUBB is not called, as it is not compatible
    !   with CAM output.

    !-----------------------------------------------------------------------

    use clubb_api_module, only:        time_precision, stats, stats_metadata_type, hm_metadata_type, &   !
                                       nvarmax_zm, stats_init_zm_api, & !
                                       nvarmax_zt, stats_init_zt_api, & !
                                       nvarmax_rad_zt, stats_init_rad_zt_api, & !
                                       nvarmax_rad_zm, stats_init_rad_zm_api, & !
                                       nvarmax_sfc, stats_init_sfc_api, & !
                                       fstderr, var_length !
    use cam_history,            only: addfld, horiz_only

    implicit none

    !----------------------- Input Variables -----------------------

    logical, intent(in) :: l_stats_in ! Stats on? T/F
  
    logical, intent(in) :: masterproc

    integer, intent(in) :: mpicom, mpi_character, mstrid, pcols, max_fieldname_len

    integer, intent(in) :: sclr_dim, edsclr_dim, hydromet_dim

    integer, intent(inout) :: clubb_init_errcode
    character(len=200), intent(inout) :: error_message

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type (hm_metadata_type), intent(inout) :: &
      hm_metadata

    real(kind=time_precision), intent(in) ::  &
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    integer, intent(in) :: nnzp     ! Grid points in the vertical [count]
    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]
    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real(kind=time_precision), intent(in) ::   delt         ! Timestep (dtmain in CLUBB)         [s]

    !----------------------- Output Variables -----------------------
    type (stats), intent(out) :: &
      stats_zt(:),      & ! stats_zt grid
      stats_zm(:),      & ! stats_zm grid
      stats_rad_zt(:),  & ! stats_rad_zt grid
      stats_rad_zm(:),  & ! stats_rad_zm grid
      stats_sfc(:)        ! stats_sfc


    !----------------------- Local Variables -----------------------

    !  Namelist Variables

    character(len=var_length), intent(in)     ::   clubb_vars_zt(:)      ! Variables on the thermodynamic levels
    character(len=var_length), intent(in)     ::   clubb_vars_zm(:)      ! Variables on the momentum levels
    character(len=var_length), intent(in) ::   clubb_vars_rad_zt(:)  ! Variables on the radiation levels
    character(len=var_length), intent(in) ::   clubb_vars_rad_zm(:)  ! Variables on the radiation levels
    character(len=var_length), intent(in)    ::   clubb_vars_sfc(:)     ! Variables at the model surface

    logical :: l_error

    character(len=200) :: temp1, sub

    integer :: i, ntot, j
    integer :: ierr

    !----------------------- Begin Code -----------------------

    !  Initialize
    l_error = .false.

    !  Set stats_variables variables with inputs from calling subroutine
    stats_metadata%l_stats = l_stats_in

    stats_metadata%stats_tsamp = stats_tsamp_in
    stats_metadata%stats_tout  = stats_tout_in

    if ( .not. stats_metadata%l_stats ) then
       stats_metadata%l_stats_samp  = .false.
       stats_metadata%l_stats_last  = .false.
       return
    end if

    !  Hardcode these for use in CAM-CLUBB, don't want either
    stats_metadata%l_netcdf = .false.
    stats_metadata%l_grads  = .false.

    !  Check sampling and output frequencies
    do j = 1, pcols

      !  The model time step length, delt (which is dtmain), should multiply
      !  evenly into the statistical sampling time step length, stats_tsamp.
      if ( abs( stats_metadata%stats_tsamp/delt - floor(stats_metadata%stats_tsamp/delt) ) > 1.e-8_kind_phys ) then
         l_error = .true.  ! This will cause the run to stop.
         write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                          'the clubb time step (delt below)'
         write(fstderr,*) 'stats_tsamp = ', stats_metadata%stats_tsamp
         write(fstderr,*) 'delt = ', delt
         error_message = "stats_init_clubb:  CLUBB stats_tsamp must be an even multiple of the timestep"
         clubb_init_errcode = 1
         return
      endif

      !  Initialize zt (mass points)

      i = 1
      do while ( ichar(clubb_vars_zt(i)(1:1)) /= 0 .and. &
                 len_trim(clubb_vars_zt(i))   /= 0 .and. &
                 i <= nvarmax_zt )
         i = i + 1
      enddo
      ntot = i - 1
      if ( ntot == nvarmax_zt ) then
         l_error = .true.
         write(fstderr,*) "There are more statistical variables listed in ",  &
                          "clubb_vars_zt than allowed for by nvarmax_zt."
         write(fstderr,*) "Check the number of variables listed for clubb_vars_zt ",  &
                          "in the stats namelist, or change nvarmax_zt."
         write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
         error_message = "stats_init_clubb:  number of zt statistical variables exceeds limit"
         clubb_init_errcode = 1
         return
      endif

      stats_zt(j)%num_output_fields = ntot
      stats_zt(j)%kk = nnzp - 1

      allocate( stats_zt(j)%z( stats_zt(j)%kk ), stat=ierr )
      if( ierr /= 0 ) then
      error_message = "stats_init_clubb: Failed to allocate stats_zt%z"
      clubb_init_errcode = 1
      return
    end if

      allocate( stats_zt(j)%accum_field_values( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) then
        error_message = "stats_init_clubb: Failed to allocate stats_zt%accum_field_values"
        clubb_init_errcode = 1
        return
      end if
      allocate( stats_zt(j)%accum_num_samples( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) then
        error_message = "stats_init_clubb: Failed to allocate stats_zt%accum_num_samples"
        clubb_init_errcode = 1
        return
      end if
      allocate( stats_zt(j)%l_in_update( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) then
        error_message = "stats_init_clubb: Failed to allocate stats_zt%l_in_update"
        clubb_init_errcode = 1
        return
      end if
      call stats_zero( stats_zt(j)%kk, stats_zt(j)%num_output_fields, stats_zt(j)%accum_field_values, &
                       stats_zt(j)%accum_num_samples, stats_zt(j)%l_in_update )

      allocate( stats_zt(j)%file%grid_avg_var( stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) then
        error_message = "stats_init_clubb: Failed to allocate stats_zt%file%grid_avg_var"
        clubb_init_errcode = 1
        return
      end if
      allocate( stats_zt(j)%file%z( stats_zt(j)%kk ), stat=ierr )
      if( ierr /= 0 ) then
        error_message = "stats_init_clubb: Failed to allocate stats_zt%file%z"
        clubb_init_errcode = 1
        return
      end if

      !  Default initialization for array indices for zt
      call stats_init_zt_api( hydromet_dim, sclr_dim, edsclr_dim, &
                              hm_metadata%hydromet_list, hm_metadata%l_mix_rat_hm, &
                              clubb_vars_zt, &
                              l_error, &
                              stats_metadata, stats_zt(j) )

      !  Initialize zm (momentum points)

      i = 1
      do while ( ichar(clubb_vars_zm(i)(1:1)) /= 0  .and. &
                 len_trim(clubb_vars_zm(i)) /= 0    .and. &
                 i <= nvarmax_zm )
         i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_zm ) then
         l_error = .true.  ! This will cause the run to stop.
         write(fstderr,*) "There are more statistical variables listed in ",  &
                          "clubb_vars_zm than allowed for by nvarmax_zm."
         write(fstderr,*) "Check the number of variables listed for clubb_vars_zm ",  &
                          "in the stats namelist, or change nvarmax_zm."
         write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
         error_message = "stats_init_clubb:  number of zm statistical variables exceeds limit"
         clubb_init_errcode = 1
         return
      endif

      stats_zm(j)%num_output_fields = ntot
      stats_zm(j)%kk = nnzp

      allocate( stats_zm(j)%z( stats_zm(j)%kk ) )

      allocate( stats_zm(j)%accum_field_values( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ) )
      allocate( stats_zm(j)%accum_num_samples( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ) )
      allocate( stats_zm(j)%l_in_update( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ) )
      call stats_zero( stats_zm(j)%kk, stats_zm(j)%num_output_fields, stats_zm(j)%accum_field_values, &
                       stats_zm(j)%accum_num_samples, stats_zm(j)%l_in_update )

      allocate( stats_zm(j)%file%grid_avg_var( stats_zm(j)%num_output_fields ) )
      allocate( stats_zm(j)%file%z( stats_zm(j)%kk ) )

      call stats_init_zm_api( hydromet_dim, sclr_dim, edsclr_dim, &
                              hm_metadata%hydromet_list, hm_metadata%l_mix_rat_hm, &
                              clubb_vars_zm, &
                              l_error, &
                              stats_metadata, stats_zm(j) )

      !  Initialize rad_zt (radiation points)

      if (stats_metadata%l_output_rad_files) then

         i = 1
         do while ( ichar(clubb_vars_rad_zt(i)(1:1)) /= 0  .and. &
                    len_trim(clubb_vars_rad_zt(i))   /= 0  .and. &
                    i <= nvarmax_rad_zt )
            i = i + 1
         end do
         ntot = i - 1
         if ( ntot == nvarmax_rad_zt ) then
            write(fstderr,*) "There are more statistical variables listed in ",  &
                             "clubb_vars_rad_zt than allowed for by nvarmax_rad_zt."
            write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zt ",  &
                             "in the stats namelist, or change nvarmax_rad_zt."
            write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
            error_message = "stats_init_clubb:  number of rad_zt statistical variables exceeds limit"
            clubb_init_errcode = 1
            return
         endif

        stats_rad_zt(j)%num_output_fields = ntot
        stats_rad_zt(j)%kk = nnrad_zt

        allocate( stats_rad_zt(j)%z( stats_rad_zt(j)%kk ) )

        allocate( stats_rad_zt(j)%accum_field_values( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ) )
        allocate( stats_rad_zt(j)%accum_num_samples( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ) )
        allocate( stats_rad_zt(j)%l_in_update( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ) )

        call stats_zero( stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields, stats_rad_zt(j)%accum_field_values, &
                       stats_rad_zt(j)%accum_num_samples, stats_rad_zt(j)%l_in_update )

        allocate( stats_rad_zt(j)%file%grid_avg_var( stats_rad_zt(j)%num_output_fields ) )
        allocate( stats_rad_zt(j)%file%z( stats_rad_zt(j)%kk ) )

         call stats_init_rad_zt_api( clubb_vars_rad_zt, &
                                     l_error, &
                                     stats_metadata, stats_rad_zt(j) )

         !  Initialize rad_zm (radiation points)

         i = 1
         do while ( ichar(clubb_vars_rad_zm(i)(1:1)) /= 0 .and. &
                    len_trim(clubb_vars_rad_zm(i))   /= 0 .and. &
                    i <= nvarmax_rad_zm )
            i = i + 1
         end do
         ntot = i - 1
         if ( ntot == nvarmax_rad_zm ) then
            l_error = .true.  ! This will cause the run to stop.
            write(fstderr,*) "There are more statistical variables listed in ",  &
                             "clubb_vars_rad_zm than allowed for by nvarmax_rad_zm."
            write(fstderr,*) "Check the number of variables listed for clubb_vars_rad_zm ",  &
                             "in the stats namelist, or change nvarmax_rad_zm."
            write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
            error_message = "stats_init_clubb:  number of rad_zm statistical variables exceeds limit"
            clubb_init_errcode = 1
            return
         endif

         stats_rad_zm(j)%num_output_fields = ntot
         stats_rad_zm(j)%kk = nnrad_zm

         allocate( stats_rad_zm(j)%z( stats_rad_zm(j)%kk ) )

         allocate( stats_rad_zm(j)%accum_field_values( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ) )
         allocate( stats_rad_zm(j)%accum_num_samples( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ) )
         allocate( stats_rad_zm(j)%l_in_update( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ) )

         call stats_zero( stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields, stats_rad_zm(j)%accum_field_values, &
                       stats_rad_zm(j)%accum_num_samples, stats_rad_zm(j)%l_in_update )

         allocate( stats_rad_zm(j)%file%grid_avg_var( stats_rad_zm(j)%num_output_fields ) )
         allocate( stats_rad_zm(j)%file%z( stats_rad_zm(j)%kk ) )

         call stats_init_rad_zm_api( clubb_vars_rad_zm, &
                                     l_error, &
                                     stats_metadata, stats_rad_zm(j) )
      end if ! l_output_rad_files


      !  Initialize sfc (surface point)
      i = 1
      do while ( ichar(clubb_vars_sfc(i)(1:1)) /= 0 .and. &
                 len_trim(clubb_vars_sfc(i))   /= 0 .and. &
                 i <= nvarmax_sfc )
         i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_sfc ) then
         l_error = .true.  ! This will cause the run to stop.
         write(fstderr,*) "There are more statistical variables listed in ",  &
                          "clubb_vars_sfc than allowed for by nvarmax_sfc."
         write(fstderr,*) "Check the number of variables listed for clubb_vars_sfc ",  &
                          "in the stats namelist, or change nvarmax_sfc."
         write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
         error_message = "stats_init_clubb:  number of sfc statistical variables exceeds limit"
         clubb_init_errcode = 1
         return
      endif

      stats_sfc(j)%num_output_fields = ntot
      stats_sfc(j)%kk = 1

      allocate( stats_sfc(j)%z( stats_sfc(j)%kk ) )

      allocate( stats_sfc(j)%accum_field_values( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ) )
      allocate( stats_sfc(j)%accum_num_samples( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ) )
      allocate( stats_sfc(j)%l_in_update( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ) )

      call stats_zero( stats_sfc(j)%kk, stats_sfc(j)%num_output_fields, stats_sfc(j)%accum_field_values, &
                       stats_sfc(j)%accum_num_samples, stats_sfc(j)%l_in_update )

      allocate( stats_sfc(j)%file%grid_avg_var( stats_sfc(j)%num_output_fields ) )
      allocate( stats_sfc(j)%file%z( stats_sfc(j)%kk ) )

      call stats_init_sfc_api( clubb_vars_sfc, &
                               l_error, &
                               stats_metadata, stats_sfc(j) )
    end do

    ! Check for errors

    if ( l_error ) then
       error_message = 'stats_init:  errors found'
       clubb_init_errcode = 1
       return
    endif

    ! Now call add fields

    do i = 1, stats_zt(1)%num_output_fields

      temp1 = trim(stats_zt(1)%file%grid_avg_var(i)%name)
      sub   = temp1
      if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)

        call addfld( trim(sub), (/ 'ilev' /), 'A', &
                     trim(stats_zt(1)%file%grid_avg_var(i)%units), &
                     trim(stats_zt(1)%file%grid_avg_var(i)%description), &
                     sampled_on_subcycle=.true. )
    enddo

    do i = 1, stats_zm(1)%num_output_fields

      temp1 = trim(stats_zm(1)%file%grid_avg_var(i)%name)
      sub   = temp1
      if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)

       call addfld( trim(sub), (/ 'ilev' /), 'A', &
                    trim(stats_zm(1)%file%grid_avg_var(i)%units), &
                    trim(stats_zm(1)%file%grid_avg_var(i)%description), &
                    sampled_on_subcycle=.true. )
    enddo

    if (stats_metadata%l_output_rad_files) then

       do i = 1, stats_rad_zt(1)%num_output_fields
          temp1 = trim(stats_rad_zt(1)%file%grid_avg_var(i)%name)
          sub   = temp1
          if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
          call addfld( trim(sub), (/ 'ilev' /), 'A', &
                       trim(stats_rad_zt(1)%file%grid_avg_var(i)%units), &
                       trim(stats_rad_zt(1)%file%grid_avg_var(i)%description), &
                       sampled_on_subcycle=.true. )
       enddo

       do i = 1, stats_rad_zm(1)%num_output_fields
          temp1 = trim(stats_rad_zm(1)%file%grid_avg_var(i)%name)
          sub   = temp1
          if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
          call addfld( trim(sub), (/ 'ilev' /), 'A', &
                       trim(stats_rad_zm(1)%file%grid_avg_var(i)%units), &
                       trim(stats_rad_zm(1)%file%grid_avg_var(i)%description), &
                       sampled_on_subcycle=.true. )
       enddo
    endif

    do i = 1, stats_sfc(1)%num_output_fields
       temp1 = trim(stats_sfc(1)%file%grid_avg_var(i)%name)
       sub   = temp1
       if (len(temp1) > max_fieldname_len) sub = temp1(1:max_fieldname_len)
       call addfld( trim(sub), horiz_only, 'A', &
                    trim(stats_sfc(1)%file%grid_avg_var(i)%units), &
                    trim(stats_sfc(1)%file%grid_avg_var(i)%description), &
                    sampled_on_subcycle=.true. )
    enddo


    return

  end subroutine stats_init_clubb

#endif


#ifdef CLUBB_SGS

    !-----------------------------------------------------------------------
  subroutine stats_zero( kk, num_output_fields, x, n, l_in_update )

    !     Description:
    !     Initialize stats to zero
    !-----------------------------------------------------------------------

    use clubb_api_module, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd


    implicit none

    !  Input
    integer, intent(in) :: kk, num_output_fields

    !  Output
    real(kind=stat_rknd), intent(out) :: x(:,:,:,:)
    integer(kind=stat_nknd), intent(out) :: n(:,:,:,:)
    logical, intent(out)                 :: l_in_update(:,:,:,:)

    !  Zero out arrays

    if ( num_output_fields > 0 ) then
       x(:,:,:,:) = 0.0_kind_phys
       n(:,:,:,:) = 0
       l_in_update(:,:,:,:) = .false.
    end if

    return

  end subroutine stats_zero

#endif

end module clubb

