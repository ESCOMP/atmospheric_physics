module clubb

  implicit none

  private

  ! Subroutines to make public
  public :: clubb_setup_init, clubb1_run, clubb2_run, clubb3_run

  contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_setup_init(pcols, pver, pverp, pcnst, masterproc, iulog, & ! in
                              sclr_dim, hydromet_dim, nzt_clubb, nzm_clubb, & ! in
                              l_implemented, l_input_fields, clubb_l_do_expldiff_rtm_thlm, & ! in
                              cnst_ndropmixed, subcol_scheme, & ! in
                              clubb_vars_zt, clubb_vars_zm, clubb_vars_sfc, & ! in
                              clubb_vars_rad_zt, clubb_vars_rad_zm, & ! in
                              edsclr_dim, clubb_params_single_col, & ! inout
                              out_zt, out_zm, out_sfc, out_radzt, out_radzm, & ! inout
                              clubb_C2rtthl, clubb_C8, clubb_c11, clubb_c11b, clubb_c14, & ! inout
                              clubb_C_wp3_pr_turb, clubb_c_K10, clubb_mult_coef, & ! inout
                              clubb_Skw_denom_coef, clubb_C2rt, clubb_C2thl, clubb_beta, & ! inout
                              clubb_c6rt, clubb_c6rtb, clubb_c6rtc, clubb_c6thl, clubb_c6thlb, & ! inout
                              clubb_c6thlc, clubb_wpxp_L_thresh, clubb_C7, clubb_C7b, & ! inout
                              clubb_gamma_coef, clubb_c_K10h, clubb_lambda0_stability_coef, & ! inout
                              clubb_lmin_coef, clubb_C8b, clubb_skw_max_mag, clubb_C1, clubb_C1b, & ! inout
                              clubb_gamma_coefb, clubb_up2_sfc_coef, clubb_C4, clubb_C_uu_shr, & ! inout
                              clubb_C_uu_buoy, clubb_c_K1, clubb_c_K2, clubb_nu2, clubb_c_K8, & ! inout
                              clubb_c_K9, clubb_nu9, clubb_C_wp2_splat, clubb_C_invrs_tau_bkgnd, & ! inout
                              clubb_C_invrs_tau_sfc, clubb_C_invrs_tau_shear, clubb_C_invrs_tau_N2, & ! inout
                              clubb_C_invrs_tau_N2_wp2, clubb_C_invrs_tau_N2_xp2, clubb_C_invrs_tau_N2_wpxp, & ! inout
                              clubb_C_invrs_tau_N2_clear_wp3, clubb_bv_efold, clubb_wpxp_Ri_exp, clubb_z_displace, & ! inout
                              lq, stats_zt, stats_zm, stats_sfc, stats_rad_zt, stats_rad_zm, & ! inout
                              stats_metadata, hm_metadata, clubb_config_flags, sclr_idx, & ! inout
                              errmsg, errflg ) ! out

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

    ! Real kinds provided by the CCPP-framework
    use ccpp_kinds, only: kind_phys

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
         time_precision, &
         stats_metadata_type, &
         clubb_config_flags_type, &
         sclr_idx_type, &
         hm_metadata_type, &
         stats, &
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

    !  Input variables, intent(in)
    integer, intent(in) :: pcols, pver, pverp
    integer, intent(in) :: iulog, pcnst
    integer, intent(in) :: sclr_dim, hydromet_dim, nzt_clubb, nzm_clubb
    logical, intent(in) :: masterproc, l_implemented
    logical, intent(in) :: cnst_ndropmixed(:)
    logical, intent(in) :: clubb_l_do_expldiff_rtm_thlm
    logical, intent(in) :: l_input_fields

    character(len=16),         intent(in) :: subcol_scheme
    character(len=var_length), intent(in) :: clubb_vars_zt(:)      ! Variables on the thermodynamic levels
    character(len=var_length), intent(in) :: clubb_vars_zm(:)      ! Variables on the momentum levels
    character(len=var_length), intent(in) :: clubb_vars_rad_zt(:)  ! Variables on the radiation levels
    character(len=var_length), intent(in) :: clubb_vars_rad_zm(:)  ! Variables on the radiation levels
    character(len=var_length), intent(in) :: clubb_vars_sfc(:)     ! Variables at the model surface


    ! Input variables, intent(inout)
    integer, intent(inout) :: edsclr_dim

    real(kind=kind_phys), intent(inout) :: &
      clubb_params_single_col(:,:)    ! Adjustable CLUBB parameters (C1, C2 ...)

    real(kind=kind_phys), allocatable, intent(inout) :: &
      out_zt(:,:,:), out_zm(:,:,:), out_radzt(:,:,:), out_radzm(:,:,:), out_sfc(:,:,:)

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

    logical, intent(inout) :: lq(:)

    ! Variables that contains all the statistics
    type (stats), intent(inout) :: &
      stats_zt(:),      & ! stats_zt grid
      stats_zm(:),      & ! stats_zm grid
      stats_rad_zt(:),  & ! stats_rad_zt grid
      stats_rad_zm(:),  & ! stats_rad_zm grid
      stats_sfc(:)        ! stats_sfc

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type (hm_metadata_type), intent(inout) :: &
      hm_metadata

    type(clubb_config_flags_type), intent(inout) :: &
      clubb_config_flags

    type (sclr_idx_type), intent(inout) :: &
      sclr_idx


    ! Input variables, intent(out)
    integer, intent(out) :: errflg
    character(len=*), intent(out) :: errmsg


    ! Local variables
    real(kind=time_precision) :: dum1, dum2, dum3

    type(err_info_type) :: &
      err_info          ! err_info struct used in CLUBB containing err_code and err_header

    integer :: j, m
    integer :: ierr = 0
    character(len=256) :: alloc_errmsg


    !----- Begin Code -----

    errmsg = ''
    errflg = 0

    if (core_rknd /= kind_phys) then
      errmsg = 'clubb_ini_cam:  CLUBB library core_rknd must match CAM kind_phys and it does not'
      errflg = 1
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
      clubb_config_flags%saturation_formula = saturation_gfdl     ! Goff & Gratch (1946) approx for Sat Vap Pres
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
!$OMP END PARALLEL

    if ( any(err_info%err_code == clubb_fatal_error) ) then
       errmsg = 'clubb_ini_cam: FATAL ERROR CALLING CHECK_CLUBB_SETTINGS_API'
       errflg = 1
       return
    end if

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
       allocate( hm_metadata%hydromet_list(1), stat=ierr, errmsg=alloc_errmsg )
       if( ierr /= 0 ) then
         errmsg = 'clubb_ini_cam: Unable to allocate hm_metadata%hydromet_list: '//trim(alloc_errmsg)
         errflg = 1
         return
       end if
       allocate( hm_metadata%l_mix_rat_hm(1), stat=ierr, errmsg=alloc_errmsg )
       if( ierr /= 0 ) then
         errmsg = 'clubb_ini_cam: Unable to allocate hm_metadata%l_mix_rat_hm: '//trim(alloc_errmsg)
         errflg = 1
         return
       end if
    end if

    !  Initialize statistics, below are dummy variables
    dum1 = 300._kind_phys
    dum2 = 1200._kind_phys
    dum3 = 300._kind_phys

    if (stats_metadata%l_stats) then

      call stats_init_clubb( pcols, .true., dum1, dum2, &
                             nzm_clubb, nzt_clubb, nzm_clubb, dum3, &
                             edsclr_dim, sclr_dim, hydromet_dim, &
                             hm_metadata, stats_metadata, &
                             clubb_vars_zt, clubb_vars_zm, clubb_vars_sfc, &
                             clubb_vars_rad_zt, clubb_vars_rad_zm, &
                             stats_zt(:), stats_zm(:), stats_sfc(:), &
                             stats_rad_zt(:), stats_rad_zm(:), &
                             errmsg, errflg)

      if (errflg /= 0) return

       allocate( out_zt(pcols,pver,stats_zt(1)%num_output_fields), stat=ierr, errmsg=alloc_errmsg )
       if( ierr /= 0 ) then
         errmsg = 'clubb_ini_cam: Unable to allocate out_zt: '//trim(alloc_errmsg)
         errflg = 1
         return
       end if
       allocate( out_zm(pcols,pverp,stats_zm(1)%num_output_fields), stat=ierr, errmsg=alloc_errmsg )
       if( ierr /= 0 ) then
         errmsg = 'clubb_ini_cam: Unable to allocate out_zm: '//trim(alloc_errmsg)
         errflg = 1
         return
       end if
       allocate( out_sfc(pcols,1,stats_sfc(1)%num_output_fields), stat=ierr, errmsg=alloc_errmsg )
       if( ierr /= 0 ) then
         errmsg = 'clubb_ini_cam: Unable to allocate out_sfc: '//trim(alloc_errmsg)
         errflg = 1
         return
       end if

       if ( stats_metadata%l_output_rad_files ) then
          allocate( out_radzt(pcols,pver,stats_rad_zt(1)%num_output_fields), &
                    stat=ierr, errmsg=alloc_errmsg )
          if( ierr /= 0 ) then
            errmsg = 'clubb_ini_cam: Unable to allocate out_radzt: '//trim(alloc_errmsg)
            errflg = 1
            return
          end if
          allocate( out_radzm(pcols,pverp,stats_rad_zm(1)%num_output_fields), &
                    stat=ierr, errmsg=alloc_errmsg )
          if( ierr /= 0 ) then
            errmsg = 'clubb_ini_cam: Unable to allocate out_radzm: '//trim(alloc_errmsg)
            errflg = 1
            return
          end if
       end if

    endif

    end subroutine clubb_setup_init

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

 subroutine clubb1_run( ncol, iam, nstep, lat, lon, hdtime, & ! in
                        pver, pverp, pcnst, clubb_timestep, & ! in
                        nzt_clubb, nzm_clubb, sclr_dim, edsclr_dim, hydromet_dim, & ! in
                        stats_metadata, hm_metadata, clubb_do_adv, first_step, first_restart_step, & ! in
                        single_column, scm_cambfb_mode, scm_clubb_iop_name, & ! in
                        shr_const_karman, shr_const_pi, shr_const_g, omega_const, theta0, & ! in
                        macmic_it, top_lev, rtpthlp_const, wpthlp_const, wprtp_const, sclr_tol, & ! in
                        ts_nudge, rtm_min, rtm_nudge_max_altitude, & ! in
                        wp3_const, cld_macmic_num_steps, clubb_params_single_col, & ! in
                        cpair, cpairv, rair, inv_p0_clubb, rairv, zvir, latvap, latice, & ! in
                        h2otrip, rga, gravit, clubb_rnevap_effic, do_cldcool, do_rainturb, & ! in
                        l_implemented, grid_type, lq, deep_scheme, & ! in
                        state_q, state_t, state_pmid, state_zm, & ! in
                        state_phis, state_pdel, state_pdeldry, state_exner, & ! in
                        state_pint, state_zi, state_omega, & ! in
                        wsx, wsy, cflx, shf, landfrac, ts, tpert, & ! in
                        sclr_idx, clubb_l_ascending_grid, clubb_do_energyfix, & ! in
                        ixq, ixcldliq, ixcldice, ixrtpthlp, ixwpthlp, & ! in
                        ixwprtp, ixwp3, ixwp2, ixthlp2, ixrtp2, ixup2, ixvp2, & ! in
                        clubb_l_intr_sfc_flux_smooth, clubb_config_flags, & ! in
                        apply_const, gr, ztodtptr, u, v, s, wprcp, & ! inout
                        ptend_q, ptend_u, ptend_v, ptend_s, & ! inout
                        pdf_params_chnk, pdf_params_zm_chnk, pdf_implicit_coefs_terms_chnk, & ! inout
                        eleak, se_dis, rho_zm, rho_zt, exner, cloud_frac, & ! inout
                        zi_g, zt_g, grid_dx, grid_dy, & ! inout
                        mf_dry_a, mf_moist_a, mf_dry_w, mf_moist_w, & ! inout
                        mf_dry_qt, mf_moist_qt, mf_dry_thl, mf_moist_thl, & ! inout
                        mf_dry_u, mf_moist_u, mf_dry_v, mf_moist_v, mf_moist_qc, & ! inout
                        s_ae, s_aw, s_awthl, s_awqt, s_awu, s_awv, & ! inout
                        mf_thlflx, mf_qtflx, & ! inout
                        ztopm1_macmic, ddcp, ddcp_macmic, cbm1, cbm1_macmic, & ! inout
                        thlm, rtm, um, vm, wm_zt, rcm, rcm_in_layer, & ! inout
                        wp2_pbuf, wp3_pbuf, wpthlp_pbuf, wprtp_pbuf, & ! inout
                        rtpthlp_pbuf, rtp2_pbuf, thlp2_pbuf, rtp3_pbuf, & ! inout
                        thlp3_pbuf, up2_pbuf, vp2_pbuf, up3_pbuf, vp3_pbuf, & ! inout
                        upwp_pbuf, vpwp_pbuf, wpthvp_pbuf, wp2thvp_pbuf, wp2up_pbuf, & ! inout
                        rtpthvp_pbuf, thlpthvp_pbuf, pdf_zm_w_1_pbuf, pdf_zm_w_2_pbuf, & ! inout
                        pdf_zm_varnce_w_1_pbuf, pdf_zm_varnce_w_2_pbuf, pdf_zm_mixt_frac_pbuf, & ! inout
                        wp2rtp_pbuf, wp2thlp_pbuf, uprcp_pbuf, vprcp_pbuf, rc_coef_zm_pbuf, & ! inout
                        wp4_pbuf, wpup2_pbuf, wpvp2_pbuf, wp2up2_pbuf, wp2vp2_pbuf, cld_pbuf, & ! inout
                        khzm_pbuf, pblh_pbuf, tke_pbuf, & ! inout
                        ice_supersat_frac_pbuf, relvar_pbuf, prer_evap_pbuf, qrl_pbuf, & ! inout
                        rtp2_mc_zt_pbuf, thlp2_mc_zt_pbuf, wprtp_mc_zt_pbuf, & ! inout
                        wpthlp_mc_zt_pbuf, rtpthlp_mc_zt_pbuf, ttend_clubb_pbuf, & ! inout
                        upwp_clubb_gw_pbuf, vpwp_clubb_gw_pbuf, thlp2_clubb_gw_pbuf, & ! inout
                        wpthlp_clubb_gw_pbuf, ttend_clubb_mc_pbuf, upwp_clubb_gw_mc_pbuf, & ! inout
                        vpwp_clubb_gw_mc_pbuf, thlp2_clubb_gw_mc_pbuf, wpthlp_clubb_gw_mc_pbuf, & ! inout
                        stats_zt, stats_zm, stats_sfc, stats_rad_zt, stats_rad_zm, & ! inout
                        out_zt, out_zm, out_sfc, out_radzt, out_radzm, & ! inout
                        invrs_cpairv, clubbtop_pbuf, & ! inout
                        mf_ztop_nadv, mf_L0_nadv, mf_cape_nadv, max_cfl_nadv, mf_freq_nadv, & ! inout
                        s_aww, s_awthlup, s_awqtup, s_awuup, s_awvup, s_awthldn, & ! inout
                        s_awqtdn, s_awudn, s_awvdn, mf_precc, mf_thlflxup, mf_qtflxup, & ! inout
                        mf_uflxup, mf_vflxup, mf_thlflxdn, mf_qtflxdn, mf_uflxdn, mf_vflxdn, & ! inout
                        mf_uflx, mf_vflx, mf_thvflx, mf_rcm, mf_ent_nadv, mf_cloudfrac_zt, mf_qc_zt, & ! inout
                        mf_upa, mf_upw, mf_dnw, mf_upmf, mf_upqt, mf_dnqt, mf_upthl, mf_dnthl, & ! inout
                        mf_upthv, mf_upth, mf_upqc, mf_upbuoy, mf_updet, mf_upent, & ! inout
                        mf_thlforcup_nadv, mf_qtforcup_nadv, mf_thlforcdn_nadv, & ! inout
                        mf_qtforcdn_nadv, mf_thlforc_nadv, mf_qtforc_nadv, & ! inout
                        ztopma, ztopmn, mf_sqtup, mf_sqtdn, prec_dp_pbuf, snow_dp_pbuf, ddcpmn, & ! inout
                        s_mfup_macmic, s_entup_macmic, s_detup_macmic, s_mfdn_macmic, s_entdn_macmic, & ! inout
                        s_detdn_macmic, s_awwup_macmic, s_awup_macmic, mf_sqtac_macmic, mf_sqtev_macmic, mf_kctop_macmic, & ! inout
                        rprddp, evapcdp, mu, md, du, eu, ed, dpg, dsubcld, jt, maxg, ideep, wup_gath, & ! inout
                        errmsg, errflg ) ! out

    use ccpp_kinds, only: kind_phys
    use clubb_mf,              only: integrate_mf, do_clubb_mf, clubb_mf_nup, &
                                     clubb_mf_Lopt, clubb_mf_ddalph, clubb_mf_up_ndt, clubb_mf_cp_ndt, &
                                     do_clubb_mf_cmt, do_clubb_mf_addtke
    use atmos_phys_pbl_utils,  only: calc_friction_velocity, calc_ideal_gas_rrho
    use wv_saturation,         only: qsat
    use clubb_api_module,      only: &
      nparams, &
      calc_derrived_params_api, &
      check_parameters_api, &
      advance_clubb_core_api, &
      zt2zm_api, &
      w_tol_sqd, &
      rt_tol, &
      thl_tol, &
      stats_begin_timestep_api, &
      cleanup_err_info_api, &
      calculate_thlp2_rad_api, update_xp2_mc_api, &
      fstderr, &
      ipdf_post_advance_fields, &
      grid, &
      stats, &
      setup_grid_api, &
      ic_K, &
      cleanup_grid_api, &
      nu_vertical_res_dep

    ! Import setup for CLUBB error messaging
    use clubb_api_module, only: &
      clubb_fatal_error,    & ! Error code value to indicate a fatal error
      stats_metadata_type,  &
      hm_metadata_type, &
      clubb_config_flags_type, &
      sclr_idx_type,        &
      err_info_type,        &
      pdf_parameter,        &
      init_err_info_api,    &
      implicit_coefs_terms, &
      zm2zt_api


    ! Input variables, intent(in)
    integer, intent(in) :: ncol, iam, ixq, ixcldliq, ixcldice, ixrtpthlp, ixwpthlp, &
                           ixwprtp, ixwp3, ixwp2, ixthlp2, ixrtp2, ixup2, ixvp2, sclr_dim, edsclr_dim, &
                           macmic_it, top_lev, cld_macmic_num_steps, grid_type, hydromet_dim, nstep, &
                           nzt_clubb, nzm_clubb, pver, pverp, pcnst
    real(kind_phys), intent(in) :: hdtime, clubb_timestep
    real(kind_phys), intent(in) :: sclr_tol(:)
    real(kind_phys), intent(in) :: shr_const_karman, shr_const_pi, shr_const_g, omega_const, cpair, &
                                   h2otrip, rga, inv_p0_clubb, rair, zvir, latvap, latice, gravit, theta0, &
                                   ts_nudge, rtm_min, rtm_nudge_max_altitude, clubb_rnevap_effic
    real(kind_phys), intent(in) :: cpairv(:,:), rairv(:,:)
    real(kind_phys), intent(in) :: rtpthlp_const, wpthlp_const, wprtp_const, wp3_const
    real(kind_phys), intent(in) :: state_q(:,:,:)
    real(kind_phys), intent(in) :: wsx(:), wsy(:), shf(:)
    real(kind_phys), intent(in) :: cflx(:,:)
    real(kind_phys), intent(in) :: clubb_params_single_col(:,:)
    real(kind_phys), intent(in) :: lat(:), lon(:), state_phis(:)
    real(kind_phys), intent(in) :: state_pint(:,:)
    real(kind_phys), intent(in) :: state_pmid(:,:)
    real(kind_phys), intent(in) :: landfrac(:)
    real(kind_phys), intent(in) :: ts(:)
    real(kind_phys), intent(in) :: state_pdel(:,:), state_pdeldry(:,:), state_omega(:,:), &
                                   state_t(:,:), state_zm(:,:), state_zi(:,:), state_exner(:,:)

    logical, intent(in) :: clubb_do_adv, first_step, first_restart_step, l_implemented, &
                           !do_clubb_mf,
                           single_column, scm_cambfb_mode, clubb_l_intr_sfc_flux_smooth, &
                           clubb_l_ascending_grid, clubb_do_energyfix, do_cldcool, do_rainturb

    logical, intent(in) :: lq(:)

    character(len=20), intent(in) :: scm_clubb_iop_name
    character(len=16), intent(in) :: deep_scheme

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type(clubb_config_flags_type), intent(in) :: &
      clubb_config_flags


    ! Input variables, intent(inout)
    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type(pdf_parameter), intent(inout) :: &
      pdf_params_chnk, &                ! PDF parameters (thermo. levs.) [units vary]
      pdf_params_zm_chnk

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms_chnk

    type (hm_metadata_type), intent(inout) :: &
      hm_metadata

    real(kind_phys), intent(inout) :: invrs_cpairv(:,:)
    real(kind_phys), intent(inout) :: zi_g(:,:), wprcp(:,:), rho_zm(:,:), rho_zt(:,:)
    real(kind_phys), intent(inout) :: ptend_u(:,:), ptend_v(:,:), ptend_s(:,:)
    real(kind_phys), intent(inout) :: eleak(:), se_dis(:)
    real(kind_phys), intent(inout) :: s(:,:), u(:,:), v(:,:)
    real(kind_phys), intent(inout) :: ptend_q(:,:,:)

    real(kind_phys), intent(inout) :: mf_ztop_nadv(:,:), mf_L0_nadv(:,:), mf_cape_nadv(:,:)
    real(kind_phys), intent(inout) :: max_cfl_nadv(:), mf_freq_nadv(:), prec_dp_pbuf(:), snow_dp_pbuf(:)
    real(kind_phys), intent(inout) :: s_aww(:,:), s_awthlup(:,:), s_awqtup(:,:), s_awuup(:,:), s_awvup(:,:), s_awthldn(:,:)
    real(kind_phys), intent(inout) :: s_awqtdn(:,:), s_awudn(:,:), s_awvdn(:,:), mf_precc(:,:), mf_thlflxup(:,:)
    real(kind_phys), intent(inout) :: mf_qtflxup(:,:), mf_uflxup(:,:), mf_vflxup(:,:), mf_thlflxdn(:,:), mf_qtflxdn(:,:)
    real(kind_phys), intent(inout) :: mf_uflxdn(:,:), mf_vflxdn(:,:), mf_uflx(:,:), mf_vflx(:,:), mf_sqtup(:,:), mf_sqtdn(:,:)
    real(kind_phys), intent(inout) :: mf_thvflx(:,:), mf_rcm(:,:), mf_ent_nadv(:,:), mf_cloudfrac_zt(:,:), mf_qc_zt(:,:)
    real(kind_phys), intent(inout) :: mf_upa(:,:,:), mf_upw(:,:,:), mf_dnw(:,:,:), ddcpmn(:,:,:)
    real(kind_phys), intent(inout) :: mf_upmf(:,:,:), mf_upqt(:,:,:), mf_dnqt(:,:,:), mf_upthl(:,:,:), mf_dnthl(:,:,:), mf_upthv(:,:,:)
    real(kind_phys), intent(inout) :: mf_upth(:,:,:), mf_upqc(:,:,:), mf_upbuoy(:,:,:), mf_updet(:,:,:), mf_upent(:,:,:), ztopmn(:,:,:)
    real(kind_phys), intent(inout) :: ztopm1_macmic(:,:), ddcp(:,:), ddcp_macmic(:,:)
    real(kind_phys), intent(inout) :: cbm1_macmic(:), cbm1(:), tpert(:)
    real(kind_phys), intent(inout) :: mf_thlforcup_nadv(:,:), mf_qtforcup_nadv(:,:), mf_thlforcdn_nadv(:,:)
    real(kind_phys), intent(inout) :: mf_qtforcdn_nadv(:,:), mf_thlforc_nadv(:,:), mf_qtforc_nadv(:,:), ztopma(:,:)
    real(kind_phys), intent(inout) :: s_mfup_macmic(:,:), s_entup_macmic(:,:), s_detup_macmic(:,:), s_mfdn_macmic(:,:)
    real(kind_phys), intent(inout) :: s_entdn_macmic(:,:), s_detdn_macmic(:,:), s_awwup_macmic(:,:), s_awup_macmic(:,:)
    real(kind_phys), intent(inout) :: mf_sqtac_macmic(:,:), mf_sqtev_macmic(:,:)
    real(kind_phys), intent(inout) :: mf_kctop_macmic(:)

    ! deep-hookup pbuf pointers
    real(kind_phys), intent(inout) :: rprddp(:,:)    ! Deep conv precip production (kg/kg/s - grid avg)
    real(kind_phys), intent(inout) :: evapcdp(:,:)   ! Deep conv precip evaporation (kg/kg/s - grid avg)
    real(kind_phys), intent(inout) :: mu(:,:)        ! Updraft mass flux (positive) (pcols,pver)
    real(kind_phys), intent(inout) :: md(:,:)        ! Downdraft mass flux (negative) (pcols,pver)
    real(kind_phys), intent(inout) :: du(:,:)        ! Mass detrain rate from updraft (pcols,pver)
    real(kind_phys), intent(inout) :: eu(:,:)        ! Mass entrain rate into updraft (pcols,pver)
    real(kind_phys), intent(inout) :: ed(:,:)        ! Mass entrain rate into downdraft (pcols,pver)
    real(kind_phys), intent(inout) :: dpg(:,:)       ! Delta pressure between interfaces, mb (pcols,pver)
    real(kind_phys), intent(inout) :: dsubcld(:)     ! Delta pressure from cloud base to sfc (pcols)
    integer, intent(inout) :: jt(:)          ! Index of cloud top for each column (pcols)
    integer, intent(inout) :: maxg(:)        ! Index of cloud base for each column (pcols)
    integer, intent(inout) :: ideep(:)       ! Gathering array (pcols)
    real(kind_phys), intent(inout) :: wup_gath(:,:)  ! Mass-flux-weighted ensemble updraft speed, gathered (pcols,pver)

    ! Variables that contains all the statistics
    type (stats), intent(inout) :: &
      stats_zt(:),      & ! stats_zt grid
      stats_zm(:),      & ! stats_zm grid
      stats_rad_zt(:),  & ! stats_rad_zt grid
      stats_rad_zm(:),  & ! stats_rad_zm grid
      stats_sfc(:)        ! stats_sfc

    real(kind=kind_phys), allocatable, intent(inout) :: &
      out_zt(:,:,:), out_zm(:,:,:), out_radzt(:,:,:), out_radzm(:,:,:), out_sfc(:,:,:)

    real(kind_phys), intent(inout) :: wp2_pbuf(:,:)
    real(kind_phys), intent(inout) :: wp3_pbuf(:,:)
    real(kind_phys), intent(inout) :: wpthlp_pbuf(:,:)
    real(kind_phys), intent(inout) :: wprtp_pbuf(:,:)
    real(kind_phys), intent(inout) :: rtpthlp_pbuf(:,:)
    real(kind_phys), intent(inout) :: rtp2_pbuf(:,:)
    real(kind_phys), intent(inout) :: thlp2_pbuf(:,:)
    real(kind_phys), intent(inout) :: rtp3_pbuf(:,:)
    real(kind_phys), intent(inout) :: thlp3_pbuf(:,:)
    real(kind_phys), intent(inout) :: up2_pbuf(:,:)
    real(kind_phys), intent(inout) :: vp2_pbuf(:,:)
    real(kind_phys), intent(inout) :: up3_pbuf(:,:)
    real(kind_phys), intent(inout) :: vp3_pbuf(:,:)
    real(kind_phys), intent(inout) :: upwp_pbuf(:,:)
    real(kind_phys), intent(inout) :: vpwp_pbuf(:,:)
    real(kind_phys), intent(inout) :: wpthvp_pbuf(:,:)
    real(kind_phys), intent(inout) :: wp2thvp_pbuf(:,:)
    real(kind_phys), intent(inout) :: wp2up_pbuf(:,:)
    real(kind_phys), intent(inout) :: rtpthvp_pbuf(:,:)
    real(kind_phys), intent(inout) :: thlpthvp_pbuf(:,:)
    real(kind_phys), intent(inout) :: pdf_zm_w_1_pbuf(:,:)
    real(kind_phys), intent(inout) :: pdf_zm_w_2_pbuf(:,:)
    real(kind_phys), intent(inout) :: pdf_zm_varnce_w_1_pbuf(:,:)
    real(kind_phys), intent(inout) :: pdf_zm_varnce_w_2_pbuf(:,:)
    real(kind_phys), intent(inout) :: pdf_zm_mixt_frac_pbuf(:,:)
    real(kind_phys), intent(inout) :: wp2rtp_pbuf(:,:)
    real(kind_phys), intent(inout) :: wp2thlp_pbuf(:,:)
    real(kind_phys), intent(inout) :: uprcp_pbuf(:,:)
    real(kind_phys), intent(inout) :: vprcp_pbuf(:,:)
    real(kind_phys), intent(inout) :: rc_coef_zm_pbuf(:,:)
    real(kind_phys), intent(inout) :: wp4_pbuf(:,:)
    real(kind_phys), intent(inout) :: wpup2_pbuf(:,:)
    real(kind_phys), intent(inout) :: wpvp2_pbuf(:,:)
    real(kind_phys), intent(inout) :: wp2up2_pbuf(:,:)
    real(kind_phys), intent(inout) :: wp2vp2_pbuf(:,:)
    real(kind_phys), intent(inout) :: cld_pbuf(:,:)
    real(kind_phys), intent(inout) :: khzm_pbuf(:,:)
    real(kind_phys), intent(inout) :: pblh_pbuf(:)
    real(kind_phys), intent(inout) :: tke_pbuf(:,:)
    real(kind_phys), intent(inout) :: ice_supersat_frac_pbuf(:,:)
    real(kind_phys), intent(inout) :: relvar_pbuf(:,:)

    real(kind_phys), intent(inout) :: prer_evap_pbuf(:,:)
    real(kind_phys), intent(inout) :: qrl_pbuf(:,:)

    real(kind_phys), intent(inout) :: rtp2_mc_zt_pbuf(:,:)
    real(kind_phys), intent(inout) :: thlp2_mc_zt_pbuf(:,:)
    real(kind_phys), intent(inout) :: wprtp_mc_zt_pbuf(:,:)
    real(kind_phys), intent(inout) :: wpthlp_mc_zt_pbuf(:,:)
    real(kind_phys), intent(inout) :: rtpthlp_mc_zt_pbuf(:,:)

    real(kind_phys), intent(inout) :: ttend_clubb_pbuf(:,:)
    real(kind_phys), intent(inout) :: upwp_clubb_gw_pbuf(:,:)
    real(kind_phys), intent(inout) :: vpwp_clubb_gw_pbuf(:,:)
    real(kind_phys), intent(inout) :: thlp2_clubb_gw_pbuf(:,:)
    real(kind_phys), intent(inout) :: wpthlp_clubb_gw_pbuf(:,:)

    real(kind_phys), intent(inout) :: ttend_clubb_mc_pbuf(:,:)
    real(kind_phys), intent(inout) :: upwp_clubb_gw_mc_pbuf(:,:)
    real(kind_phys), intent(inout) :: vpwp_clubb_gw_mc_pbuf(:,:)
    real(kind_phys), intent(inout) :: thlp2_clubb_gw_mc_pbuf(:,:)
    real(kind_phys), intent(inout) :: wpthlp_clubb_gw_mc_pbuf(:,:)

    real(kind_phys), intent(inout) :: &
      grid_dx(:), grid_dy(:)                    ! CAM grid [m]

    real(kind_phys), intent(inout) :: &
      rtm(:,:),                            & ! mean moisture mixing ratio                                  [kg/kg]
      thlm(:,:),                           & ! mean temperature                                                    [K]
      rcm(:,:),                            & ! CLUBB cloud water mixing ratio                [kg/kg]
      um(:,:),                             & ! mean east-west wind                                               [m/s]
      vm(:,:),                             & ! mean north-south wind                                     [m/s]
      rcm_in_layer(:,:), &
      wm_zt(:,:), &
      exner(:,:), &
      cloud_frac(:,:), &
      zt_g(:,:)

    real(kind_phys), intent(inout) :: &  ! MF plume
      mf_dry_a(:,:),   mf_moist_a(:,:),    &
      mf_dry_w(:,:),   mf_moist_w(:,:),    &
      mf_dry_qt(:,:),  mf_moist_qt(:,:),   &
      mf_dry_thl(:,:), mf_moist_thl(:,:),  &
      mf_dry_u(:,:),   mf_moist_u(:,:),    &
      mf_dry_v(:,:),   mf_moist_v(:,:),    &
                  mf_moist_qc(:,:),   &
      s_ae(:,:),       s_aw(:,:),          &
      s_awthl(:,:),    s_awqt(:,:),        &
      s_awu(:,:),      s_awv(:,:),         &
      mf_thlflx(:,:),  mf_qtflx(:,:)

    real(kind_phys), intent(inout) :: apply_const

    integer, intent(inout) :: clubbtop_pbuf(:)

    real(kind_phys), intent(inout) :: ztodtptr

    type(grid), intent(inout) :: &
      gr

    ! Input variables, intent(out)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! ------------- Local variables --------------
    real(kind_phys) :: rad2deg

    real(kind_phys) :: &
      ustar2(ncol), &
      rrho(ncol), &
      deltaz(ncol), &
      fcor(ncol), &                             ! Coriolis forcing                                          [s^-1]
      fcor_y(ncol), &                           ! Non-traditional coriolis forcing                          [s^-1]
      sfc_elevation(ncol), &                              ! Elevation of ground                                           [m AMSL][m]
      wpthlp_sfc(ncol), &                       ! w' theta_l' at surface                      [(m K)/s]
      wprtp_sfc(ncol), &                        ! w' r_t' at surface                          [(kg m)/( kg s)]
      upwp_sfc(ncol), &                         ! u'w' at surface                             [m^2/s^2]
      vpwp_sfc(ncol), &                         ! v'w' at surface                             [m^2/s^2]
      p_sfc(ncol), &
      upwp_sfc_pert(ncol), &                    ! perturbed u'w' at surface                   [m^2/s^2]
      vpwp_sfc_pert(ncol)                    ! perturbed v'w' at surface                   [m^2/s^2]

    type(err_info_type) :: &
      err_info

    type(nu_vertical_res_dep) :: &
      nu_vert_res_dep

    real(kind_phys), dimension(ncol,sclr_dim) :: &
      wpsclrp_sfc            ! Scalar flux at surface                        [{units vary} m/s]

    real(kind_phys), dimension(ncol,edsclr_dim) :: &
      wpedsclrp_sfc        ! Eddy-scalar flux at surface                   [{units vary} m/s]

    real(kind_phys), dimension(ncol,nzt_clubb) :: &
      thlm_forcing,                   & ! theta_l forcing (thermodynamic levels)        [K/s]
      rtm_forcing,                    & ! r_t forcing (thermodynamic levels)            [(kg/kg)/s]
      um_forcing,                     & ! u wind forcing (thermodynamic levels)           [m/s/s]
      vm_forcing,                     & ! v wind forcing (thermodynamic levels)           [m/s/s]
      rtm_ref,                        & ! Initial profile of rtm                        [kg/kg]
      thlm_ref,                       & ! Initial profile of thlm                       [K]
      um_ref,                         & ! Initial profile of um                         [m/s]
      vm_ref,                         & ! Initial profile of vm                         [m/s]
      ug,                             & ! U geostrophic wind                            [m/s]
      vg,                             & ! V geostrophic wind                            [m/s]
      p_in_Pa,                        & ! Air pressure (thermodynamic levels)             [Pa]
      rho_ds_zt,                      & ! Dry, static density on thermodynamic levels   [kg/m^3]
      invrs_rho_ds_zt,                & ! Inv. dry, static density on thermo. levels    [m^3/kg]
      thv_ds_zt,                      & ! Dry, base-state theta_v on thermo. levels     [K]
      rfrzm,                          &
      rvm,                            & ! water vapor mixing ratio                      [kg/kg]
      um_pert,                        & ! Perturbed U wind                              [m/s]
      vm_pert,                        & ! Perturbed V wind                              [m/s]
      khzt,                           & ! eddy diffusivity on thermo grids              [m^2/s]
      w_up_in_cloud,                  &
      w_down_in_cloud,                &
      cloudy_updraft_frac,            &
      cloudy_downdraft_frac,          &
      cloud_cover,                    & ! CLUBB output of in-cloud cloud fraction       [fraction]
      pre,                            & ! input for precip evaporation
      qrl_clubb,                      &
      qclvar,                         & ! cloud water variance                          [kg^2/kg^2]
      Lscale,                         &
      dz_g,                           & ! thickness of layer                            [m]
      invrs_dz_g,                     & ! Inverse of layer thickness                    [1/m]
      invrs_exner_zt,                 & ! thermodynamic grid

      !MF Plume Variables on the thermodynamic grid.
      kappa_zt,                       & ! thermodynamic grid
      qc_zt,                          & ! thermodynamic grid
      th_zt,                          & ! thermodynamic grid
      qv_zt,                          & ! thermodynamic grid
      mf_thlforcup, mf_qtforcup,      & ! variables that represents a forcing or a tendency ($d/dt$)
      mf_thlforcdn, mf_qtforcdn,      & ! applied to a mean state variable (like $\theta_l$ or $q_t$)
      mf_sthlup,                      &
      mf_sthldn,                      &
      mf_sthl,                        &
      mf_sqt,                         &
      mf_sqtac, mf_sqtev,             &
      mf_sqtac_nadv, mf_sqtev_nadv

    real(kind_phys), dimension(ncol,nzm_clubb) :: &
      thlp2_rad,                &
      wprtp_forcing,            &
      wpthlp_forcing,           &
      rtp2_forcing,             &
      thlp2_forcing,            &
      rtpthlp_forcing,          &
      wm_zm,                    & ! w mean wind component on momentum levels              [m/s]
      rho_ds_zm,                & ! Dry, static density on momentum levels                      [kg/m^3]
      invrs_rho_ds_zm,          & ! Inv. dry, static density on momentum levels                 [m^3/kg]
      thv_ds_zm,                & ! Dry, base-state theta_v on momentum levels                  [K]
      upwp_pert,                & ! Perturbed u'w'                                        [m^2/s^2]
      vpwp_pert,                & ! Perturbed v'w'                                        [m^2/s^2]
      khzm,                     & ! Eddy diffusivity of heat/moisture on momentum levels  [m^2/s]
      thlprcp,                  &
      invrs_tau_zm,             & ! CLUBB output of 1 divided by time-scale               [1/s]
      rtp2_mc,                  & ! total water tendency from rain evap
      thlp2_mc,                 & ! thetal tendency from rain evap
      wprtp_mc,                 &
      wpthlp_mc,                &
      rtpthlp_mc,               &
      Lscale_zm,                &

      ! MF Plume variables on momentum levels.
      s_ac, s_aup, s_adn, s_awup, s_awdn, s_awwup, s_awwdn, &
      mf_qc, mf_qc_nadv, mf_cloudfrac, mf_cloudfrac_nadv, mf_rcm_nadv, &
      s_mfup,     s_entup,   s_detup,  &
      s_mfdn,     s_entdn,   s_detdn,  &
      s_mfup_nadv, s_entup_nadv, s_detup_nadv, &
      s_mfdn_nadv, s_entdn_nadv, s_detdn_nadv, &
      s_awwup_nadv, s_awup_nadv

    real(kind_phys), dimension(ncol,nzm_clubb) :: &
      ! MF local momentum vars
      rtm_zm,     thlm_zm,       & ! momentum grid
      th_zm,      qv_zm,         & ! momentum grid
      qc_zm,      tke_zm,        & ! momentum grid
      kappa_zm,   p_in_Pa_zm,    & ! momentum grid
      invrs_exner_zm               ! momentum grid

    real(kind_phys), dimension(ncol,nzm_clubb,clubb_mf_nup) ::     &
      mf_dna, mf_dnthv, mf_dnth, mf_dnqc

    real(kind_phys), dimension(ncol,nzt_clubb,sclr_dim) :: &
      sclrm_forcing,  & ! Passive scalar forcing                        [{units vary}/s]
      sclrm,          & ! Passive scalar mean (thermo. levels)          [units vary]
      sclrp3            ! sclr'^3 (thermo. levels)                      [{units vary}^3]

    real(kind_phys), dimension(ncol,nzm_clubb,sclr_dim) :: &
      sclrp2,         & ! sclr'^2 (momentum levels)                     [{units vary}^2]
      sclrprtp,       & ! sclr'rt' (momentum levels)                    [{units vary} (kg/kg)]
      sclrpthlp,      & ! sclr'thlp' (momentum levels)                  [{units vary} (K)]
      wpsclrp,        & ! w'sclr' (momentum levels)                     [{units vary} m/s]
      sclrpthvp         ! sclr'th_v' (momentum levels)                  [{units vary} (K)]

    real(kind_phys), dimension(ncol,nzt_clubb,edsclr_dim) :: &
      edsclrm_forcing,  & ! Eddy passive scalar forcing                 [{units vary}/s]
      edsclr                 ! Scalars to be diffused through CLUBB      [units vary]

    real(kind_phys), dimension(ncol,nzt_clubb,hydromet_dim) :: &
      wp2hmp,       &
      rtphmp_zt,    &
      thlphmp_zt

    real(kind_phys), dimension(ncol,nzm_clubb,hydromet_dim) :: &
      wphydrometp

    real(kind_phys) :: clubb_s(ncol,nzt_clubb)         ! diagnosed dry static energy from clubb

    real(kind_phys) :: &
      dum1,                     & ! dummy variable                                [units vary]
      invrs_hdtime,             &
      invrs_macmic_num_steps,   &
      lmin,                     &
      mixt_frac_max_mag,        &
      dtime,                    & ! CLUBB time step                               [s]
      ubar,                     & ! surface wind                                  [m/s]
      ustar,                    & ! surface stress                                                      [m/s]
      bflx22,                   & ! Variable for buoyancy flux for pbl            [K m/s]
      zo,                       & ! roughness height                              [m]
      relvarmax,                &
      ! Variables below are needed to compute energy integrals for conservation
      te_a, se_a, ke_a, wv_a, wl_a, &
      te_b, se_b, ke_b, wv_b, wl_b

    real(kind_phys), dimension(ncol,clubb_mf_nup) :: mf_ztop,    &
                                                   mf_ztopm1,  mf_ztopm1_nadv, &
                                                   mf_L0,       &
                                                   mf_ddcp,    mf_ddcp_nadv,   &
                                                   mf_cape

    real(kind_phys), dimension(ncol)       :: mf_precc_nadv, mf_snow_nadv,&
                                              mf_cbm1,       mf_cbm1_nadv

    real(kind_phys), dimension(ncol) :: mf_kctop, mf_kctop_nadv
    integer :: lengath, tmptop, kcam, icol
    real(kind_phys) :: tmpcfl, tmpnet
    real(kind_phys), dimension(ncol,pver) :: thlcam
    real(kind_phys), dimension(ncol) :: thl700, thl1000, mf_lts

    intrinsic :: max

    real(kind_phys), dimension(ncol,nparams) :: &
      clubb_params

    integer :: &
      i, k, tt, ixind, nadv, n,      & ! Loop variables
      k_cam, k_clubb, sclr, iedsclr, & ! Loop variables
      icnt, &
      stats_nsamp, stats_nout         ! Stats sampling and output intervals for CLUBB [timestep]

    ! -------------------------------------------------
    ! ----------------- Start -------------------------

    errmsg = ''
    errflg = 0

    rad2deg = 180.0_kind_phys / shr_const_pi

    ! Initialize err_info with parallelization and geographical info
    ! BAS putting a -1 here in place of "lchnk"
    call init_err_info_api(ncol, -1, iam, lat*rad2deg, lon*rad2deg, err_info)

    ! the rhinv preprocessing (RH at 500 hPa / column RH) has moved into
    ! integrate_mf, which computes it locally from its own column inputs
    ! when do_clubb_mf_rhtke or clubb_mf_Lopt 7/8 need it.
    !
    ! lower-tropospheric stability thl700 - thl1000 (K) for the marine-Sc
    ! plume inhibition in integrate_mf (do_clubb_mf_invswitch).  Mirrors the
    ! CLUBB-core expldiff criterion (advance_clubb_core_module, pvertinterp
    ! of thlm to 700 and 1000 hPa), including pvertinterp's boundary clamping.
    ! NOTE deliberately kept verbatim from the first implementation so runs are
    ! bit-for-bit reproducible against it: state%exner is SURFACE-referenced
    ! ((pint(pverp)/pmid)**cappav, dp_coupling), not (p0/pmid)**cappa, so this
    ! thl is scaled by ~(psfc/p0)**cappa relative to the standard definition
    ! (~+0.4% over ocean).  The p0-based correction is a known follow-up: swap
    ! state_loc%exner for (1.e5_kind_phys/state_loc%pmid)**cappa when a deliberate
    ! answer change is acceptable.
    if (do_clubb_mf) then
       ! Saturation Vapor Pressure calculation
       do k = 1, pver
          thlcam(:ncol,k) = ( state_t(:ncol,k) &
                              - (latvap/cpair)*state_q(:ncol,k,ixcldliq) ) &
                            * state_exner(:ncol,k)
       end do
       call mf_pinterp(ncol, pver, state_pmid(:ncol,:),  70000._kind_phys, thlcam(:ncol,:), thl700)
       call mf_pinterp(ncol, pver, state_pmid(:ncol,:), 100000._kind_phys, thlcam(:ncol,:), thl1000)
       mf_lts(:ncol) = thl700(:ncol) - thl1000(:ncol)
    end if

    !--------------------- Scalar Setting --------------------

    !  Set the ztodt timestep in pbuf for SILHS, this is needed because hdtime is not input to silhs
    ztodtptr = 1.0_kind_phys * hdtime

    !  Determine CLUBB time step and make it sub-step friendly
    !  For now we want CLUBB time step to be 5 min since that is
    !  what has been scientifically validated.  However, there are certain
    !  instances when a 5 min time step will not be possible (based on
    !  host model time step or on macro-micro sub-stepping
    dtime = clubb_timestep

    !  Now check to see if dtime is greater than the host model
    !    (or sub stepped) time step.  If it is, then simply
    !    set it equal to the host (or sub step) time step.
    !    This section is mostly to deal with small host model
    !    time steps (or small sub-steps)
    if (dtime > hdtime) then
      dtime = hdtime
    endif

    !  Now check to see if CLUBB time step divides evenly into
    !    the host model time step.  If not, force it to divide evenly.
    !    We also want it to be 5 minutes or less.  This section is
    !    mainly for host model time steps that are not evenly divisible
    !    by 5 minutes
    if (mod(hdtime,dtime) .ne. 0) then
      dtime = hdtime/2._kind_phys
      do while (dtime > clubb_timestep)
        dtime = dtime/2._kind_phys
      end do
    endif

    !  If resulting host model time step and CLUBB time step do not divide evenly
    !    into each other, have model throw a fit.
    if (mod(hdtime,dtime) .ne. 0) then
      errmsg = 'clubb1_run:  CLUBB time step and HOST time step NOT compatible'
      errflg = 1
      return
    endif

    !  determine number of timesteps CLUBB core should be advanced,
    !  host time step divided by CLUBB time step
    nadv = max(hdtime/dtime,1._kind_phys)

    ! Precalculte the hdtime inverse
    invrs_hdtime = 1._kind_phys / hdtime

    !  Set stats output and increment equal to CLUBB and host dt
    stats_metadata%stats_tsamp = dtime
    stats_metadata%stats_tout  = hdtime

    stats_nsamp = nint(stats_metadata%stats_tsamp/dtime)
    stats_nout = nint(stats_metadata%stats_tout/dtime)

    if (clubb_do_adv) then
      apply_const = 1._kind_phys  ! Initialize to one, only if CLUBB's moments are advected
    else
      apply_const = 0._kind_phys  ! Never want this if CLUBB's moments are not advected
    endif

    ! Initialize the apply_const variable (note special logic is due to eulerian backstepping)
    if (clubb_do_adv .and. (first_step .or. all(wpthlp_pbuf(1:ncol,:)  ==  0._kind_phys))) then
      apply_const = 0._kind_phys  ! On first time through do not remove constant
                           !  from moments since it has not been added yet
    endif

    !----------------------------------------- BEGIN GPU SECTION -----------------------------------------
    ! everything within should be functional with the OpenACC code, or be prevented from running
    ! with using OpenACC, see the "ifdef _OPENACC" section above for restriction examples

    !$acc data copyin( pdf_params_chnk, pdf_params_zm_chnk, sclr_idx, &
    !$acc              state_q, u, v, state_t, state_pmid, &
    !$acc              state_zm, state_phis, state_pdel, state_pdeldry, s, &
    !$acc              state_pint, state_zi, state_omega, lat, &
    !$acc              wsx, wsy, cflx, shf, &
    !$acc              err_info, err_info%err_header, &
    !$acc              cpairv, rairv, se_dis, eleak, cld_pbuf, clubb_params_single_col, grid_dx, grid_dy ) &
    !$acc     copyout( clubb_s, clubbtop_pbuf, &
    !$acc              qclvar, wprcp, rcm_in_layer, rcm, cloud_frac, thlm, rtm, qc_zt, qv_zt, th_zt, &
    !$acc              um, vm, wm_zt, exner, zt_g, zi_g, invrs_cpairv, &
    !$acc              rho_zm, rho_zt, &
    !$acc              pdf_params_chnk%rt_1,                pdf_params_chnk%rt_2,  &
    !$acc              pdf_params_chnk%varnce_rt_1,         pdf_params_chnk%varnce_rt_2, &
    !$acc              pdf_params_chnk%mixt_frac ) &
    !$acc        copy( khzm_pbuf, upwp_pbuf, vpwp_pbuf, up2_pbuf, vp2_pbuf, up3_pbuf, vp3_pbuf, wprtp_pbuf, &
    !$acc              wpthlp_pbuf, wp2_pbuf, wp3_pbuf, rtp2_pbuf, rtp3_pbuf, thlp2_pbuf, thlp3_pbuf, &
    !$acc              rtpthlp_pbuf, wpthvp_pbuf, wp2thvp_pbuf, wp2up_pbuf, ice_supersat_frac_pbuf, &
    !$acc              rtpthvp_pbuf, thlpthvp_pbuf, wp2rtp_pbuf, wp2thlp_pbuf, uprcp_pbuf, vprcp_pbuf, &
    !$acc              rc_coef_zm_pbuf, wp4_pbuf, wpup2_pbuf, wpvp2_pbuf, wp2up2_pbuf, wp2vp2_pbuf ) &
    !$acc      create( um_pert, vm_pert, upwp_pert, vpwp_pert, khzm, &
    !$acc              khzt, thlprcp, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac, &
    !$acc              cloudy_downdraft_frac, cloud_cover, invrs_tau_zm, Lscale, &
    !$acc              invrs_exner_zt, fcor, fcor_y, sfc_elevation, thlm_forcing, rtm_forcing, um_forcing, &
    !$acc              vm_forcing, wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
    !$acc              rtpthlp_forcing, wm_zm, wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, invrs_dz_g, &
    !$acc              p_sfc, upwp_sfc_pert, vpwp_sfc_pert, rtm_ref, thlm_ref, um_ref, vm_ref, &
    !$acc              ug, vg, p_in_Pa, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
    !$acc              invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, rfrzm, clubb_params, deltaz, err_info%err_code, &
    !$acc              pdf_params_chnk%w_1,                 pdf_params_chnk%w_2, &
    !$acc              pdf_params_chnk%varnce_w_1,          pdf_params_chnk%varnce_w_2, &
    !$acc              pdf_params_chnk%thl_1,               pdf_params_chnk%thl_2, &
    !$acc              pdf_params_chnk%varnce_thl_1,        pdf_params_chnk%varnce_thl_2, &
    !$acc              pdf_params_chnk%corr_w_rt_1,         pdf_params_chnk%corr_w_rt_2,  &
    !$acc              pdf_params_chnk%corr_w_thl_1,        pdf_params_chnk%corr_w_thl_2, &
    !$acc              pdf_params_chnk%corr_rt_thl_1,       pdf_params_chnk%corr_rt_thl_2,&
    !$acc              pdf_params_chnk%alpha_thl,           pdf_params_chnk%alpha_rt, &
    !$acc              pdf_params_chnk%crt_1,               pdf_params_chnk%crt_2, &
    !$acc              pdf_params_chnk%cthl_1,              pdf_params_chnk%cthl_2, &
    !$acc              pdf_params_chnk%chi_1,               pdf_params_chnk%chi_2, &
    !$acc              pdf_params_chnk%stdev_chi_1,         pdf_params_chnk%stdev_chi_2, &
    !$acc              pdf_params_chnk%stdev_eta_1,         pdf_params_chnk%stdev_eta_2, &
    !$acc              pdf_params_chnk%covar_chi_eta_1,     pdf_params_chnk%covar_chi_eta_2, &
    !$acc              pdf_params_chnk%corr_w_chi_1,        pdf_params_chnk%corr_w_chi_2, &
    !$acc              pdf_params_chnk%corr_w_eta_1,        pdf_params_chnk%corr_w_eta_2, &
    !$acc              pdf_params_chnk%corr_chi_eta_1,      pdf_params_chnk%corr_chi_eta_2, &
    !$acc              pdf_params_chnk%rsatl_1,             pdf_params_chnk%rsatl_2, &
    !$acc              pdf_params_chnk%rc_1,                pdf_params_chnk%rc_2, &
    !$acc              pdf_params_chnk%cloud_frac_1,        pdf_params_chnk%cloud_frac_2,  &
    !$acc              pdf_params_chnk%ice_supersat_frac_1, pdf_params_chnk%ice_supersat_frac_2 )

    !$acc data if( clubb_config_flags%l_call_pdf_closure_twice ) &
    !$acc        copy( pdf_zm_w_1_pbuf, pdf_zm_w_2_pbuf, pdf_zm_varnce_w_1_pbuf, pdf_zm_varnce_w_2_pbuf, pdf_zm_mixt_frac_pbuf, &
    !$acc              pdf_params_zm_chnk%w_1, pdf_params_zm_chnk%w_2, &
    !$acc              pdf_params_zm_chnk%varnce_w_1, pdf_params_zm_chnk%varnce_w_2, &
    !$acc              pdf_params_zm_chnk%mixt_frac ) &
    !$acc      create( pdf_params_zm_chnk%rt_1,                pdf_params_zm_chnk%rt_2, &
    !$acc              pdf_params_zm_chnk%varnce_rt_1,         pdf_params_zm_chnk%varnce_rt_2, &
    !$acc              pdf_params_zm_chnk%thl_1,               pdf_params_zm_chnk%thl_2, &
    !$acc              pdf_params_zm_chnk%varnce_thl_1,        pdf_params_zm_chnk%varnce_thl_2, &
    !$acc              pdf_params_zm_chnk%corr_w_rt_1,         pdf_params_zm_chnk%corr_w_rt_2, &
    !$acc              pdf_params_zm_chnk%corr_w_thl_1,        pdf_params_zm_chnk%corr_w_thl_2, &
    !$acc              pdf_params_zm_chnk%corr_rt_thl_1,       pdf_params_zm_chnk%corr_rt_thl_2, &
    !$acc              pdf_params_zm_chnk%alpha_thl,           pdf_params_zm_chnk%alpha_rt, &
    !$acc              pdf_params_zm_chnk%crt_1,               pdf_params_zm_chnk%crt_2, &
    !$acc              pdf_params_zm_chnk%cthl_1,              pdf_params_zm_chnk%cthl_2, &
    !$acc              pdf_params_zm_chnk%chi_1,               pdf_params_zm_chnk%chi_2, &
    !$acc              pdf_params_zm_chnk%stdev_chi_1,         pdf_params_zm_chnk%stdev_chi_2, &
    !$acc              pdf_params_zm_chnk%stdev_eta_1,         pdf_params_zm_chnk%stdev_eta_2, &
    !$acc              pdf_params_zm_chnk%covar_chi_eta_1,     pdf_params_zm_chnk%covar_chi_eta_2, &
    !$acc              pdf_params_zm_chnk%corr_w_chi_1,        pdf_params_zm_chnk%corr_w_chi_2, &
    !$acc              pdf_params_zm_chnk%corr_w_eta_1,        pdf_params_zm_chnk%corr_w_eta_2, &
    !$acc              pdf_params_zm_chnk%corr_chi_eta_1,      pdf_params_zm_chnk%corr_chi_eta_2, &
    !$acc              pdf_params_zm_chnk%rsatl_1,             pdf_params_zm_chnk%rsatl_2, &
    !$acc              pdf_params_zm_chnk%rc_1,                pdf_params_zm_chnk%rc_2, &
    !$acc              pdf_params_zm_chnk%cloud_frac_1,        pdf_params_zm_chnk%cloud_frac_2, &
    !$acc              pdf_params_zm_chnk%ice_supersat_frac_1, pdf_params_zm_chnk%ice_supersat_frac_2 )

    !$acc data if( sclr_dim > 0 ) &
    !$acc      create( wpsclrp_sfc, sclrm_forcing, sclrm, wpsclrp, sclrp2, sclrp3, sclrprtp, sclrpthlp, sclrpthvp ) &
    !$acc      copyin( sclr_tol )

    !$acc data if( edsclr_dim > 0 ) &
    !$acc     copyout( edsclr ) &
    !$acc      create( wpedsclrp_sfc, edsclrm_forcing )

    !$acc data if( hydromet_dim > 0 ) &
    !$acc      copyin( hm_metadata, hm_metadata%l_mix_rat_hm ) &
    !$acc      create( wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt )

    !----------------------------------------- Zeroing -----------------------------------------

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol

        !  Define forcings from CAM to CLUBB as zero for momentum and thermo,
        !  forcings already applied through CAM
        thlm_forcing(i,k)   = 0._kind_phys
        rtm_forcing(i,k)    = 0._kind_phys
        um_forcing(i,k)     = 0._kind_phys
        vm_forcing(i,k)     = 0._kind_phys

        rtm_ref(i,k)        = 0.0_kind_phys
        thlm_ref(i,k)       = 0.0_kind_phys
        um_ref(i,k)         = 0.0_kind_phys
        vm_ref(i,k)         = 0.0_kind_phys
        ug(i,k)             = 0.0_kind_phys
        vg(i,k)             = 0.0_kind_phys

        ! Perturbed winds are not used in CAM
        um_pert(i,k)        = 0.0_kind_phys
        vm_pert(i,k)        = 0.0_kind_phys
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm_clubb
      do i = 1, ncol
        ! Perturbed winds are not used in CAM
        upwp_pert(i,k)      = 0.0_kind_phys
        vpwp_pert(i,k)      = 0.0_kind_phys
      end do
    end do

    !$acc parallel loop gang vector default(present)
    do i = 1, ncol
      ! Perturbed winds are not used in CAM
      upwp_sfc_pert(i) = 0.0_kind_phys
      vpwp_sfc_pert(i) = 0.0_kind_phys

      ! When run in host models, CLUBB does not apply Coriolis tendencies to the
      ! mean horizontal wind components (this is controlled by the `l_implemented`
      ! flag, which should be hardcoded to .true. in this file).
      !
      ! However, enabling `clubb_l_ho_nontrad_coriolis` or `clubb_l_ho_trad_coriolis`
      ! introduces Coriolis effects in higher-order moments (e.g., wp2up).
      ! Therefore, we still compute the Coriolis parameters here for potential
      ! use by those higher-order terms.
      fcor(i)   = 2._kind_phys * omega_const * sin( lat(i) )
      fcor_y(i) = 2._kind_phys * omega_const * cos( lat(i) )
    end do

    if ( sclr_dim > 0 ) then
      !  higher order scalar stuff, put to zero
      !$acc parallel loop gang vector collapse(3) default(present)
      do sclr = 1, sclr_dim
        do k = 1, nzt_clubb
          do i = 1, ncol
            sclrm(i,k,sclr)           = 0._kind_phys
            sclrp3(i,k,sclr)          = 0._kind_phys
            sclrm_forcing(i,k,sclr)   = 0._kind_phys
          end do
        end do
      end do

      !  higher order scalar stuff, put to zero
      !$acc parallel loop gang vector collapse(3) default(present)
      do sclr = 1, sclr_dim
        do k = 1, nzm_clubb
          do i = 1, ncol
            wpsclrp(i,k,sclr)         = 0._kind_phys
            sclrp2(i,k,sclr)          = 0._kind_phys
            sclrprtp(i,k,sclr)        = 0._kind_phys
            sclrpthlp(i,k,sclr)       = 0._kind_phys
            sclrpthvp(i,k,sclr)       = 0._kind_phys
          end do
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do sclr = 1, sclr_dim
        do i = 1, ncol
          wpsclrp_sfc(i,sclr) = 0._kind_phys
        end do
      end do
    end if

    if ( edsclr_dim > 0 ) then
      !$acc parallel loop gang vector collapse(3) default(present)
      do iedsclr = 1, edsclr_dim
        do k = 1, nzt_clubb
          do i = 1, ncol
            edsclrm_forcing(i,k,iedsclr) = 0._kind_phys
          end do
        end do
      end do

      !  Define surface sources for transported variables for diffusion, will
      !  be zero as these tendencies are done in vertical_diffusion
      !$acc parallel loop gang vector collapse(2) default(present)
      do iedsclr = 1, edsclr_dim
        do i = 1, ncol
          wpedsclrp_sfc(i,iedsclr) = 0._kind_phys
        end do
      end do
    end if

    if ( hydromet_dim > 0 ) then

      !$acc parallel loop gang vector collapse(3) default(present)
      do ixind = 1, hydromet_dim
        do k = 1, nzt_clubb
          do i = 1, ncol
            wp2hmp(i,k,ixind)      = 0._kind_phys
            rtphmp_zt(i,k,ixind)   = 0._kind_phys
            thlphmp_zt(i,k,ixind)  = 0._kind_phys
          end do
        end do
      end do

      !$acc parallel loop gang vector collapse(3) default(present)
      do ixind = 1, hydromet_dim
        do k = 1, nzm_clubb
          do i = 1, ncol
            wphydrometp(i,k,ixind) = 0._kind_phys
          end do
        end do
      end do

    end if

   if (do_clubb_mf) then
      mf_L0          = 0._kind_phys
      mf_L0_nadv     = 0._kind_phys
      mf_ztop        = 0._kind_phys
      mf_ztop_nadv   = 0._kind_phys
      mf_ztopm1      = 0._kind_phys
      mf_ztopm1_nadv = 0._kind_phys
      mf_cape        = 0._kind_phys
      mf_cape_nadv   = 0._kind_phys
      mf_ddcp_nadv   = 0._kind_phys
      mf_cbm1        = 0._kind_phys
      mf_cbm1_nadv   = 0._kind_phys
      mf_freq_nadv   = 0._kind_phys

      if (macmic_it==1) then
         ztopm1_macmic(:ncol,:) = 0._kind_phys
         ddcp_macmic(:ncol,:) = 0._kind_phys
         cbm1_macmic(:ncol) = 0._kind_phys
         s_mfup_macmic(:ncol,:)   = 0._kind_phys
         s_entup_macmic(:ncol,:)  = 0._kind_phys
         s_detup_macmic(:ncol,:)  = 0._kind_phys
         s_mfdn_macmic(:ncol,:)   = 0._kind_phys
         s_entdn_macmic(:ncol,:)  = 0._kind_phys
         s_detdn_macmic(:ncol,:)  = 0._kind_phys
         mf_sqtac_macmic(:ncol,:) = 0._kind_phys
         mf_sqtev_macmic(:ncol,:) = 0._kind_phys
         mf_kctop_macmic(:ncol)   = 0._kind_phys
         s_awwup_macmic(:ncol,:)  = 0._kind_phys
         s_awup_macmic(:ncol,:)   = 0._kind_phys
      end if

      mf_precc_nadv(:ncol)      = 0._kind_phys
      mf_snow_nadv(:ncol)       = 0._kind_phys

      mf_qc(:ncol,:nzm_clubb)             = 0._kind_phys
      mf_rcm(:ncol,:nzm_clubb)            = 0._kind_phys
      mf_cloudfrac(:ncol,:nzm_clubb)      = 0._kind_phys
      mf_qc_nadv(:ncol,:nzm_clubb)        = 0._kind_phys
      mf_rcm_nadv(:ncol,:nzm_clubb)       = 0._kind_phys
      mf_ent_nadv(:ncol,:nzm_clubb)       = 0._kind_phys
      mf_cloudfrac_nadv(:ncol,:nzm_clubb) = 0._kind_phys

      mf_thlforcup_nadv(:ncol,:nzt_clubb) = 0._kind_phys
      mf_qtforcup_nadv(:ncol,:nzt_clubb)  = 0._kind_phys
      mf_thlforcdn_nadv(:ncol,:nzt_clubb) = 0._kind_phys
      mf_qtforcdn_nadv(:ncol,:nzt_clubb)  = 0._kind_phys
      mf_thlforc_nadv(:ncol,:nzt_clubb)   = 0._kind_phys
      mf_qtforc_nadv(:ncol,:nzt_clubb)    = 0._kind_phys

      s_mfup_nadv(:ncol,:nzm_clubb)   = 0._kind_phys
      s_entup_nadv(:ncol,:nzm_clubb)  = 0._kind_phys
      s_detup_nadv(:ncol,:nzm_clubb)  = 0._kind_phys
      s_mfdn_nadv(:ncol,:nzm_clubb)   = 0._kind_phys
      s_entdn_nadv(:ncol,:nzm_clubb)  = 0._kind_phys
      s_detdn_nadv(:ncol,:nzm_clubb)  = 0._kind_phys
      s_awwup_nadv(:ncol,:nzm_clubb)  = 0._kind_phys
      s_awup_nadv(:ncol,:nzm_clubb)   = 0._kind_phys
      mf_sqtac_nadv(:ncol,:nzt_clubb) = 0._kind_phys
      mf_sqtev_nadv(:ncol,:nzt_clubb) = 0._kind_phys
      mf_kctop_nadv(:ncol)            = 0._kind_phys

      max_cfl_nadv(:ncol) = 0._kind_phys
   end if

    !----------------------------------------- Initializing arrays -----------------------------------------

    if ( clubb_do_adv ) then

      if (macmic_it  ==  1) then

        !  Note that some of the moments below can be positive or negative.
        !    Remove a constant that was added to prevent dynamics from clipping
        !    them to prevent dynamics from making them positive.
        do k = 1, nzt_clubb
          do i = 1, ncol
            k_cam = top_lev - 1 + k
            rtpthlp_pbuf(i,k) = state_q(i,k_cam,ixrtpthlp) - ( rtpthlp_const * apply_const )
            wpthlp_pbuf(i,k)  = state_q(i,k_cam, ixwpthlp) - ( wpthlp_const  * apply_const )
            wprtp_pbuf(i,k)   = state_q(i,k_cam,  ixwprtp) - ( wprtp_const   * apply_const )
            wp3_pbuf(i,k)     = state_q(i,k_cam,    ixwp3) - ( wp3_const     * apply_const )
            wp2_pbuf(i,k)     = max(  w_tol_sqd, state_q(i,k_cam,    ixwp2) )
            thlp2_pbuf(i,k)   = max( thl_tol**2, state_q(i,k_cam,  ixthlp2) )
            rtp2_pbuf(i,k)    = max(  rt_tol**2, state_q(i,k_cam,   ixrtp2) )
            up2_pbuf(i,k)     = max(  w_tol_sqd, state_q(i,k_cam,    ixup2) )
            vp2_pbuf(i,k)     = max(  w_tol_sqd, state_q(i,k_cam,    ixvp2) )
          enddo
        enddo

      endif

      ! If not last step of macmic loop then set apply_const back to
      !   zero to prevent output from being corrupted.
      if (macmic_it  ==  cld_macmic_num_steps) then
        apply_const = 1._kind_phys
      else
        apply_const = 0._kind_phys
      endif

    endif

    !$acc parallel loop gang vector collapse(2) default(present)
    do n = 1, nparams
      do i = 1, ncol
        clubb_params(i,n) = clubb_params_single_col(1,n)
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, pver
      do i = 1, ncol
        invrs_cpairv(i,k) = 1._kind_phys / cpairv(i,k)
      end do
    end do

    !  Compute thermodynamic stuff needed for CLUBB on thermo levels.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol

        k_cam = top_lev - 1 + k

        ! Define the CLUBB thermodynamic grid (in units of m)
        zt_g(i,k) = state_zm(i,k_cam) - state_zi(i,pverp)

        invrs_dz_g(i,k) = 1._kind_phys / ( state_zi(i,k_cam) - state_zi(i,k_cam+1) )  ! compute thickness

        rho_zt(i,k)          = rga * state_pdel(i,k_cam)    * invrs_dz_g(i,k)

        rho_ds_zt(i,k)       = rga * state_pdeldry(i,k_cam) * invrs_dz_g(i,k)

        invrs_rho_ds_zt(i,k) = 1._kind_phys / rho_ds_zt(i,k)

      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol

        k_cam = top_lev - 1 + k

        p_in_Pa(i,k) = state_pmid(i,k_cam)

        !  Compute inverse exner function consistent with CLUBB's definition, which uses a constant
        !  surface pressure.  CAM's exner (in state) does not.  Therefore, for consistent
        !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables
        !  (such as thlm), use "invrs_exner_zt" otherwise use the exner in state
        exner(i,k) = ( p_in_Pa(i,k) * inv_p0_clubb )**( rairv(i,k_cam) * invrs_cpairv(i,k_cam) )

        invrs_exner_zt(i,k) = 1._kind_phys / exner(i,k)

        ! exception - setting this to moist thv_ds_zt
        thv_ds_zt(i,k) = state_t(i,k_cam) * invrs_exner_zt(i,k)  &
                         * (1._kind_phys + zvir * state_q(i,k_cam,ixq) - state_q(i,k_cam,ixcldliq))

        rcm(i,k)    = state_q(i,k_cam,ixcldliq)
        rtm(i,k)    = state_q(i,k_cam,ixq) + state_q(i,k_cam,ixcldliq)

        thlm(i,k)   = ( state_t(i,k_cam) - ( latvap * invrs_cpairv(i,k_cam) ) &
                                               * state_q(i,k_cam,ixcldliq) ) * invrs_exner_zt(i,k)
        qc_zt(i,k) = state_q(i,k_cam,ixcldliq)
        qv_zt(i,k) = state_q(i,k_cam,ixq)
        th_zt(i,k) = state_t(i,k_cam)*invrs_exner_zt(i,k)
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol

        k_cam = top_lev - 1 + k

        !  Compute mean w wind on thermo grid, convert from omega to w
        wm_zt(i,k) = -1._kind_phys * ( state_omega(i,k_cam) - state_omega(i,pver) ) / ( rho_zt(i,k) * gravit )

        cloud_frac(i,k)       = cld_pbuf(i,k_cam)

        um(i,k) = u(i,k_cam)
        vm(i,k) = v(i,k_cam)

        rfrzm(i,k)  = state_q(i,k_cam,ixcldice)

      end do
    end do

    !$acc parallel loop gang vector default(present)
    do i = 1, ncol

      deltaz(i)         = state_zi(i,pverp-1) - state_zi(i,pverp)

      !  Set the surface pressure
      p_sfc(i)          = state_pint(i,pverp)

      !  Set the elevation of the surface
      sfc_elevation(i)  = state_zi(i,pverp)

    end do

    ! Define the CLUBB momentum grid (in height, units of m)
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm_clubb
      do i = 1, ncol
        k_cam = top_lev - 1 + k
        zi_g(i,k) = state_zi(i,k_cam) - state_zi(i,pverp)
      end do
    end do

    if (do_clubb_mf) then

      do k = 1, nzt_clubb
        do i = 1, ncol
          k_cam = top_lev - 1 + k
          kappa_zt(i,k)   = rairv(i,k_cam) * invrs_cpairv(i,k_cam)
          dz_g(i,k)       = state_zi(i,k_cam) - state_zi(i,k_cam+1)  ! compute thickness
        end do
      end do

      ! pressure on momentum grid needed for mass flux calc.
      do k = 1, nzm_clubb
        do i = 1, ncol
          k_cam = top_lev - 1 + k
          p_in_Pa_zm(i,k)     = state_pint(i,k_cam)
          tke_zm(i,k)     = tke_pbuf(i,k_cam)
        end do
      end do

    end if

    !----------------------------------------- Initializing CLUBB grid -----------------------------------------
    ! Note: these few routines, setup_grid_api, calc_derrived_params_api, and check_parameters_api are not
    !       GPUized yet, so we need to copy data to and from the GPU.

    !  Heights need to be set at each timestep.  Therefore, recall
    !  setup_grid and calc_derrived_params for this.
    !  IMPORTANT NOTE:  do not make any calls that use CLUBB grid-height
    !                   operators (such as zt2zm_api, etc.) until AFTER the
    !                   call to setup_grid_heights_api.
    !$acc update host( deltaz, zi_g, zt_g, clubb_params, sfc_elevation )

    ! Calculate grid assuming a descending grid (cam grid), since we want to
    ! confine ascending behavior to advance_clubb_core
    call setup_grid_api( nzm_clubb, ncol, sfc_elevation, l_implemented,   & ! intent(in)
                         .false., grid_type,                              & ! intent(in)
                         deltaz, zi_g(:,nzm_clubb), zi_g(:,1),            & ! intent(in)
                         zi_g, zt_g,                                      & ! intent(in)
                         gr, err_info )                                     ! intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) ) then
       errmsg = 'clubb1_run: '//err_info%err_header_global//NEW_LINE('a')// &
                   'in CLUBB setup_grid'
       errflg = 1
       goto 200
    end if

    call calc_derrived_params_api( gr, ncol, grid_type, deltaz,                 & ! Intent(in)
                                   clubb_params,                                & ! Intent(in)
                                   clubb_config_flags%l_prescribed_avg_deltaz,  & ! Intent(in)
                                   nu_vert_res_dep, lmin,                       & ! intent(inout)
                                   mixt_frac_max_mag )                            ! intent(inout)

    call check_parameters_api( ncol, clubb_params, lmin, & ! Intent(in)
                               err_info )                  ! Intent(inout)

    if ( any(err_info%err_code == clubb_fatal_error) ) then
       errmsg = 'clubb1_run: '//err_info%err_header_global//NEW_LINE('a')// &
                   'in CLUBB check_parameters_api'
       errflg = 1
       goto 200
    end if

    ! CLUBB's grid data structure (gr) and nu_vert_res_dep contain arrays that need to
    ! be copied to the GPU
    !$acc data copyin( gr, gr%zm, gr%zt, gr%dzm, gr%dzt, gr%invrs_dzt, gr%invrs_dzm, &
    !$acc              gr%weights_zt2zm, gr%weights_zm2zt, &
    !$acc              nu_vert_res_dep, nu_vert_res_dep%nu2, nu_vert_res_dep%nu9, &
    !$acc              nu_vert_res_dep%nu1, nu_vert_res_dep%nu8, nu_vert_res_dep%nu10, &
    !$acc              nu_vert_res_dep%nu6)
    !----------------------------------------- END CLUBB grid initialization -----------------------------------------

#ifdef SILHS
    ! Add forcings for SILHS covariance contributions
    rtp2_forcing    = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,    rtp2_mc_zt_pbuf(1:ncol,:) )
    thlp2_forcing   = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,   thlp2_mc_zt_pbuf(1:ncol,:) )
    wprtp_forcing   = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,   wprtp_mc_zt_pbuf(1:ncol,:) )
    wpthlp_forcing  = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,  wpthlp_mc_zt_pbuf(1:ncol,:) )
    rtpthlp_forcing = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, rtpthlp_mc_zt_pbuf(1:ncol,:) )

    ! Zero out SILHS covariance contribution terms
    rtp2_mc_zt_pbuf     = 0.0_kind_phys
    thlp2_mc_zt_pbuf    = 0.0_kind_phys
    wprtp_mc_zt_pbuf    = 0.0_kind_phys
    wpthlp_mc_zt_pbuf   = 0.0_kind_phys
    rtpthlp_mc_zt_pbuf  = 0.0_kind_phys
#else
    ! Set forcings to zero if not using SILHS
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm_clubb
      do i = 1, ncol
        rtp2_forcing(i,k)    = 0._kind_phys
        thlp2_forcing(i,k)   = 0._kind_phys
        wprtp_forcing(i,k)   = 0._kind_phys
        wpthlp_forcing(i,k)  = 0._kind_phys
        rtpthlp_forcing(i,k) = 0._kind_phys
      end do
    end do
#endif

    ! Compute some inputs from the thermodynamic grid to the momentum grid
    rho_ds_zm       = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, rho_ds_zt )
    invrs_rho_ds_zm = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, invrs_rho_ds_zt )
    rho_zm          = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, rho_zt )
    thv_ds_zm       = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, thv_ds_zt )
    wm_zm           = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, wm_zt )

    ! Surface fluxes provided by host model
    !$acc parallel loop gang vector default(present)
    do i = 1, ncol
      wpthlp_sfc(i) = shf(i) / ( cpairv(i,pver) * rho_ds_zm(i,nzm_clubb) ) & ! Sensible heat flux
                      * invrs_exner_zt(i,nzt_clubb)
      wprtp_sfc(i)  = cflx(i,1) / rho_ds_zm(i,nzm_clubb)                           ! Moisture flux
    end do


    ! ------------------------------------------------- !
    ! Begin case specific code for SCAM cases.          !
    ! This section of code block is NOT called in       !
    ! global simulations                                !
    ! ------------------------------------------------- !
    if (single_column .and. .not. scm_cambfb_mode) then

      !  Initialize zo if variable ustar is used
      if (landfrac(1) >= 0.5_kind_phys) then
        zo = 0.035_kind_phys
      else
        zo = 0.0001_kind_phys
      endif

      !  Compute surface wind (ubar)
      ubar = sqrt(um(1,nzt_clubb)**2+vm(1,nzt_clubb)**2)
      if (ubar <  0.25_kind_phys) ubar = 0.25_kind_phys

      !  Below denotes case specifics for surface momentum
      !  and thermodynamic fluxes, depending on the case

      !  Define ustar (based on case, if not variable)
      ustar = 0.25_kind_phys   ! Initialize ustar in case no case

      if(trim(scm_clubb_iop_name)  ==  'BOMEX_5day') then
        ustar = 0.28_kind_phys
      endif

      if(trim(scm_clubb_iop_name)  ==  'ATEX_48hr') then
        ustar = 0.30_kind_phys
      endif

      if(trim(scm_clubb_iop_name)  ==  'RICO_3day') then
        ustar      = 0.28_kind_phys
      endif

      if(trim(scm_clubb_iop_name)  ==  'arm97' .or. trim(scm_clubb_iop_name)  ==  'gate' .or. &
         trim(scm_clubb_iop_name)  ==  'toga' .or. trim(scm_clubb_iop_name)  ==  'mpace' .or. &
         trim(scm_clubb_iop_name)  ==  'ARM_CC') then

          bflx22 = (gravit/theta0)*wpthlp_sfc(1)
          ustar  = diag_ustar(zt_g(1,nzt_clubb),bflx22,ubar,zo,shr_const_karman,shr_const_pi,shr_const_g)
      endif

      !  Compute the surface momentum fluxes, if this is a SCAM simulation
      upwp_sfc(1) = -um(1,nzt_clubb)*ustar**2/ubar
      vpwp_sfc(1) = -vm(1,nzt_clubb)*ustar**2/ubar

    end if

    ! Implementation after Thomas Toniazzo (NorESM) and Colin Zarzycki (PSU)
    !  Other Surface fluxes provided by host model
    if( (cld_macmic_num_steps > 1) .and. clubb_l_intr_sfc_flux_smooth ) then

      !$acc update host( u, v, state_t, state_pmid, wsx, wsy )

      ! Adjust surface stresses using winds from the prior macmic iteration
      do i = 1, ncol
        ubar = sqrt(u(i,pver)**2+v(i,pver)**2)
        if (ubar <  0.25_kind_phys) ubar = 0.25_kind_phys

        rrho(i)     = calc_ideal_gas_rrho(rair, state_t(i,pver), state_pmid(i,pver))
        ustar    = calc_friction_velocity(wsx(i), wsy(i), rrho(i))

        upwp_sfc(i) = -u(i,pver)*ustar**2/ubar
        vpwp_sfc(i) = -v(i,pver)*ustar**2/ubar
      end do

      !$acc update device( upwp_sfc, vpwp_sfc )

    else

      !$acc parallel loop gang vector default(present)
      do i = 1, ncol
        upwp_sfc(i)   = wsx(i) / rho_ds_zm(i,nzm_clubb)               ! Surface meridional momentum flux
        vpwp_sfc(i)   = wsy(i) / rho_ds_zm(i,nzm_clubb)               ! Surface zonal momentum flux
      end do

    endif

    ! We only need to copy pdf_params from pbuf if this is a restart, we're calling pdf_closure
    ! at the end of advance_clubb_core, and calling it twice for pdf_params_zm as well
    if ( first_restart_step &
         .and. clubb_config_flags%l_call_pdf_closure_twice &
         .and. clubb_config_flags%ipdf_call_placement .eq. ipdf_post_advance_fields ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm_clubb
        do i = 1, ncol
          pdf_params_zm_chnk%w_1(i,k)        = pdf_zm_w_1_pbuf(i,k)
          pdf_params_zm_chnk%w_2(i,k)        = pdf_zm_w_2_pbuf(i,k)
          pdf_params_zm_chnk%varnce_w_1(i,k) = pdf_zm_varnce_w_1_pbuf(i,k)
          pdf_params_zm_chnk%varnce_w_2(i,k) = pdf_zm_varnce_w_2_pbuf(i,k)
          pdf_params_zm_chnk%mixt_frac(i,k)  = pdf_zm_mixt_frac_pbuf(i,k)
        end do
      end do

    end if

    if ( edsclr_dim > 0 ) then

      !  Copy the cam version of the tracers to the clubb version
      ! NOTE: if clubb_l_do_expldiff_rtm_thlm=.true., then the last two
      !       tracers are thlm and rtm, which are added inside clubb
      icnt=0
      do ixind = 1, pcnst
        if (lq(ixind))  then

          icnt = icnt+1

          !$acc parallel loop gang vector collapse(2) default(present)
          do k = 1, nzt_clubb
            do i = 1, ncol
              k_cam = top_lev - 1 + k
              edsclr(i,k,icnt)       = state_q(i,k_cam,ixind)
            end do
          end do

        end if
      end do

    end if

    !----------------------------------------- Substepping loop -----------------------------------------
    do tt = 1, nadv    ! do needed number of "sub" timesteps for each CAM step

      !  Increment the statistics then begin stats timestep
      if (stats_metadata%l_stats) then
        call stats_begin_timestep_api( tt, stats_nsamp, stats_nout, &
                                       stats_metadata )
      endif

      if ( clubb_l_ascending_grid ) then

        ! CLUBB is to be run in ascending mode, which has the surface at k=1, which is
        ! the opposite of the cam grid that the rest of clubb_intr uses, so
        ! we need to flip the fields (in the vertical dimensions) before calling advance_clubb_core
        !
        ! NOTE: We do not neccesarily flip all arrays, only ones that are used within this
        !       subroutine (advance_clubb_core). For example, only the pdf_params fields that
        !       are used within this subroutine (or used in a subroutine we call) need to
        !       be flipped.

        wm_zt                     =                     wm_zt(:,nzt_clubb:1:-1)
        rho_zt                    =                    rho_zt(:,nzt_clubb:1:-1)
        rho_ds_zt                 =                 rho_ds_zt(:,nzt_clubb:1:-1)
        invrs_rho_ds_zt           =           invrs_rho_ds_zt(:,nzt_clubb:1:-1)
        thv_ds_zt                 =                 thv_ds_zt(:,nzt_clubb:1:-1)
        khzt                      =                      khzt(:,nzt_clubb:1:-1)
        rtm_ref                   =                   rtm_ref(:,nzt_clubb:1:-1)
        thlm_ref                  =                  thlm_ref(:,nzt_clubb:1:-1)
        um_ref                    =                    um_ref(:,nzt_clubb:1:-1)
        vm_ref                    =                    vm_ref(:,nzt_clubb:1:-1)
        ug                        =                        ug(:,nzt_clubb:1:-1)
        vg                        =                        vg(:,nzt_clubb:1:-1)
        p_in_Pa                   =                   p_in_Pa(:,nzt_clubb:1:-1)
        exner                     =                     exner(:,nzt_clubb:1:-1)
        rfrzm                     =                     rfrzm(:,nzt_clubb:1:-1)
        um                        =                        um(:,nzt_clubb:1:-1)
        vm                        =                        vm(:,nzt_clubb:1:-1)
        up3_pbuf                  =                  up3_pbuf(:,nzt_clubb:1:-1)
        vp3_pbuf                  =                  vp3_pbuf(:,nzt_clubb:1:-1)
        wp3_pbuf                  =                  wp3_pbuf(:,nzt_clubb:1:-1)
        rtp3_pbuf                 =                 rtp3_pbuf(:,nzt_clubb:1:-1)
        thlp3_pbuf                =                thlp3_pbuf(:,nzt_clubb:1:-1)
        rcm                       =                       rcm(:,nzt_clubb:1:-1)
        cloud_frac                =                cloud_frac(:,nzt_clubb:1:-1)
        wpup2_pbuf                =                wpup2_pbuf(:,nzt_clubb:1:-1)
        wpvp2_pbuf                =                wpvp2_pbuf(:,nzt_clubb:1:-1)
        wp2rtp_pbuf               =               wp2rtp_pbuf(:,nzt_clubb:1:-1)
        wp2thlp_pbuf              =              wp2thlp_pbuf(:,nzt_clubb:1:-1)
        ice_supersat_frac_pbuf    =    ice_supersat_frac_pbuf(:,nzt_clubb:1:-1)
        um_pert                   =                   um_pert(:,nzt_clubb:1:-1)
        vm_pert                   =                   vm_pert(:,nzt_clubb:1:-1)
        wp2thvp_pbuf              =              wp2thvp_pbuf(:,nzt_clubb:1:-1)
        wp2up_pbuf                =                wp2up_pbuf(:,nzt_clubb:1:-1)
        rtm                       =                       rtm(:,nzt_clubb:1:-1)
        thlm                      =                      thlm(:,nzt_clubb:1:-1)

        ! flipping grid definition
        zt_g                      =                      zt_g(:,nzt_clubb:1:-1)

        if (do_clubb_mf) then
          dz_g                    =                      dz_g(:,nzt_clubb:1:-1)
          invrs_dz_g              =                invrs_dz_g(:,nzt_clubb:1:-1)
          th_zt                   =                     th_zt(:,nzt_clubb:1:-1)
          qv_zt                   =                     qv_zt(:,nzt_clubb:1:-1)
          qc_zt                   =                     qc_zt(:,nzt_clubb:1:-1)
          kappa_zt                =                  kappa_zt(:,nzt_clubb:1:-1)
          invrs_exner_zt          =            invrs_exner_zt(:,nzt_clubb:1:-1)
        end if

        wprtp_forcing             =             wprtp_forcing(:,nzm_clubb:1:-1)
        wpthlp_forcing            =            wpthlp_forcing(:,nzm_clubb:1:-1)
        rtp2_forcing              =              rtp2_forcing(:,nzm_clubb:1:-1)
        thlp2_forcing             =             thlp2_forcing(:,nzm_clubb:1:-1)
        rtpthlp_forcing           =           rtpthlp_forcing(:,nzm_clubb:1:-1)
        wm_zm                     =                     wm_zm(:,nzm_clubb:1:-1)
        rho_zm                    =                    rho_zm(:,nzm_clubb:1:-1)
        rho_ds_zm                 =                 rho_ds_zm(:,nzm_clubb:1:-1)
        invrs_rho_ds_zm           =           invrs_rho_ds_zm(:,nzm_clubb:1:-1)
        thv_ds_zm                 =                 thv_ds_zm(:,nzm_clubb:1:-1)
        upwp_pbuf                 =                 upwp_pbuf(:,nzm_clubb:1:-1)
        vpwp_pbuf                 =                 vpwp_pbuf(:,nzm_clubb:1:-1)

        ! flipping grid definition
        zi_g                      =                      zi_g(:,nzm_clubb:1:-1)

        if (do_clubb_mf) then
          p_in_Pa_zm              =                p_in_Pa_zm(:,nzm_clubb:1:-1)
          tke_zm                  =                    tke_zm(:,nzm_clubb:1:-1)
        end if

        up2_pbuf                  =                  up2_pbuf(:,nzm_clubb:1:-1)
        vp2_pbuf                  =                  vp2_pbuf(:,nzm_clubb:1:-1)
        wprtp_pbuf                =                wprtp_pbuf(:,nzm_clubb:1:-1)
        wpthlp_pbuf               =               wpthlp_pbuf(:,nzm_clubb:1:-1)
        wp2_pbuf                  =                  wp2_pbuf(:,nzm_clubb:1:-1)
        rtp2_pbuf                 =                 rtp2_pbuf(:,nzm_clubb:1:-1)
        thlp2_pbuf                =                thlp2_pbuf(:,nzm_clubb:1:-1)
        rtpthlp_pbuf              =              rtpthlp_pbuf(:,nzm_clubb:1:-1)
        wpthvp_pbuf               =               wpthvp_pbuf(:,nzm_clubb:1:-1)
        rtpthvp_pbuf              =              rtpthvp_pbuf(:,nzm_clubb:1:-1)
        thlpthvp_pbuf             =             thlpthvp_pbuf(:,nzm_clubb:1:-1)
        uprcp_pbuf                =                uprcp_pbuf(:,nzm_clubb:1:-1)
        vprcp_pbuf                =                vprcp_pbuf(:,nzm_clubb:1:-1)
        rc_coef_zm_pbuf           =           rc_coef_zm_pbuf(:,nzm_clubb:1:-1)
        wp4_pbuf                  =                  wp4_pbuf(:,nzm_clubb:1:-1)
        wp2up2_pbuf               =               wp2up2_pbuf(:,nzm_clubb:1:-1)
        wp2vp2_pbuf               =               wp2vp2_pbuf(:,nzm_clubb:1:-1)
        upwp_pert                 =                 upwp_pert(:,nzm_clubb:1:-1)
        vpwp_pert                 =                 vpwp_pert(:,nzm_clubb:1:-1)

        if ( edsclr_dim > 0 ) then
          edsclr          =          edsclr(:,nzt_clubb:1:-1,:)
          edsclrm_forcing = edsclrm_forcing(:,nzt_clubb:1:-1,:)
        end if

        if ( sclr_dim > 0 ) then

          sclrm_forcing    =    sclrm_forcing(:,nzt_clubb:1:-1,:)
          sclrm            =            sclrm(:,nzt_clubb:1:-1,:)
          sclrp3           =           sclrp3(:,nzt_clubb:1:-1,:)

          sclrp2           =           sclrp2(:,nzm_clubb:1:-1,:)
          sclrprtp         =         sclrprtp(:,nzm_clubb:1:-1,:)
          sclrpthlp        =        sclrpthlp(:,nzm_clubb:1:-1,:)
          wpsclrp          =          wpsclrp(:,nzm_clubb:1:-1,:)
          sclrpthvp        =        sclrpthvp(:,nzm_clubb:1:-1,:)
        end if

        ! These are flipped, ensuring these are stored in descending mode, regardless of clubb_l_ascending_grid.
        ! only because these are need to be stored for restarts
        if ( clubb_config_flags%l_call_pdf_closure_twice ) then
          pdf_params_zm_chnk%w_1        = pdf_params_zm_chnk%w_1       (:,nzm_clubb:1:-1)
          pdf_params_zm_chnk%w_2        = pdf_params_zm_chnk%w_2       (:,nzm_clubb:1:-1)
          pdf_params_zm_chnk%varnce_w_1 = pdf_params_zm_chnk%varnce_w_1(:,nzm_clubb:1:-1)
          pdf_params_zm_chnk%varnce_w_2 = pdf_params_zm_chnk%varnce_w_2(:,nzm_clubb:1:-1)
          pdf_params_zm_chnk%mixt_frac  = pdf_params_zm_chnk%mixt_frac (:,nzm_clubb:1:-1)
        end if

        ! These are flipped, ensuring these are stored in descending mode, regardless of clubb_l_ascending_grid.
        ! only for pdfp_rtp2_output calc
        pdf_params_chnk%mixt_frac    = pdf_params_chnk%mixt_frac  (:,nzt_clubb:1:-1)
        pdf_params_chnk%rt_1         = pdf_params_chnk%rt_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk%rt_2         = pdf_params_chnk%rt_2       (:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_rt_1  = pdf_params_chnk%varnce_rt_1(:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_rt_2  = pdf_params_chnk%varnce_rt_2(:,nzt_clubb:1:-1)

        ! These are flipped, ensuring these are stored in descending mode, regardless of clubb_l_ascending_grid.
        ! only for update_xp2_mc_api call
        pdf_params_chnk%w_1          = pdf_params_chnk%w_1         (:,nzt_clubb:1:-1)
        pdf_params_chnk%w_2          = pdf_params_chnk%w_2         (:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_w_1   = pdf_params_chnk%varnce_w_1  (:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_w_2   = pdf_params_chnk%varnce_w_2  (:,nzt_clubb:1:-1)
        pdf_params_chnk%thl_1        = pdf_params_chnk%thl_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk%thl_2        = pdf_params_chnk%thl_2       (:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_thl_1 = pdf_params_chnk%varnce_thl_1(:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_thl_2 = pdf_params_chnk%varnce_thl_2(:,nzt_clubb:1:-1)

        ! These are flipped for silhs, which uses a cam grid
        pdf_params_chnk%rc_1                 = pdf_params_chnk%rc_1               (:,nzt_clubb:1:-1)
        pdf_params_chnk%rc_2                 = pdf_params_chnk%rc_2               (:,nzt_clubb:1:-1)
        pdf_params_chnk%cloud_frac_1         = pdf_params_chnk%cloud_frac_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk%cloud_frac_2         = pdf_params_chnk%cloud_frac_2       (:,nzt_clubb:1:-1)
        pdf_params_chnk%chi_1                = pdf_params_chnk%chi_1              (:,nzt_clubb:1:-1)
        pdf_params_chnk%chi_2                = pdf_params_chnk%chi_2              (:,nzt_clubb:1:-1)
        pdf_params_chnk%stdev_chi_1          = pdf_params_chnk%stdev_chi_1        (:,nzt_clubb:1:-1)
        pdf_params_chnk%stdev_chi_2          = pdf_params_chnk%stdev_chi_2        (:,nzt_clubb:1:-1)
        pdf_params_chnk%crt_1                = pdf_params_chnk%crt_1              (:,nzt_clubb:1:-1)
        pdf_params_chnk%crt_2                = pdf_params_chnk%crt_2              (:,nzt_clubb:1:-1)
        pdf_params_chnk%cthl_1               = pdf_params_chnk%cthl_1             (:,nzt_clubb:1:-1)
        pdf_params_chnk%cthl_2               = pdf_params_chnk%cthl_2             (:,nzt_clubb:1:-1)
        pdf_params_chnk%ice_supersat_frac_1  = pdf_params_chnk%ice_supersat_frac_1(:,nzt_clubb:1:-1)
        pdf_params_chnk%ice_supersat_frac_2  = pdf_params_chnk%ice_supersat_frac_2(:,nzt_clubb:1:-1)
        pdf_params_chnk%corr_chi_eta_1       = pdf_params_chnk%corr_chi_eta_1     (:,nzt_clubb:1:-1)
        pdf_params_chnk%corr_chi_eta_2       = pdf_params_chnk%corr_chi_eta_2     (:,nzt_clubb:1:-1)
        pdf_params_chnk%corr_w_chi_1         = pdf_params_chnk%corr_w_chi_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk%corr_w_chi_2         = pdf_params_chnk%corr_w_chi_2       (:,nzt_clubb:1:-1)


        call cleanup_grid_api( gr )

        ! we are in ascending mode, need to recalculate gr in ascending mode
        call setup_grid_api( nzm_clubb, ncol, sfc_elevation, l_implemented,    & ! intent(in)
                             clubb_l_ascending_grid, grid_type,                & ! intent(in)
                             deltaz, zi_g(:,1), zi_g(:,nzm_clubb),             & ! intent(in)
                             zi_g, zt_g,                                       & ! intent(in)
                             gr, err_info )                                      ! intent(inout)

        if ( any(err_info%err_code == clubb_fatal_error) ) then
           errmsg = 'clubb1_run: '//err_info%err_header_global//NEW_LINE('a')// &
                       'in CLUBB setup_grid'
           errflg = 1
           goto 100
        end if

      end if

      ! CALL MF DIAGNOSTIC PLUMES
      if (do_clubb_mf) then
        call compute_mf_plumes()
        if (errflg /= 0) goto 100
      end if

      !  Advance CLUBB CORE one timestep in the future

      ! These updates are required because the pbuf variables are dimensioned with pcols, when
      ! we only need ncol. This requires us to slice the arrays when inputting to advance_clubb_core_api,
      ! which happens on the CPU, so we need the CPU version of these to be correct.
      ! REMOVECAM: This will be unnecessary once pbuf is gone and these are dimensioned ncol.
      !$acc update host(  upwp_pbuf, vpwp_pbuf, up2_pbuf, vp2_pbuf, up3_pbuf, vp3_pbuf, wprtp_pbuf, &
      !$acc               wpthlp_pbuf, wp2_pbuf, wp3_pbuf, rtp2_pbuf, rtp3_pbuf, thlp2_pbuf, thlp3_pbuf, &
      !$acc               rtpthlp_pbuf, wpthvp_pbuf, wp2thvp_pbuf, wp2up_pbuf, rtpthvp_pbuf, thlpthvp_pbuf, wp2rtp_pbuf, &
      !$acc               wp2thlp_pbuf, uprcp_pbuf, vprcp_pbuf, rc_coef_zm_pbuf, wp4_pbuf, wpup2_pbuf, wpvp2_pbuf, &
      !$acc               wp2up2_pbuf, wp2vp2_pbuf, ice_supersat_frac_pbuf )

      call advance_clubb_core_api( gr, nzm_clubb, nzt_clubb, ncol, &        ! Inputs
          l_implemented, dtime, fcor, fcor_y, sfc_elevation, &
          hydromet_dim, &
          sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &
          thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
          sclrm_forcing, edsclrm_forcing, wprtp_forcing, &
          wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
          rtpthlp_forcing, wm_zm, wm_zt, &
          wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
          wpsclrp_sfc, wpedsclrp_sfc, &
          upwp_sfc_pert, vpwp_sfc_pert, &
          rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &
          p_in_Pa, rho_zm, rho_zt, exner, &
          rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
          invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &
          hm_metadata%l_mix_rat_hm, &
          rfrzm, &
          wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &
          grid_dx, grid_dy, &
          clubb_params, nu_vert_res_dep, lmin, &
          mixt_frac_max_mag, theta0, ts_nudge, &
          rtm_min, rtm_nudge_max_altitude, &
          clubb_config_flags, &
          stats_metadata, &
          stats_zt(:ncol), stats_zm(:ncol), stats_sfc(:ncol), &                 ! InOuts
          um, vm, upwp_pbuf(:ncol,:), vpwp_pbuf(:ncol,:), &
          up2_pbuf(:ncol,:), vp2_pbuf(:ncol,:), up3_pbuf(:ncol,:), vp3_pbuf(:ncol,:), &
          thlm, rtm, wprtp_pbuf(:ncol,:), wpthlp_pbuf(:ncol,:), &
          wp2_pbuf(:ncol,:), wp3_pbuf(:ncol,:), rtp2_pbuf(:ncol,:), rtp3_pbuf(:ncol,:), &
          thlp2_pbuf(:ncol,:), thlp3_pbuf(:ncol,:), rtpthlp_pbuf(:ncol,:), &
          sclrm, &
          sclrp2, sclrp3, sclrprtp, sclrpthlp, &
          wpsclrp, edsclr, err_info, &
          rcm, cloud_frac, &
          wpthvp_pbuf(:ncol,:), wp2thvp_pbuf(:ncol,:), wp2up_pbuf(:ncol,:), rtpthvp_pbuf(:ncol,:), thlpthvp_pbuf(:ncol,:), &
          sclrpthvp, &
          wp2rtp_pbuf(:ncol,:), wp2thlp_pbuf(:ncol,:), uprcp_pbuf(:ncol,:), &
          vprcp_pbuf(:ncol,:), rc_coef_zm_pbuf(:ncol,:), &
          wp4_pbuf(:ncol,:), wpup2_pbuf(:ncol,:), wpvp2_pbuf(:ncol,:), &
          wp2up2_pbuf(:ncol,:), wp2vp2_pbuf(:ncol,:), ice_supersat_frac_pbuf(:ncol,:), &
          um_pert, vm_pert, upwp_pert, vpwp_pert, &
          pdf_params_chnk, pdf_params_zm_chnk, &
          pdf_implicit_coefs_terms_chnk, &
          khzm, khzt, &                                                          ! Outputs
          qclvar, thlprcp, &
          wprcp, w_up_in_cloud, w_down_in_cloud, &
          cloudy_updraft_frac, cloudy_downdraft_frac, &
          rcm_in_layer, cloud_cover, invrs_tau_zm, &
          Lscale )

      ! The "unslice" copyback step updates the CPU (host) variables, so we need to copy those back to GPU.
      ! REMOVECAM: This will be unnecessary once pbuf is gone and these are dimensioned ncol.
      !$acc update device( upwp_pbuf, vpwp_pbuf, up2_pbuf, vp2_pbuf, up3_pbuf, vp3_pbuf, wprtp_pbuf, &
      !$acc                wpthlp_pbuf, wp2_pbuf, wp3_pbuf, rtp2_pbuf, rtp3_pbuf, thlp2_pbuf, thlp3_pbuf, &
      !$acc                rtpthlp_pbuf, wpthvp_pbuf, wp2thvp_pbuf, wp2up_pbuf, rtpthvp_pbuf, thlpthvp_pbuf, wp2rtp_pbuf, &
      !$acc                wp2thlp_pbuf, uprcp_pbuf, vprcp_pbuf, rc_coef_zm_pbuf, wp4_pbuf, wpup2_pbuf, wpvp2_pbuf, &
      !$acc                wp2up2_pbuf, wp2vp2_pbuf, ice_supersat_frac_pbuf )

      if ( clubb_l_ascending_grid ) then

        ! If running in ascending mode, we flip the arrays before calling advance_clubb_core
        ! so we need to flip them back. This section should flip every array that was flipped
        ! before the advance_clubb_core call.

        thlm_forcing               =               thlm_forcing(:,nzt_clubb:1:-1)
        rtm_forcing                =                rtm_forcing(:,nzt_clubb:1:-1)
        um_forcing                 =                 um_forcing(:,nzt_clubb:1:-1)
        vm_forcing                 =                 vm_forcing(:,nzt_clubb:1:-1)

        if (do_clubb_mf) then
          mf_qtforcup                =                mf_qtforcup(:,nzt_clubb:1:-1)
          mf_thlforcup               =               mf_thlforcup(:,nzt_clubb:1:-1)
          mf_qtforcdn                =                mf_qtforcdn(:,nzt_clubb:1:-1)
          mf_thlforcdn               =               mf_thlforcdn(:,nzt_clubb:1:-1)
        end if

        wm_zt                      =                      wm_zt(:,nzt_clubb:1:-1)
        rho_zt                     =                     rho_zt(:,nzt_clubb:1:-1)
        rho_ds_zt                  =                  rho_ds_zt(:,nzt_clubb:1:-1)
        invrs_rho_ds_zt            =            invrs_rho_ds_zt(:,nzt_clubb:1:-1)
        thv_ds_zt                  =                  thv_ds_zt(:,nzt_clubb:1:-1)
        khzt                       =                       khzt(:,nzt_clubb:1:-1)
        rtm_ref                    =                    rtm_ref(:,nzt_clubb:1:-1)
        thlm_ref                   =                   thlm_ref(:,nzt_clubb:1:-1)
        um_ref                     =                     um_ref(:,nzt_clubb:1:-1)
        vm_ref                     =                     vm_ref(:,nzt_clubb:1:-1)
        ug                         =                         ug(:,nzt_clubb:1:-1)
        vg                         =                         vg(:,nzt_clubb:1:-1)
        p_in_Pa                    =                    p_in_Pa(:,nzt_clubb:1:-1)
        exner                      =                      exner(:,nzt_clubb:1:-1)
        rfrzm                      =                      rfrzm(:,nzt_clubb:1:-1)
        um                         =                         um(:,nzt_clubb:1:-1)
        vm                         =                         vm(:,nzt_clubb:1:-1)
        up3_pbuf                   =                   up3_pbuf(:,nzt_clubb:1:-1)
        vp3_pbuf                   =                   vp3_pbuf(:,nzt_clubb:1:-1)
        wp3_pbuf                   =                   wp3_pbuf(:,nzt_clubb:1:-1)
        rtp3_pbuf                  =                  rtp3_pbuf(:,nzt_clubb:1:-1)
        thlp3_pbuf                 =                 thlp3_pbuf(:,nzt_clubb:1:-1)
        rcm                        =                        rcm(:,nzt_clubb:1:-1)
        cloud_frac                 =                 cloud_frac(:,nzt_clubb:1:-1)
        wpup2_pbuf                 =                 wpup2_pbuf(:,nzt_clubb:1:-1)
        wpvp2_pbuf                 =                 wpvp2_pbuf(:,nzt_clubb:1:-1)
        wp2rtp_pbuf                =                wp2rtp_pbuf(:,nzt_clubb:1:-1)
        wp2thlp_pbuf               =               wp2thlp_pbuf(:,nzt_clubb:1:-1)
        qclvar                     =                     qclvar(:,nzt_clubb:1:-1)
        cloud_cover                =                cloud_cover(:,nzt_clubb:1:-1)
        w_up_in_cloud              =              w_up_in_cloud(:,nzt_clubb:1:-1)
        w_down_in_cloud            =            w_down_in_cloud(:,nzt_clubb:1:-1)
        cloudy_updraft_frac        =        cloudy_updraft_frac(:,nzt_clubb:1:-1)
        cloudy_downdraft_frac      =      cloudy_downdraft_frac(:,nzt_clubb:1:-1)
        rcm_in_layer               =               rcm_in_layer(:,nzt_clubb:1:-1)
        ice_supersat_frac_pbuf     =     ice_supersat_frac_pbuf(:,nzt_clubb:1:-1)
        um_pert                    =                    um_pert(:,nzt_clubb:1:-1)
        vm_pert                    =                    vm_pert(:,nzt_clubb:1:-1)
        wp2thvp_pbuf               =               wp2thvp_pbuf(:,nzt_clubb:1:-1)
        wp2up_pbuf                 =                 wp2up_pbuf(:,nzt_clubb:1:-1)
        rtm                        =                        rtm(:,nzt_clubb:1:-1)
        thlm                       =                       thlm(:,nzt_clubb:1:-1)
        Lscale                     =                     Lscale(:,nzt_clubb:1:-1)

        ! flipping grid definition
        zt_g                       =                       zt_g(:,nzt_clubb:1:-1)

        if (do_clubb_mf) then
          dz_g                     =                       dz_g(:,nzt_clubb:1:-1)
          invrs_dz_g               =                 invrs_dz_g(:,nzt_clubb:1:-1)
          th_zt                    =                      th_zt(:,nzt_clubb:1:-1)
          qv_zt                    =                      qv_zt(:,nzt_clubb:1:-1)
          qc_zt                    =                      qc_zt(:,nzt_clubb:1:-1)
          kappa_zt                 =                   kappa_zt(:,nzt_clubb:1:-1)
          invrs_exner_zt           =             invrs_exner_zt(:,nzt_clubb:1:-1)
        end if

        wprtp_forcing              =              wprtp_forcing(:,nzm_clubb:1:-1)
        wpthlp_forcing             =             wpthlp_forcing(:,nzm_clubb:1:-1)
        rtp2_forcing               =               rtp2_forcing(:,nzm_clubb:1:-1)
        thlp2_forcing              =              thlp2_forcing(:,nzm_clubb:1:-1)
        rtpthlp_forcing            =            rtpthlp_forcing(:,nzm_clubb:1:-1)
        wm_zm                      =                      wm_zm(:,nzm_clubb:1:-1)
        rho_zm                     =                     rho_zm(:,nzm_clubb:1:-1)
        rho_ds_zm                  =                  rho_ds_zm(:,nzm_clubb:1:-1)
        invrs_rho_ds_zm            =            invrs_rho_ds_zm(:,nzm_clubb:1:-1)
        thv_ds_zm                  =                  thv_ds_zm(:,nzm_clubb:1:-1)
        upwp_pbuf                  =                  upwp_pbuf(:,nzm_clubb:1:-1)
        vpwp_pbuf                  =                  vpwp_pbuf(:,nzm_clubb:1:-1)
        up2_pbuf                   =                   up2_pbuf(:,nzm_clubb:1:-1)
        vp2_pbuf                   =                   vp2_pbuf(:,nzm_clubb:1:-1)
        wprtp_pbuf                 =                 wprtp_pbuf(:,nzm_clubb:1:-1)
        wpthlp_pbuf                =                wpthlp_pbuf(:,nzm_clubb:1:-1)
        wp2_pbuf                   =                   wp2_pbuf(:,nzm_clubb:1:-1)
        rtp2_pbuf                  =                  rtp2_pbuf(:,nzm_clubb:1:-1)
        thlp2_pbuf                 =                 thlp2_pbuf(:,nzm_clubb:1:-1)
        rtpthlp_pbuf               =               rtpthlp_pbuf(:,nzm_clubb:1:-1)
        wpthvp_pbuf                =                wpthvp_pbuf(:,nzm_clubb:1:-1)
        rtpthvp_pbuf               =               rtpthvp_pbuf(:,nzm_clubb:1:-1)
        thlpthvp_pbuf              =              thlpthvp_pbuf(:,nzm_clubb:1:-1)
        uprcp_pbuf                 =                 uprcp_pbuf(:,nzm_clubb:1:-1)
        vprcp_pbuf                 =                 vprcp_pbuf(:,nzm_clubb:1:-1)
        rc_coef_zm_pbuf            =            rc_coef_zm_pbuf(:,nzm_clubb:1:-1)
        wp4_pbuf                   =                   wp4_pbuf(:,nzm_clubb:1:-1)
        wp2up2_pbuf                =                wp2up2_pbuf(:,nzm_clubb:1:-1)
        wp2vp2_pbuf                =                wp2vp2_pbuf(:,nzm_clubb:1:-1)
        upwp_pert                  =                  upwp_pert(:,nzm_clubb:1:-1)
        vpwp_pert                  =                  vpwp_pert(:,nzm_clubb:1:-1)
        khzm                       =                       khzm(:,nzm_clubb:1:-1)
        thlprcp                    =                    thlprcp(:,nzm_clubb:1:-1)
        wprcp                      =                      wprcp(:,nzm_clubb:1:-1)
        invrs_tau_zm               =               invrs_tau_zm(:,nzm_clubb:1:-1)

        ! flipping grid definition
        zi_g                       =                       zi_g(:,nzm_clubb:1:-1)

        if (do_clubb_mf) then
          p_in_Pa_zm               =                 p_in_Pa_zm(:,nzm_clubb:1:-1)
          tke_zm                   =                     tke_zm(:,nzm_clubb:1:-1)
        end if

        if ( edsclr_dim > 0 ) then
          edsclr           =           edsclr(:,nzt_clubb:1:-1,:)
          edsclrm_forcing  =  edsclrm_forcing(:,nzt_clubb:1:-1,:)
        end if

        if ( sclr_dim > 0 ) then

          sclrm_forcing   =   sclrm_forcing(:,nzt_clubb:1:-1,:)
          sclrm           =           sclrm(:,nzt_clubb:1:-1,:)
          sclrp3          =          sclrp3(:,nzt_clubb:1:-1,:)

          sclrp2          =          sclrp2(:,nzm_clubb:1:-1,:)
          sclrprtp        =        sclrprtp(:,nzm_clubb:1:-1,:)
          sclrpthlp       =       sclrpthlp(:,nzm_clubb:1:-1,:)
          wpsclrp         =         wpsclrp(:,nzm_clubb:1:-1,:)
          sclrpthvp       =       sclrpthvp(:,nzm_clubb:1:-1,:)
        end if

        ! These are flipped, ensuring these are stored in descending mode, regardless of clubb_l_ascending_grid
        ! only because these are need to be stored for restarts
        if ( clubb_config_flags%l_call_pdf_closure_twice ) then
          pdf_params_zm_chnk%w_1        = pdf_params_zm_chnk%w_1       (:,nzm_clubb:1:-1)
          pdf_params_zm_chnk%w_2        = pdf_params_zm_chnk%w_2       (:,nzm_clubb:1:-1)
          pdf_params_zm_chnk%varnce_w_1 = pdf_params_zm_chnk%varnce_w_1(:,nzm_clubb:1:-1)
          pdf_params_zm_chnk%varnce_w_2 = pdf_params_zm_chnk%varnce_w_2(:,nzm_clubb:1:-1)
          pdf_params_zm_chnk%mixt_frac  = pdf_params_zm_chnk%mixt_frac (:,nzm_clubb:1:-1)
        end if

        ! These are flipped, ensuring these are stored in descending mode, regardless of clubb_l_ascending_grid
        ! only for pdfp_rtp2_output calc
        pdf_params_chnk%mixt_frac    = pdf_params_chnk%mixt_frac  (:,nzt_clubb:1:-1)
        pdf_params_chnk%rt_1         = pdf_params_chnk%rt_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk%rt_2         = pdf_params_chnk%rt_2       (:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_rt_1  = pdf_params_chnk%varnce_rt_1(:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_rt_2  = pdf_params_chnk%varnce_rt_2(:,nzt_clubb:1:-1)

        ! These are flipped, ensuring these are stored in descending mode, regardless of clubb_l_ascending_grid
        ! only for update_xp2_mc_api call
        pdf_params_chnk%w_1          = pdf_params_chnk%w_1         (:,nzt_clubb:1:-1)
        pdf_params_chnk%w_2          = pdf_params_chnk%w_2         (:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_w_1   = pdf_params_chnk%varnce_w_1  (:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_w_2   = pdf_params_chnk%varnce_w_2  (:,nzt_clubb:1:-1)
        pdf_params_chnk%thl_1        = pdf_params_chnk%thl_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk%thl_2        = pdf_params_chnk%thl_2       (:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_thl_1 = pdf_params_chnk%varnce_thl_1(:,nzt_clubb:1:-1)
        pdf_params_chnk%varnce_thl_2 = pdf_params_chnk%varnce_thl_2(:,nzt_clubb:1:-1)

        ! These are flipped for silhs, which uses a cam grid
        pdf_params_chnk%rc_1                 = pdf_params_chnk%rc_1               (:,nzt_clubb:1:-1)
        pdf_params_chnk%rc_2                 = pdf_params_chnk%rc_2               (:,nzt_clubb:1:-1)
        pdf_params_chnk%cloud_frac_1         = pdf_params_chnk%cloud_frac_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk%cloud_frac_2         = pdf_params_chnk%cloud_frac_2       (:,nzt_clubb:1:-1)
        pdf_params_chnk%chi_1                = pdf_params_chnk%chi_1              (:,nzt_clubb:1:-1)
        pdf_params_chnk%chi_2                = pdf_params_chnk%chi_2              (:,nzt_clubb:1:-1)
        pdf_params_chnk%stdev_chi_1          = pdf_params_chnk%stdev_chi_1        (:,nzt_clubb:1:-1)
        pdf_params_chnk%stdev_chi_2          = pdf_params_chnk%stdev_chi_2        (:,nzt_clubb:1:-1)
        pdf_params_chnk%crt_1                = pdf_params_chnk%crt_1              (:,nzt_clubb:1:-1)
        pdf_params_chnk%crt_2                = pdf_params_chnk%crt_2              (:,nzt_clubb:1:-1)
        pdf_params_chnk%cthl_1               = pdf_params_chnk%cthl_1             (:,nzt_clubb:1:-1)
        pdf_params_chnk%cthl_2               = pdf_params_chnk%cthl_2             (:,nzt_clubb:1:-1)
        pdf_params_chnk%ice_supersat_frac_1  = pdf_params_chnk%ice_supersat_frac_1(:,nzt_clubb:1:-1)
        pdf_params_chnk%ice_supersat_frac_2  = pdf_params_chnk%ice_supersat_frac_2(:,nzt_clubb:1:-1)
        pdf_params_chnk%corr_chi_eta_1       = pdf_params_chnk%corr_chi_eta_1     (:,nzt_clubb:1:-1)
        pdf_params_chnk%corr_chi_eta_2       = pdf_params_chnk%corr_chi_eta_2     (:,nzt_clubb:1:-1)
        pdf_params_chnk%corr_w_chi_1         = pdf_params_chnk%corr_w_chi_1       (:,nzt_clubb:1:-1)
        pdf_params_chnk%corr_w_chi_2         = pdf_params_chnk%corr_w_chi_2       (:,nzt_clubb:1:-1)

        call cleanup_grid_api( gr )

        ! recalculate descending grid
        call setup_grid_api( nzm_clubb, ncol, sfc_elevation, l_implemented,   & ! intent(in)
                             .false., grid_type,                              & ! intent(in)
                             deltaz, zi_g(:,nzm_clubb), zi_g(:,1),            & ! intent(in)
                             zi_g, zt_g,                                      & ! intent(in)
                             gr, err_info )                                     ! intent(inout)

      end if

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) "Fatal error in CLUBB advance_clubb_core: at timestep ", nstep
        errmsg = 'clubb1_run: '//err_info%err_header_global//NEW_LINE('a')//'Fatal error in CLUBB advance_clubb_core'
        errflg = 1
        goto 100
      end if

      if ( do_rainturb ) then

        do k = 1, nzt_clubb
          do i = 1, ncol
            k_cam = top_lev - 1 + k
            rvm(i,k) = rtm(i,k) - rcm(i,k)
            pre(i,k) = prer_evap_pbuf(i,k_cam)
          end do
        end do

        call update_xp2_mc_api( gr, nzm_clubb, nzt_clubb, ncol, dtime, cloud_frac, &
                                rcm(:ncol,:), rvm, thlm(:ncol,:), wm_zt, &
                                exner, pre, pdf_params_chnk, &
                                rtp2_mc, thlp2_mc, &
                                wprtp_mc, wpthlp_mc, &
                                rtpthlp_mc)

        do k = 1, nzm_clubb
          do i = 1, ncol
            dum1 = (1._kind_phys - landfrac(i))

            ! update turbulent moments based on rain evaporation
            rtp2_pbuf(i,k)   = rtp2_pbuf(i,k)   + clubb_rnevap_effic * dum1 * rtp2_mc(i,k)   * dtime
            thlp2_pbuf(i,k)  = thlp2_pbuf(i,k)  + clubb_rnevap_effic * dum1 * thlp2_mc(i,k)  * dtime
            wprtp_pbuf(i,k)  = wprtp_pbuf(i,k)  + clubb_rnevap_effic * dum1 * wprtp_mc(i,k)  * dtime
            wpthlp_pbuf(i,k) = wpthlp_pbuf(i,k) + clubb_rnevap_effic * dum1 * wpthlp_mc(i,k) * dtime
          end do
        end do

      end if

      if (do_cldcool) then

        thlp2_rad(:,:) = 0._kind_phys

        do k = 1, nzt_clubb
          do i = 1, ncol
            k_cam = top_lev - 1 + k
            qrl_clubb(i,k) = qrl_pbuf(i,k_cam) / ( cpairv(i,k_cam) * state_pdeldry(i,k_cam) )
          end do
        end do

        call calculate_thlp2_rad_api( ncol, nzm_clubb, nzt_clubb, gr, &
                                      rcm(:ncol,:), thlprcp, qrl_clubb, clubb_params, &
                                      thlp2_rad )

        do k = 1, nzm_clubb
          do i = 1, ncol
            thlp2_pbuf(i,k) = max( thl_tol**2, thlp2_pbuf(i,k) + thlp2_rad(i,k) * dtime )
          end do
        end do

      end if

      !  Check to see if stats should be output, here stats are read into
      !  output arrays to make them conformable to CAM output
      if (stats_metadata%l_stats) then
        do i = 1, ncol
          call stats_end_timestep_clubb(i, stats_zt(i), stats_zm(i), stats_rad_zt(i), stats_rad_zm(i), stats_sfc(i), &
                                        out_zt, out_zm, out_radzt, out_radzm, out_sfc, &
                                        stats_metadata, clubb_l_ascending_grid, pver, pverp, top_lev, &
                                        errmsg, errflg )
          if ( errflg /= 0 ) goto 100
        end do
      end if

    end do  ! end time loop
    !----------------------------------------- END substepping loop -----------------------------------------

    if ( clubb_l_ascending_grid .and. do_clubb_mf ) then

      ! integrate_mf leaves its output in ascending order, but the diagnostic
      ! mapping further down (mf_*_output(i,k) = mf_*(i,k_clubb), with
      ! k_clubb = k + 1 - top_lev) assumes the descending host-model convention.
      ! Restore descending order here so that everything downstream of the
      ! substepping loop can assume the host convention, as the rest of this
      ! module does.
      !
      ! This is done once per macmic step rather than inside nadv
      ! flip-back block, because only the final substep's values reach the
      ! output mapping. Note that some of those expressions pair an MF flux with
      ! rho_zm, which is already restored to descending, so flipping here also
      ! fixes that mixed-orientation product.

      ! MF plume ensemble means on momentum levels
      mf_dry_a     =     mf_dry_a(:,nzm_clubb:1:-1)
      mf_moist_a   =   mf_moist_a(:,nzm_clubb:1:-1)
      mf_dry_w     =     mf_dry_w(:,nzm_clubb:1:-1)
      mf_moist_w   =   mf_moist_w(:,nzm_clubb:1:-1)
      mf_dry_qt    =    mf_dry_qt(:,nzm_clubb:1:-1)
      mf_moist_qt  =  mf_moist_qt(:,nzm_clubb:1:-1)
      mf_dry_thl   =   mf_dry_thl(:,nzm_clubb:1:-1)
      mf_moist_thl = mf_moist_thl(:,nzm_clubb:1:-1)
      mf_dry_u     =     mf_dry_u(:,nzm_clubb:1:-1)
      mf_moist_u   =   mf_moist_u(:,nzm_clubb:1:-1)
      mf_dry_v     =     mf_dry_v(:,nzm_clubb:1:-1)
      mf_moist_v   =   mf_moist_v(:,nzm_clubb:1:-1)
      mf_moist_qc  =  mf_moist_qc(:,nzm_clubb:1:-1)

      ! Sub-plume moment sums on momentum levels
      s_ae         =         s_ae(:,nzm_clubb:1:-1)
      s_aw         =         s_aw(:,nzm_clubb:1:-1)
      s_aww        =        s_aww(:,nzm_clubb:1:-1)
      s_awthl      =      s_awthl(:,nzm_clubb:1:-1)
      s_awthlup    =    s_awthlup(:,nzm_clubb:1:-1)
      s_awthldn    =    s_awthldn(:,nzm_clubb:1:-1)
      s_awqt       =       s_awqt(:,nzm_clubb:1:-1)
      s_awqtup     =     s_awqtup(:,nzm_clubb:1:-1)
      s_awqtdn     =     s_awqtdn(:,nzm_clubb:1:-1)
      s_awu        =        s_awu(:,nzm_clubb:1:-1)
      s_awuup      =      s_awuup(:,nzm_clubb:1:-1)
      s_awudn      =      s_awudn(:,nzm_clubb:1:-1)
      s_awv        =        s_awv(:,nzm_clubb:1:-1)
      s_awvup      =      s_awvup(:,nzm_clubb:1:-1)
      s_awvdn      =      s_awvdn(:,nzm_clubb:1:-1)

      ! deep-hookup ensemble mf/ent/det profiles on momentum levels: the
      ! nadv accumulators are flipped (the raw per-subcycle arrays hold only
      ! the last subcycle and are not used downstream).  mf_kctop_nadv needs
      ! no flip -- kctop is an orientation-independent offset from the surface.
      s_mfup_nadv  =  s_mfup_nadv(:,nzm_clubb:1:-1)
      s_entup_nadv = s_entup_nadv(:,nzm_clubb:1:-1)
      s_detup_nadv = s_detup_nadv(:,nzm_clubb:1:-1)
      s_mfdn_nadv  =  s_mfdn_nadv(:,nzm_clubb:1:-1)
      s_entdn_nadv = s_entdn_nadv(:,nzm_clubb:1:-1)
      s_detdn_nadv = s_detdn_nadv(:,nzm_clubb:1:-1)
      s_awwup_nadv = s_awwup_nadv(:,nzm_clubb:1:-1)
      s_awup_nadv  =  s_awup_nadv(:,nzm_clubb:1:-1)

      ! MF fluxes on momentum levels
      mf_thlflx    =    mf_thlflx(:,nzm_clubb:1:-1)
      mf_thlflxup  =  mf_thlflxup(:,nzm_clubb:1:-1)
      mf_thlflxdn  =  mf_thlflxdn(:,nzm_clubb:1:-1)
      mf_qtflx     =     mf_qtflx(:,nzm_clubb:1:-1)
      mf_qtflxup   =   mf_qtflxup(:,nzm_clubb:1:-1)
      mf_qtflxdn   =   mf_qtflxdn(:,nzm_clubb:1:-1)
      mf_uflx      =      mf_uflx(:,nzm_clubb:1:-1)
      mf_uflxup    =    mf_uflxup(:,nzm_clubb:1:-1)
      mf_uflxdn    =    mf_uflxdn(:,nzm_clubb:1:-1)
      mf_vflx      =      mf_vflx(:,nzm_clubb:1:-1)
      mf_vflxup    =    mf_vflxup(:,nzm_clubb:1:-1)
      mf_vflxdn    =    mf_vflxdn(:,nzm_clubb:1:-1)
      mf_thvflx    =    mf_thvflx(:,nzm_clubb:1:-1)

      mf_precc     =     mf_precc(:,nzm_clubb:1:-1)

      ! Plume microphysics source terms on thermodynamic levels, needed for
      ! edmf_sqtup / edmf_sqtdn.
      mf_sqtup     =     mf_sqtup(:,nzt_clubb:1:-1)
      mf_sqtdn     =     mf_sqtdn(:,nzt_clubb:1:-1)

      ! deep-hookup plume autoconversion / rain evaporation on thermo
      ! levels (nadv accumulators, as above)
      mf_sqtac_nadv = mf_sqtac_nadv(:,nzt_clubb:1:-1)
      mf_sqtev_nadv = mf_sqtev_nadv(:,nzt_clubb:1:-1)

      ! NOTE: mf_rcm needs no flip here -- it is overwritten from mf_rcm_nadv
      ! below, and that accumulator is flipped at the end of this block.

      ! Per-plume profiles on momentum levels
      mf_upa       =       mf_upa(:,nzm_clubb:1:-1,:)
      mf_upw       =       mf_upw(:,nzm_clubb:1:-1,:)
      mf_upmf      =      mf_upmf(:,nzm_clubb:1:-1,:)
      mf_upqt      =      mf_upqt(:,nzm_clubb:1:-1,:)
      mf_upthl     =     mf_upthl(:,nzm_clubb:1:-1,:)
      mf_upthv     =     mf_upthv(:,nzm_clubb:1:-1,:)
      mf_upth      =      mf_upth(:,nzm_clubb:1:-1,:)
      mf_upqc      =      mf_upqc(:,nzm_clubb:1:-1,:)
      mf_upbuoy    =    mf_upbuoy(:,nzm_clubb:1:-1,:)
      mf_upent     =     mf_upent(:,nzm_clubb:1:-1,:)
      mf_updet     =     mf_updet(:,nzm_clubb:1:-1,:)
      mf_dnw       =       mf_dnw(:,nzm_clubb:1:-1,:)
      mf_dnthl     =     mf_dnthl(:,nzm_clubb:1:-1,:)
      mf_dnqt      =      mf_dnqt(:,nzm_clubb:1:-1,:)

      ! Substep accumulators for the cloud fields. These are flipped so that the
      ! zm2zt_api call below sees descending data against the descending gr that
      ! was restored at the end of the last substep; without this the
      ! interpolation weights would not match the data.
      mf_qc_nadv        =        mf_qc_nadv(:,nzm_clubb:1:-1)
      mf_rcm_nadv       =       mf_rcm_nadv(:,nzm_clubb:1:-1)
      mf_cloudfrac_nadv = mf_cloudfrac_nadv(:,nzm_clubb:1:-1)

    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol
        qclvar(i,k)        = min( 1._kind_phys, qclvar(i,k) ) ! We should move this clipping inside clubb
      end do
    end do

    if (do_clubb_mf .and. do_clubb_mf_addtke) then
      Lscale_zm = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, Lscale )
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm_clubb
        do i = 1, ncol
          k_cam = top_lev - 1 + k
          wp2_pbuf(i,k) = wp2_pbuf(i,k) + s_aww(i,k)
          khzm_pbuf(i,k_cam) = khzm(i,k) + &
               clubb_params(i,ic_K) * Lscale_zm(i,k) * sqrt(0.5_kind_phys * s_aww(i,k))
        end do
      end do
    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm_clubb
        do i = 1, ncol
          k_cam = top_lev - 1 + k
          khzm_pbuf(i,k_cam)         = khzm(i,k)
        end do
      end do
    end if

    ! pdf_params_zm_chnk is already persistent across calls, but we
    ! save a pbuf version for restarts
    if ( clubb_config_flags%l_call_pdf_closure_twice ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm_clubb
        do i = 1, ncol
          pdf_zm_w_1_pbuf(i,k)        = pdf_params_zm_chnk%w_1(i,k)
          pdf_zm_w_2_pbuf(i,k)        = pdf_params_zm_chnk%w_2(i,k)
          pdf_zm_varnce_w_1_pbuf(i,k) = pdf_params_zm_chnk%varnce_w_1(i,k)
          pdf_zm_varnce_w_2_pbuf(i,k) = pdf_params_zm_chnk%varnce_w_2(i,k)
          pdf_zm_mixt_frac_pbuf(i,k)  = pdf_params_zm_chnk%mixt_frac(i,k)
        end do
      end do
    end if

    ! Compute static energy using CLUBB's variables
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = top_lev, pver
      do i = 1, ncol
        k_clubb = k + 1 - top_lev
        clubb_s(i,k_clubb) = cpairv(i,k) * thlm(i,k_clubb) / invrs_exner_zt(i,k_clubb) &
                       + latvap * rcm(i,k_clubb) &
                       + gravit * state_zm(i,k) + state_phis(i)
      end do
    end do

    ! Section below is concentrated on energy fixing for conservation.
    !   because CLUBB and CAM's thermodynamic variables are different.

    ! Initialize clubbtop_pbuf to top_lev, for finding the highlest level CLUBB is
    !  active for informing where to apply the energy fixer.
    !$acc parallel loop gang vector default(present)
    do i = 1, ncol
      clubbtop_pbuf(i) = top_lev
      k_clubb     = clubbtop_pbuf(i) + 1 - top_lev
      do while ((rtp2_pbuf(i,k_clubb) <= 1.e-15_kind_phys .and. rcm(i,k_clubb)  ==  0._kind_phys) .and. clubbtop_pbuf(i) <  pver)
        clubbtop_pbuf(i) = clubbtop_pbuf(i) + 1
        k_clubb          = clubbtop_pbuf(i) + 1 - top_lev
      end do
    end do

    ! process clubb-mf accumulated arrays
    if (do_clubb_mf) then
      ! average over nadv
      mf_L0_nadv   = mf_L0_nadv/REAL(nadv,kind_phys)
      mf_ztop_nadv = mf_ztop_nadv/REAL(nadv,kind_phys)
      mf_ztopm1_nadv = mf_ztopm1_nadv/REAL(nadv,kind_phys)
      mf_cape_nadv = mf_cape_nadv/REAL(nadv,kind_phys)
      mf_ddcp_nadv = mf_ddcp_nadv/REAL(nadv,kind_phys)
      mf_cbm1_nadv = mf_cbm1_nadv/REAL(nadv,kind_phys)
      mf_freq_nadv = mf_freq_nadv/REAL(nadv,kind_phys)

      mf_qc(:ncol,:nzm_clubb)        = mf_qc_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      mf_rcm(:ncol,:nzm_clubb)       = mf_rcm_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      mf_cloudfrac(:ncol,:nzm_clubb) = mf_cloudfrac_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      prec_dp_pbuf(:ncol)            = mf_precc_nadv(:ncol)/REAL(nadv,kind_phys)
      snow_dp_pbuf(:ncol)            = mf_snow_nadv(:ncol)/REAL(nadv,kind_phys)

      mf_thlforcup_nadv(:ncol,:nzt_clubb) = mf_thlforcup_nadv(:ncol,:nzt_clubb)/REAL(nadv,kind_phys)
      mf_qtforcup_nadv(:ncol,:nzt_clubb)  = mf_qtforcup_nadv(:ncol,:nzt_clubb)/REAL(nadv,kind_phys)
      mf_thlforcdn_nadv(:ncol,:nzt_clubb) = mf_thlforcdn_nadv(:ncol,:nzt_clubb)/REAL(nadv,kind_phys)
      mf_qtforcdn_nadv(:ncol,:nzt_clubb)  = mf_qtforcdn_nadv(:ncol,:nzt_clubb)/REAL(nadv,kind_phys)
      mf_thlforc_nadv(:ncol,:nzt_clubb)   = mf_thlforc_nadv(:ncol,:nzt_clubb)/REAL(nadv,kind_phys)
      mf_qtforc_nadv(:ncol,:nzt_clubb)    = mf_qtforc_nadv(:ncol,:nzt_clubb)/REAL(nadv,kind_phys)
      mf_ent_nadv(:ncol,:nzm_clubb)       = mf_ent_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)

      s_mfup_nadv(:ncol,:nzm_clubb)   = s_mfup_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      s_entup_nadv(:ncol,:nzm_clubb)  = s_entup_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      s_detup_nadv(:ncol,:nzm_clubb)  = s_detup_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      s_mfdn_nadv(:ncol,:nzm_clubb)   = s_mfdn_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      s_entdn_nadv(:ncol,:nzm_clubb)  = s_entdn_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      s_detdn_nadv(:ncol,:nzm_clubb)  = s_detdn_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      s_awwup_nadv(:ncol,:nzm_clubb)  = s_awwup_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      s_awup_nadv(:ncol,:nzm_clubb)   = s_awup_nadv(:ncol,:nzm_clubb)/REAL(nadv,kind_phys)
      mf_sqtac_nadv(:ncol,:nzt_clubb) = mf_sqtac_nadv(:ncol,:nzt_clubb)/REAL(nadv,kind_phys)
      mf_sqtev_nadv(:ncol,:nzt_clubb) = mf_sqtev_nadv(:ncol,:nzt_clubb)/REAL(nadv,kind_phys)

      ! accumulate in buffer
      ztopm1_macmic(:ncol,:) = ztopm1_macmic(:ncol,:) + mf_ztopm1_nadv(:ncol,:)
      ddcp_macmic(:ncol,:) = ddcp_macmic(:ncol,:) + mf_ddcp_nadv(:ncol,:)
      cbm1_macmic(:ncol) = cbm1_macmic(:ncol) + mf_cbm1_nadv(:ncol)

      ! accumulate the nadv-averaged deep-hookup variables in buffer
      s_mfup_macmic(:ncol,:)   = s_mfup_macmic(:ncol,:)   + s_mfup_nadv(:ncol,:nzm_clubb)
      s_entup_macmic(:ncol,:)  = s_entup_macmic(:ncol,:)  + s_entup_nadv(:ncol,:nzm_clubb)
      s_detup_macmic(:ncol,:)  = s_detup_macmic(:ncol,:)  + s_detup_nadv(:ncol,:nzm_clubb)
      s_mfdn_macmic(:ncol,:)   = s_mfdn_macmic(:ncol,:)   + s_mfdn_nadv(:ncol,:nzm_clubb)
      s_entdn_macmic(:ncol,:)  = s_entdn_macmic(:ncol,:)  + s_entdn_nadv(:ncol,:nzm_clubb)
      s_detdn_macmic(:ncol,:)  = s_detdn_macmic(:ncol,:)  + s_detdn_nadv(:ncol,:nzm_clubb)
      s_awwup_macmic(:ncol,:)  = s_awwup_macmic(:ncol,:)  + s_awwup_nadv(:ncol,:nzm_clubb)
      s_awup_macmic(:ncol,:)   = s_awup_macmic(:ncol,:)   + s_awup_nadv(:ncol,:nzm_clubb)
      mf_sqtac_macmic(:ncol,:) = mf_sqtac_macmic(:ncol,:) + mf_sqtac_nadv(:ncol,:nzt_clubb)
      mf_sqtev_macmic(:ncol,:) = mf_sqtev_macmic(:ncol,:) + mf_sqtev_nadv(:ncol,:nzt_clubb)
      mf_kctop_macmic(:ncol)   = max(mf_kctop_macmic(:ncol), mf_kctop_nadv(:ncol))

      if (macmic_it == cld_macmic_num_steps) then

        cbm1(:ncol) = cbm1_macmic(:ncol)/REAL(cld_macmic_num_steps,kind_phys)

        if (clubb_mf_up_ndt == 1) then
          ztopma(:ncol,:) = ztopm1_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        else
          ztopmn(2:clubb_mf_up_ndt,:ncol,:) = ztopmn(1:clubb_mf_up_ndt-1,:ncol,:)
          ztopmn(1,:ncol,:) = ztopm1_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
          ztopma(:ncol,:) = 0._kind_phys
          do tt=1,clubb_mf_up_ndt
            ztopma(:ncol,:) = ztopma(:ncol,:) + ztopmn(tt,:ncol,:)
          end do
          ztopma(:ncol,:) = ztopma(:ncol,:)/REAL(clubb_mf_up_ndt,kind_phys)
        end if

        if (clubb_mf_cp_ndt == 1) then
          ddcp(:ncol,:) = ddcp_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        else
          ddcpmn(2:clubb_mf_cp_ndt,:ncol,:) = ddcpmn(1:clubb_mf_cp_ndt-1,:ncol,:)
          ddcpmn(1,:ncol,:) = ddcp_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
          ddcp(:ncol,:) = 0._kind_phys
          do tt=1,clubb_mf_cp_ndt
            ddcp(:ncol,:) = ddcp(:ncol,:) + ddcpmn(tt,:ncol,:)
          end do
          ddcp(:ncol,:) = ddcp(:ncol,:)/REAL(clubb_mf_cp_ndt,kind_phys)
        end if ! if (clubb_mf_cp_ndt == 1) then

        ddcp(:ncol,:) = clubb_mf_ddalph*ddcp(:ncol,:)

      end if ! if (macmic_it == cld_macmic_num_steps)

      ! Need moist_qc and cloudfrac on thermo grid for output
      mf_qc_zt(:,:) = 0._kind_phys
      mf_cloudfrac_zt(:,:) = 0._kind_phys
      mf_qc_zt = zm2zt_api( nzm_clubb, nzt_clubb, ncol, gr, mf_qc)
      mf_cloudfrac_zt = zm2zt_api( nzm_clubb, nzt_clubb, ncol, gr, mf_cloudfrac)

      ! ------------------------------------------------------------------- !
      ! Deep-convection hookup:                                             !
      ! populate the deep-convection pbuf arrays from the MF plume ensemble !
      ! so that aero_convproc (aerosol transport/activation/scavenging,     !
      ! convproc_do_aer=.true.) and zm_conv_convtran_run via                !
      ! convect_deep_tend_2 (tracer transport) operate on the plumes.       !
      ! both consumers run later in tphysac this same timestep.             !
      !                                                                     !
      ! the pbuf arrays are only populated on the LAST macmic iteration,    !
      ! from the deep-hookup variables averaged over all nadv subcycles and !
      ! macmic iterations.  a column is gathered (ideep) if it had an       !
      ! active plume ensemble at ANY point in the sub-cycles, i.e. if the   !
      ! highest plume top (kctop, a max over the sub-cycles) cleared the    !
      ! first momentum interface.                                           !
      ! ------------------------------------------------------------------- !
      if (macmic_it == cld_macmic_num_steps) then

        ! average the macmic accumulators over the macmic iterations, back
        ! into the working arrays the population code below reads
        s_mfup(:ncol,:)   = s_mfup_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        s_entup(:ncol,:)  = s_entup_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        s_detup(:ncol,:)  = s_detup_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        s_mfdn(:ncol,:)   = s_mfdn_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        s_entdn(:ncol,:)  = s_entdn_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        s_detdn(:ncol,:)  = s_detdn_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        mf_sqtac(:ncol,:) = mf_sqtac_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        mf_sqtev(:ncol,:) = mf_sqtev_macmic(:ncol,:)/REAL(cld_macmic_num_steps,kind_phys)
        ! plume top: highest kctop over all nadv/macmic sub-cycles (max, no averaging)
        mf_kctop(:ncol)   = mf_kctop_macmic(:ncol)

        mu(:,:)  = 0._kind_phys
        eu(:,:)  = 0._kind_phys
        du(:,:)  = 0._kind_phys
        md(:,:)  = 0._kind_phys
        ed(:,:)  = 0._kind_phys
        dpg(:,:) = 0._kind_phys
        dsubcld(:)   = 0._kind_phys
        wup_gath(:,:) = 0._kind_phys
        jt(:)    = pver
        maxg(:)  = pver
        ideep(:) = 0

        lengath = 0
        do i=1, ncol
          ! mf_kctop counts momentum interfaces from the surface, so the CAM
          ! interface index of the ensemble plume top is pverp-mf_kctop+1; the
          ! deep-convection cloud-top LEVEL is one above the interface.
          ! mf_kctop here is the HIGHEST plume top over the nadv/macmic sub-cycles.
          tmptop = min(max(pverp - int(mf_kctop(i)), top_lev), pver)
          ! gather every column that was active at ANY point in the
          ! nadv/macmic sub-cycles (kctop is the max over the sub-cycles,
          ! so int(kctop) > 1 means the plume top cleared the first
          ! momentum interface at some sub-cycle)
          if (int(mf_kctop(i)) > 1) then
            lengath = lengath + 1
            ideep(lengath) = i
            jt(lengath)    = tmptop
            maxg(lengath)  = pver               ! surface-launched plumes
            dpg(lengath,:pver) = state_pdel(i,:pver)/100._kind_phys
          end if
        end do

        rprddp(:,:)  = 0._kind_phys
        evapcdp(:,:) = 0._kind_phys
        do kcam = top_lev, pver

          k_clubb = kcam + 1 - top_lev

          do i=1, lengath
            icol = ideep(i)
            ! plume MFs are in kg/m2/s; ZM wants hPa/s.  mu/md at CAM level kcam
            ! hold the ensemble mass flux at the interface at the layer's top.
            mu(i,kcam) = s_mfup(icol,k_clubb)*gravit/100._kind_phys
            ! mass-flux-weighted ensemble updraft speed at the same
            ! interface as mu (sum a*w^2 / sum a*w over all nadv and
            ! macmic sub-cycles. deliberately NOT Courant-scaled: aero_convproc's
            ! activation should see the physical plume w, not the
            ! transport-limited effective velocity. bounded for safety.
            wup_gath(i,kcam) = min(50._kind_phys, max(0._kind_phys, &
                 s_awwup_macmic(icol,k_clubb) &
                 /max(s_awup_macmic(icol,k_clubb), 1.e-10_kind_phys)))
            md(i,kcam) = s_mfdn(icol,k_clubb)*gravit/100._kind_phys
            ! s_entup carries the sub-cycle-averaged entrainment MASS
            ! FLUX (fractional entrainment x plume mass flux, kg/m3/s)
            eu(i,kcam) = s_entup(icol,k_clubb)*(gravit/100._kind_phys)*dz_g(icol,k_clubb)/dpg(i,kcam)
          end do
          do i=1, ncol
            rprddp(i,kcam)  = -1._kind_phys*mf_sqtac(i,k_clubb)
            evapcdp(i,kcam) = mf_sqtev(i,k_clubb)
          end do
        end do

        ! convtran/convproc REQUIRE (a) a Courant limit mu*dt/dp <= 1 (ZM enforces
        ! this on its cloud-base closure; violating it lets the flux-form limiters
        ! remove more tracer than a layer holds -> negative gas concentrations ->
        ! rrtmgp set_vmr abort) and (b) EXACT discrete continuity
        ! mu(k) - mu(k+1) = (eu(k)-du(k))*dp(k) with mu(pver+1)=0 (and the
        ! analogous downdraft relation), which the recursion for the in-updraft
        ! concentration assumes.  scale the plume mass fluxes to the Courant
        ! target, then rebuild du/ed from continuity (keeping the plume-derived
        ! eu wherever a non-negative du can absorb it).
        do i=1, lengath
          ! Courant check against DRY pressure thickness: convtran scales eu/du/ed
          ! by dp/dpdry and integrates dry-mixing-ratio constituents (the GHGs)
          ! against dpdry, so the effective Courant number is mu*dt/dpdry.  Target
          ! 0.5 (rather than the ZM-style 1.0) to suppress the residual flux-form
          ! limiter overshoots that drive convtran2 QNEG3 negatives.
          ! the aero_convproc courant number also carries the entrainment
          ! flux term (eudp = eu*dpdry), so include eu*dt in the cap --
          ! otherwise its internal substep count ntsub = 1+int(courant)
          ! is unbounded by the mass-flux cap alone
          tmpcfl = 0._kind_phys
          do kcam=1, pver
            tmpcfl = max(tmpcfl, max(mu(i,kcam), abs(md(i,kcam)))*hdtime &
                                  /(state_pdeldry(ideep(i),kcam)/100._kind_phys) &
                                 + eu(i,kcam)*hdtime)
          end do
          if (tmpcfl > 0.5_kind_phys) then
            mu(i,:pver) = mu(i,:pver)*0.5_kind_phys/tmpcfl
            md(i,:pver) = md(i,:pver)*0.5_kind_phys/tmpcfl
            eu(i,:pver) = eu(i,:pver)*0.5_kind_phys/tmpcfl
          end if
          do kcam=pver, 1, -1
            ! updraft: net = d(mu)/dp across layer kcam (mu below surface = 0)
            if (kcam == pver) then
              tmpnet = mu(i,kcam)/dpg(i,kcam)
            else
              tmpnet = (mu(i,kcam) - mu(i,kcam+1))/dpg(i,kcam)
            end if
            eu(i,kcam) = max(eu(i,kcam), tmpnet, 0._kind_phys)
            du(i,kcam) = eu(i,kcam) - tmpnet
            ! downdraft: md(k+1) = md(k) - ed(k)*dp(k); ed >= 0
            if (kcam < pver) then
              ed(i,kcam) = max((md(i,kcam) - md(i,kcam+1))/dpg(i,kcam), 0._kind_phys)
            else
              ed(i,kcam) = 0._kind_phys
            end if
          end do
        end do

      end if ! macmic_it == cld_macmic_num_steps (deep hookup)

    end if !clubbmf

    !$acc parallel loop gang vector default(present)
    do i = 1, ncol

      se_a = 0._kind_phys
      ke_a = 0._kind_phys
      wv_a = 0._kind_phys
      wl_a = 0._kind_phys

      se_b = 0._kind_phys
      ke_b = 0._kind_phys
      wv_b = 0._kind_phys
      wl_b = 0._kind_phys

      ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
      ! after CLUBB is called.  This is for energy conservation purposes.
      do k = top_lev, pver
        k_clubb     = k + 1 - top_lev
        se_a = se_a + clubb_s(i,k_clubb)*state_pdel(i,k)*rga
        ke_a = ke_a + 0.5_kind_phys*(um(i,k_clubb)**2+vm(i,k_clubb)**2)*state_pdel(i,k)*rga
        wv_a = wv_a + (rtm(i,k_clubb)-rcm(i,k_clubb))*state_pdeldry(i,k)*rga
        wl_a = wl_a + (rcm(i,k_clubb))*state_pdeldry(i,k)*rga
      end do

      ! Based on these integrals, compute the total energy after CLUBB call
      te_a = se_a + ke_a + (latvap+latice) * wv_a + latice * wl_a

      do k = top_lev, pver
        ! Do the same as above, but for before CLUBB was called.
        se_b = se_b + s(i,k)*state_pdel(i,k)*rga
        ke_b = ke_b + 0.5_kind_phys*(u(i,k)**2+v(i,k)**2)*state_pdel(i,k)*rga
        wv_b = wv_b + state_q(i,k,ixq)*state_pdeldry(i,k)*rga
        wl_b = wl_b + state_q(i,k,ixcldliq)*state_pdeldry(i,k)*rga
      end do

      ! Based on these integrals, compute the total energy before CLUBB call
      te_b = se_b + ke_b + (latvap+latice) * wv_b + latice * wl_b

      ! Take into account the surface fluxes of heat and moisture
      !  Use correct qflux from cam_in, not lhf/latvap as was done previously
      te_b = te_b + (shf(i)+cflx(i,1)*(latvap+latice)) * hdtime

      if (do_clubb_mf) then
         ! subtract enthalpy of falling precip from tb
         te_b = te_b - prec_dp_pbuf(i)*1000._kind_phys*latice*hdtime
      end if

      ! Compute the disbalance of total energy, over depth where CLUBB is active
      se_dis(i) = ( te_a - te_b ) / ( state_pint(i,pverp) - state_pint(i,clubbtop_pbuf(i)) )

      eleak(i) = ( te_a - te_b ) * invrs_hdtime

    end do

    ! Fix the total energy coming out of CLUBB so it achieves energy conservation.
    ! Apply this fixer throughout the column evenly, but only at layers where
    ! CLUBB is active.
    !
    ! NOTE: The energy fixer seems to cause the climate to change significantly
    ! when using specified dynamics, so allow this to be turned off via a namelist
    ! variable.
    if (clubb_do_energyfix) then

      !$acc parallel loop gang vector default(present)
      do i = 1, ncol

        do k = clubbtop_pbuf(i), pver
          k_clubb = k + 1 - top_lev
          clubb_s(i,k_clubb) = clubb_s(i,k_clubb) - se_dis(i) * gravit
        end do
        ! convert to units of +ve [K]
        se_dis(i) = -1._kind_phys * se_dis(i) * gravit * invrs_cpairv(i,pver)

      end do

    endif

    100 continue
    !$acc end data      ! closes gr/nu_vert_res_dep region
    200 continue
    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data

    if (errflg /= 0) return

    ! ------------------------------------------------- !
    ! Diagnose relative cloud water variance            !
    ! ------------------------------------------------- !

    if (deep_scheme  ==  'CLUBB_SGS') then
      relvarmax = 2.0_kind_phys
    else
      relvarmax = 10.0_kind_phys
    endif

    do k = 1, pver
      do i = 1, ncol
        relvar_pbuf(i,k) = relvarmax  ! default
      end do
    end do

    if (deep_scheme .ne. 'CLUBB_SGS') then
      do k = top_lev, pver
        do i = 1, ncol
          k_clubb = k + 1 - top_lev
          if ( rcm(i,k_clubb) /= 0 .and. qclvar(i,k_clubb) /= 0 ) then
            relvar_pbuf(i,k) = min( relvarmax, max(0.001_kind_phys, rcm(i,k_clubb)**2 / qclvar(i,k_clubb) ) )
          end if
        end do
      end do
    endif

    !  turbulent kinetic energy
    do k = top_lev, pverp
      do i = 1, ncol
        k_clubb     = k + 1 - top_lev
        tke_pbuf(i,k) = 0.5_kind_phys * ( up2_pbuf(i,k_clubb) + vp2_pbuf(i,k_clubb) + wp2_pbuf(i,k_clubb) )
      enddo
    enddo

    do k = top_lev, pver
      do i = 1, ncol
        k_clubb     = k + 1 - top_lev
        ptend_u(i,k)          = ( um(i,k_clubb) - u(i,k))           * invrs_hdtime ! east-west wind
        ptend_v(i,k)          = ( vm(i,k_clubb) - v(i,k))           * invrs_hdtime ! north-south wind
        ptend_q(i,k,ixq)      = ( rtm(i,k_clubb) - rcm(i,k_clubb) &
                                      -state_q(i,k,ixq) )                     * invrs_hdtime ! water vapor
        ptend_q(i,k,ixcldliq) = ( rcm(i,k_clubb) - state_q(i,k,ixcldliq)) * invrs_hdtime ! Tendency of liquid water
        ptend_s(i,k)          = ( clubb_s(i,k_clubb) - s(i,k))      * invrs_hdtime ! Tendency of static energy
      end do
    end do

    invrs_macmic_num_steps = 1.0_kind_phys / REAL(cld_macmic_num_steps,kind_phys)

    do k = top_lev, pver
      do i = 1, ncol

        k_clubb = k + 1 - top_lev

        ! need to initialize macmic coupling to zero
        if ( macmic_it == 1 ) then
          ttend_clubb_mc_pbuf(i,k_clubb)     = 0._kind_phys
        end if

        !  Accumulate vars through macmic subcycle for Gravity Wave parameterization
        ttend_clubb_mc_pbuf(i,k_clubb) = ttend_clubb_mc_pbuf(i,k_clubb) + ptend_s(i,k) / cpair

        ! And average at last macmic step
        if (macmic_it == cld_macmic_num_steps) then
          ttend_clubb_pbuf(i,k)  = ttend_clubb_mc_pbuf(i,k_clubb) * invrs_macmic_num_steps
        end if

      end do
    end do

    do k = top_lev, pverp
      do i = 1, ncol

        k_clubb = k + 1 - top_lev

        ! need to initialize macmic coupling to zero
        if ( macmic_it == 1 ) then
          upwp_clubb_gw_mc_pbuf(i,k_clubb)   = 0._kind_phys
          vpwp_clubb_gw_mc_pbuf(i,k_clubb)   = 0._kind_phys
          thlp2_clubb_gw_mc_pbuf(i,k_clubb)  = 0._kind_phys
          wpthlp_clubb_gw_mc_pbuf(i,k_clubb) = 0._kind_phys
        end if

        !  Accumulate vars through macmic subcycle for Gravity Wave parameterization
        upwp_clubb_gw_mc_pbuf  (i,k_clubb) =   upwp_clubb_gw_mc_pbuf(i,k_clubb) + upwp_pbuf  (i,k_clubb)
        vpwp_clubb_gw_mc_pbuf  (i,k_clubb) =   vpwp_clubb_gw_mc_pbuf(i,k_clubb) + vpwp_pbuf  (i,k_clubb)
        thlp2_clubb_gw_mc_pbuf (i,k_clubb) =  thlp2_clubb_gw_mc_pbuf(i,k_clubb) + thlp2_pbuf (i,k_clubb)
        wpthlp_clubb_gw_mc_pbuf(i,k_clubb) = wpthlp_clubb_gw_mc_pbuf(i,k_clubb) + wpthlp_pbuf(i,k_clubb)

        ! And average at last macmic step
        if (macmic_it == cld_macmic_num_steps) then
          upwp_clubb_gw_pbuf  (i,k) =   upwp_clubb_gw_mc_pbuf(i,k_clubb) * invrs_macmic_num_steps
          vpwp_clubb_gw_pbuf  (i,k) =   vpwp_clubb_gw_mc_pbuf(i,k_clubb) * invrs_macmic_num_steps
          thlp2_clubb_gw_pbuf (i,k) =  thlp2_clubb_gw_mc_pbuf(i,k_clubb) * invrs_macmic_num_steps
          wpthlp_clubb_gw_pbuf(i,k) = wpthlp_clubb_gw_mc_pbuf(i,k_clubb) * invrs_macmic_num_steps
        end if

      end do
    end do

    if (clubb_do_adv) then
      if (macmic_it == cld_macmic_num_steps) then

        do k = top_lev, pver
          do i = 1, ncol

            k_clubb = k + 1 - top_lev

            thlp2_pbuf(i,k_clubb) = max( thl_tol**2, thlp2_pbuf(i,k_clubb) )
            rtp2_pbuf (i,k_clubb) = max(  rt_tol**2,  rtp2_pbuf(i,k_clubb) )
            wp2_pbuf  (i,k_clubb) = max(  w_tol_sqd,   wp2_pbuf(i,k_clubb) )
            up2_pbuf  (i,k_clubb) = max(  w_tol_sqd,   up2_pbuf(i,k_clubb) )
            vp2_pbuf  (i,k_clubb) = max(  w_tol_sqd,   vp2_pbuf(i,k_clubb) )

            ! Here add a constant to moments which can be either positive or
            !  negative.  This is to prevent clipping when dynamics tries to
            !  make all constituents positive
            wp3_pbuf    (i,k_clubb) =     wp3_pbuf(i,k_clubb) +     wp3_const
            rtpthlp_pbuf(i,k_clubb) = rtpthlp_pbuf(i,k_clubb) + rtpthlp_const
            wpthlp_pbuf (i,k_clubb) =  wpthlp_pbuf(i,k_clubb) +  wpthlp_const
            wprtp_pbuf  (i,k_clubb) =   wprtp_pbuf(i,k_clubb) +   wprtp_const

            ptend_q(i,k,ixrtpthlp) = (rtpthlp_pbuf(i,k_clubb) - state_q(i,k,ixrtpthlp) ) * invrs_hdtime ! RTP THLP covariance
            ptend_q(i,k,ixwpthlp)  = ( wpthlp_pbuf(i,k_clubb) - state_q(i,k,ixwpthlp)  ) * invrs_hdtime ! WPTHLP
            ptend_q(i,k,ixwprtp)   = (  wprtp_pbuf(i,k_clubb) - state_q(i,k,ixwprtp)   ) * invrs_hdtime ! WPRTP
            ptend_q(i,k,ixwp3)     = (    wp3_pbuf(i,k_clubb) - state_q(i,k,ixwp3)     ) * invrs_hdtime ! WP3
            ptend_q(i,k,ixwp2)     = (    wp2_pbuf(i,k_clubb) - state_q(i,k,ixwp2)     ) * invrs_hdtime ! WP2
            ptend_q(i,k,ixthlp2)   = (  thlp2_pbuf(i,k_clubb) - state_q(i,k,ixthlp2)   ) * invrs_hdtime ! THLP Variance
            ptend_q(i,k,ixrtp2)    = (   rtp2_pbuf(i,k_clubb) - state_q(i,k,ixrtp2)    ) * invrs_hdtime ! RTP Variance
            ptend_q(i,k,ixup2)     = (    up2_pbuf(i,k_clubb) - state_q(i,k,ixup2)     ) * invrs_hdtime ! UP2
            ptend_q(i,k,ixvp2)     = (    vp2_pbuf(i,k_clubb) - state_q(i,k,ixvp2)     ) * invrs_hdtime ! VP2

          end do
        end do

      end if
    end if

    !  Apply tendencies to ice mixing ratio, liquid and ice number, and aerosol constituents that aren't mixed by ndrop
    !  Loading up this array doesn't mean the tendencies are applied.
    ! edsclr is compressed with just the constituents being used, ptend and state are not compressed
    icnt=0
    do ixind = 1, pcnst
      if (lq(ixind)) then
        icnt=icnt+1
        if ((ixind /= ixq)       .and. (ixind /= ixcldliq) .and.&
            (ixind /= ixthlp2)   .and. (ixind /= ixrtp2)   .and.&
            (ixind /= ixrtpthlp) .and. (ixind /= ixwpthlp) .and.&
            (ixind /= ixwprtp)   .and. (ixind /= ixwp2)    .and.&
            (ixind /= ixwp3)     .and. (ixind /= ixup2)    .and. (ixind /= ixvp2) ) then

          ! Zero out levels above top_lev
          do k = 1, top_lev-1
            do i = 1, ncol
              ptend_q(i,k,ixind) = 0._kind_phys
            end do
          end do

          ! Copy CLUBB's edsclr values
          do k = top_lev, pver
            do i = 1, ncol
              k_clubb = k + 1 - top_lev
              ptend_q(i,k,ixind) = (edsclr(i,k_clubb,icnt)-state_q(i,k,ixind)) / hdtime ! transported constituents
            end do
          end do

        end if
      end if
    end do

    ! Cleanup err_info
    call cleanup_err_info_api(err_info)

  contains

    subroutine compute_mf_plumes()
      implicit none

      ! Local variables restricted to the MF diagnostic plumes block
      real(kind_phys), parameter             :: cflval = 1._kind_phys    ! Upper limit on Courant number for stability
      real(kind_phys)                        :: lambda
      character(len=256)                     :: lambda_errmsg
      real(kind_phys), dimension(ncol) :: cflfac, max_cfl, th_sfc
      logical                                :: cfllim
      integer                                :: i, k

      ! fields on momentum grid needed for mass flux calc.
      rtm_zm     = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,  rtm(:ncol,:) )
      thlm_zm    = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, thlm(:ncol,:) )
      th_zm      = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, th_zt(:ncol,:) )
      qv_zm      = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, qv_zt(:ncol,:) )
      qc_zm      = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, qc_zt(:ncol,:) )
      kappa_zm   = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, kappa_zt(:ncol,:) )

      do k = 1, nzm_clubb
         do i = 1, ncol
            invrs_exner_zm(i,k) = 1._kind_phys / ( (p_in_Pa_zm(i,k) * inv_p0_clubb)**kappa_zm(i,k) )
         end do
      end do

      ! theta at surface needed for plume initialization.
      ! gr%k_lb_zm is the surface (lower boundary) momentum index: 1 for an
      ! ascending grid, nzm_clubb for a descending grid.
      th_sfc(1:ncol) = ts(1:ncol)*invrs_exner_zm(1:ncol,gr%k_lb_zm)

      ! friction velocity needed for clubb-mf
      ustar2(1:ncol) = calc_friction_velocity(wsx(1:ncol), wsy(1:ncol), rrho(1:ncol))

      ! update thv if clubb is subcycled
      if (tt>1) then
         do i = 1, ncol
            thv_ds_zt(i,:) = thlm(i,:) &
                 + latvap*rcm(i,:)*invrs_exner_zt(i,:)/cpair
            thv_ds_zt(i,:) = thv_ds_zt(i,:) &
                 * (1._kind_phys+zvir*(rtm(i,:)-rcm(i,:)) - rcm(i,:))
         end do
         thv_ds_zm(1:ncol,:) = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,   thv_ds_zt(1:ncol,:) )
      end if

      ! set memory terms
      mf_ztopm1(1:ncol,:) = ztopma(1:ncol,:)
      mf_ddcp(1:ncol,:) = ddcp(1:ncol,:)
      mf_cbm1(1:ncol) = cbm1(1:ncol)

      !--------------------------------------- integrate_mf call ---------------------------------------
      ! integrate_mf expects arguments of individual columns.
      ! If the column loop gets pushed into it, we can also avoid the array slicing.

      do i = 1, ncol

         call integrate_mf( nzm_clubb,      nzt_clubb,       dtime,                                   & ! input
              cpair,          gravit,          latvap,          latice,          zvir,                & ! input
              h2otrip,                         shr_const_pi,                                          & ! input
              rho_zm(i,:),    zi_g(i,:),                        p_in_Pa_zm(i,:), invrs_exner_zm(i,:), & ! input
              rho_zt(i,:),    dz_g(i,:),       zt_g(i,:),       p_in_Pa(i,:),    invrs_exner_zt(i,:), & ! input
              um(i,:),        vm(i,:),         thlm(i,:),       rtm(i,:),        thv_ds_zt(i,:),      & ! input
              wm_zm(i,:),      th_zt(i,:),      qv_zt(i,:),      qc_zt(i,:),          & ! input
              thlm_zm(i,:),   rtm_zm(i,:),     thv_ds_zm(i,:),                                        & ! input
              th_zm(i,:),     qv_zm(i,:),      qc_zm(i,:),                                            & ! input
              ustar2(i),      th_sfc(i),       wpthlp_sfc(i),   wprtp_sfc(i),    pblh_pbuf(i),        & ! input
              tke_zm(i,:),    tpert(i),        mf_lts(i),                                             & ! input
              wpthlp_pbuf(i,:),                wpthvp_pbuf(i,:),                 wprtp_pbuf(i,:),     & ! input
              mf_ztopm1(i,:),                  mf_ddcp(i,:),                     mf_cbm1(i),          & ! in-out
              mf_cape(i,:),                                                                           & ! output
              mf_upa(i,:,:),    mf_dna(i,:,:),                                                        & ! output
              mf_upw(i,:,:),    mf_dnw(i,:,:),                                                        & ! output
              mf_upmf(i,:,:),                                                                         & ! output
              mf_upqt(i,:,:),   mf_dnqt(i,:,:),                                                       & ! output
              mf_upthl(i,:,:),  mf_dnthl(i,:,:),                                                      & ! output
              mf_upthv(i,:,:),  mf_dnthv(i,:,:),                                                      & ! output
              mf_upth(i,:,:),   mf_dnth(i,:,:),                                                       & ! output
              mf_upqc(i,:,:),   mf_dnqc(i,:,:),                                                       & ! output
              mf_upbuoy(i,:,:),                                                                       & ! output
              mf_upent(i,:,:),                                                                        & ! output
              mf_updet(i,:,:),                                                                        & ! output
              mf_dry_a(i,:),    mf_moist_a(i,:),                                                      & ! output
              mf_dry_w(i,:),    mf_moist_w(i,:),                                                      & ! output
              mf_dry_qt(i,:),   mf_moist_qt(i,:),                                                     & ! output
              mf_dry_thl(i,:),  mf_moist_thl(i,:),                                                    & ! output
              mf_dry_u(i,:),    mf_moist_u(i,:),                                                      & ! output
              mf_dry_v(i,:),    mf_moist_v(i,:),                                                      & ! output
              mf_moist_qc(i,:),                                                     & ! output
              s_ae(i,:),                                                                              & ! output
              s_ac(i,:),        s_aup(i,:),      s_adn(i,:),                                          & ! output
              s_aw(i,:),        s_awup(i,:),     s_awdn(i,:),                                         & ! output
              s_aww(i,:),       s_awwup(i,:),    s_awwdn(i,:),                                        & ! output
              s_awthlup(i,:),   s_awqtup(i,:),   s_awuup(i,:),   s_awvup(i,:),                        & ! output
              s_awthldn(i,:),   s_awqtdn(i,:),   s_awudn(i,:),   s_awvdn(i,:),                        & ! output
              s_awthl(i,:),     s_awqt(i,:),                                                          & ! output
              s_awu(i,:),       s_awv(i,:),                                                           & ! output
              mf_thlflxup(i,:), mf_qtflxup(i,:), mf_uflxup(i,:), mf_vflxup(i,:),                      & ! output
              mf_thlflxdn(i,:), mf_qtflxdn(i,:), mf_uflxdn(i,:), mf_vflxdn(i,:),                      & ! output
              mf_thlflx(i,:),   mf_qtflx(i,:),   mf_uflx(i,:),   mf_vflx(i,:),                        & ! output
              mf_thvflx(i,:),                                                                         & ! output
              mf_sqtup(i,:),    mf_sthlup(i,:),                                                       & ! output
              mf_sqtdn(i,:),    mf_sthldn(i,:),                                                       & ! output
              mf_sqt(i,:),      mf_sthl(i,:),                                                         & ! output
              mf_sqtac(i,:),    mf_sqtev(i,:),                                                        & ! output
              mf_precc(i,:),                                                                          & ! output
              mf_ztop(i,:),     mf_L0(i,:),                                                           & ! output
              s_mfup(i,:),      s_entup(i,:),    s_detup(i,:),                                        & ! output
              s_mfdn(i,:),      s_entdn(i,:),    s_detdn(i,:),                                        & ! output
              mf_kctop(i) )                                                                             ! output

      end do

      !--------------------------------------- END integrate_mf call ---------------------------------------

      ! CFL limiter
      cfllim = .true.
      cflfac(:ncol) = 1._kind_phys
      ! no mass flux through the ground; gr%k_lb_zm is the surface momentum
      ! index (1 ascending, nzm_clubb descending)
      s_aw(:ncol,gr%k_lb_zm) = 0._kind_phys
      max_cfl(:ncol)= 0._kind_phys
      do i=1,ncol
         do k=1,nzt_clubb
            max_cfl(i) = max(max_cfl(i),dtime*invrs_dz_g(i,k)*max(abs(s_aw(i,k)),abs(s_aw(i,k+1))))
         end do
         if (max_cfl(i).gt.cflval.and.cfllim) cflfac(i) = cflval/max_cfl(i)
      end do

      ! check for mass conservation
      do k=1,nzt_clubb
         do i=1,ncol
            if ((-1._kind_phys*mf_sqt(i,k)*dtime) > rtm(i,k)) then
               lambda = -1._kind_phys*rtm(i,k)/(mf_sqt(i,k)*dtime)
               write(lambda_errmsg,'(a,3i6,4es16.8)') &
                    'MF drying exceeds rtm. i,k,lambda,mf_sqt,rtm,dtime = ', &
                     i, k, lambda, mf_sqt(i,k), rtm(i,k), dtime
               errmsg = 'clubb_tend_cam: '//trim(lambda_errmsg)
               errflg = 1
               return
            end if
         end do
      end do

      ! pass MF turbulent advection term as CLUBB explicit forcing term
      rtm_forcing(:ncol,:)  = 0._kind_phys
      thlm_forcing(:ncol,:) = 0._kind_phys
      mf_qtforcup(:ncol,:)  = 0._kind_phys
      mf_thlforcup(:ncol,:) = 0._kind_phys
      mf_qtforcdn(:ncol,:)  = 0._kind_phys
      mf_thlforcdn(:ncol,:) = 0._kind_phys

      ! pass MF turbulent advection term as CLUBB explicit forcing term
      do k = 1, nzt_clubb
         do i = 1, ncol
            rtm_forcing(i,k)  = rtm_forcing(i,k) - invrs_rho_ds_zt(i,k) * gr%grid_dir * invrs_dz_g(i,k) * cflfac(i) * &
                 ((rho_ds_zm(i,k+1) * mf_qtflx(i,k+1)) - (rho_ds_zm(i,k) * mf_qtflx(i,k))) &
                 + mf_sqt(i,k)

            thlm_forcing(i,k) = thlm_forcing(i,k) - invrs_rho_ds_zt(i,k) * gr%grid_dir * invrs_dz_g(i,k) * cflfac(i) * &
                 ((rho_ds_zm(i,k+1) * mf_thlflx(i,k+1)) - (rho_ds_zm(i,k) * mf_thlflx(i,k))) &
                 + mf_sthl(i,k)

            mf_qtforcup(i,k)  = mf_qtforcup(i,k) - invrs_rho_ds_zt(i,k) * gr%grid_dir * invrs_dz_g(i,k) * cflfac(i) * &
                 ((rho_ds_zm(i,k+1) * mf_qtflxup(i,k+1)) - (rho_ds_zm(i,k) * mf_qtflxup(i,k))) &
                 + mf_sqtup(i,k)

            mf_thlforcup(i,k) = mf_thlforcup(i,k) - invrs_rho_ds_zt(i,k) * gr%grid_dir * invrs_dz_g(i,k) * cflfac(i) * &
                 ((rho_ds_zm(i,k+1) * mf_thlflxup(i,k+1)) - (rho_ds_zm(i,k) * mf_thlflxup(i,k))) &
                 + mf_sthlup(i,k)

            mf_qtforcdn(i,k)  = mf_qtforcdn(i,k) - invrs_rho_ds_zt(i,k) * gr%grid_dir * invrs_dz_g(i,k) * cflfac(i) * &
                 ((rho_ds_zm(i,k+1) * mf_qtflxdn(i,k+1)) - (rho_ds_zm(i,k) * mf_qtflxdn(i,k))) &
                 + mf_sqtdn(i,k)

            mf_thlforcdn(i,k) = mf_thlforcdn(i,k) - invrs_rho_ds_zt(i,k) * gr%grid_dir * invrs_dz_g(i,k) * cflfac(i) * &
                 ((rho_ds_zm(i,k+1) * mf_thlflxdn(i,k+1)) - (rho_ds_zm(i,k) * mf_thlflxdn(i,k))) &
                 + mf_sthldn(i,k)
         end do
      end do

      ! convective momentum transport
      if (do_clubb_mf_cmt) then
         um_forcing(:ncol,:) = 0._kind_phys
         vm_forcing(:ncol,:) = 0._kind_phys
         do k = 1, nzt_clubb
            do i = 1, ncol
               um_forcing(i,k)   = um_forcing(i,k) - invrs_rho_ds_zt(i,k) * gr%grid_dir * invrs_dz_g(i,k) * cflfac(i) * &
                    ((rho_ds_zm(i,k+1) * mf_uflx(i,k+1)) - (rho_ds_zm(i,k) * mf_uflx(i,k)))

               vm_forcing(i,k)   = vm_forcing(i,k) - invrs_rho_ds_zt(i,k) * gr%grid_dir * invrs_dz_g(i,k) * cflfac(i) * &
                    ((rho_ds_zm(i,k+1) * mf_vflx(i,k+1)) - (rho_ds_zm(i,k) * mf_vflx(i,k)))
            end do
         end do
      end if

      do i=1,ncol
         ! compute ensemble cloud properties
         mf_qc_nadv(i,:nzm_clubb)        = mf_qc_nadv(i,:nzm_clubb) + mf_moist_qc(i,:nzm_clubb)
         mf_rcm_nadv(i,:nzm_clubb)       = mf_rcm_nadv(i,:nzm_clubb) + mf_moist_a(i,:nzm_clubb)*mf_moist_qc(i,:nzm_clubb)
         mf_cloudfrac_nadv(i,:nzm_clubb) = mf_cloudfrac_nadv(i,:nzm_clubb) + mf_moist_a(i,:nzm_clubb)

         ! surface precip, [kg/m2/s]->[m/s]
         mf_precc_nadv(i) = mf_precc_nadv(i) + mf_precc(i,gr%k_lb_zm)/1000._kind_phys
         mf_snow_nadv(i)  = 0._kind_phys

         ! accumulate the deep-hookup variables over the nadv subcycles.
         ! the fractional entrainment/detrainment are RATIO diagnostics
         ! (mass-flux-weighted per subcycle, denominator = mfup): near
         ! plume tops the per-subcycle ratio can be immense while the
         ! same subcycle's mass flux vanishes -- averaging the ratio
         ! separately from the mass flux breaks that cancellation and
         ! produced enormous ZM_EU values (aero_convproc courant
         ! substep count exploded -> nstep-0 hang).  accumulate the
         ! entrainment/detrainment MASS FLUX (rate x mf, bounded)
         ! instead; the population block consumes it directly.
         s_mfup_nadv(i,:nzm_clubb)   = s_mfup_nadv(i,:nzm_clubb)   + s_mfup(i,:nzm_clubb)
         s_entup_nadv(i,:nzm_clubb)  = s_entup_nadv(i,:nzm_clubb)  + s_entup(i,:nzm_clubb)*s_mfup(i,:nzm_clubb)
         s_detup_nadv(i,:nzm_clubb)  = s_detup_nadv(i,:nzm_clubb)  + s_detup(i,:nzm_clubb)*s_mfup(i,:nzm_clubb)
         s_mfdn_nadv(i,:nzm_clubb)   = s_mfdn_nadv(i,:nzm_clubb)   + s_mfdn(i,:nzm_clubb)
         s_entdn_nadv(i,:nzm_clubb)  = s_entdn_nadv(i,:nzm_clubb)  + s_entdn(i,:nzm_clubb)*abs(s_mfdn(i,:nzm_clubb))
         s_detdn_nadv(i,:nzm_clubb)  = s_detdn_nadv(i,:nzm_clubb)  + s_detdn(i,:nzm_clubb)*abs(s_mfdn(i,:nzm_clubb))
         ! accumulate sum(a*w^2) and sum(a*w) of the updraft
         ! ensemble so the ratio gives the mass-flux-weighted plume w
         s_awwup_nadv(i,:nzm_clubb)  = s_awwup_nadv(i,:nzm_clubb)  + s_awwup(i,:nzm_clubb)
         s_awup_nadv(i,:nzm_clubb)   = s_awup_nadv(i,:nzm_clubb)   + s_awup(i,:nzm_clubb)
         mf_sqtac_nadv(i,:nzt_clubb) = mf_sqtac_nadv(i,:nzt_clubb) + mf_sqtac(i,:nzt_clubb)
         mf_sqtev_nadv(i,:nzt_clubb) = mf_sqtev_nadv(i,:nzt_clubb) + mf_sqtev(i,:nzt_clubb)
         ! kctop is not averaged: keep the HIGHEST plume top over the
         ! subcycles (it doubles as the gathering/activity criterion:
         ! int(kctop) > 1 means the column was active at some sub-cycle)
         mf_kctop_nadv(i)            = max(mf_kctop_nadv(i), mf_kctop(i))

         mf_L0_nadv(i,:)     = mf_L0_nadv(i,:) + mf_L0(i,:)
         mf_ztop_nadv(i,:)   = mf_ztop_nadv(i,:) + mf_ztop(i,:)
         mf_ztopm1_nadv(i,:) = mf_ztopm1_nadv(i,:) + mf_ztopm1(i,:)
         mf_cape_nadv(i,:)   = mf_cape_nadv(i,:) + mf_cape(i,:)
         mf_ddcp_nadv(i,:)   = mf_ddcp_nadv(i,:) + mf_ddcp(i,:)
         mf_cbm1_nadv(i)     = mf_cbm1_nadv(i) + mf_cbm1(i)

         if (ANY(mf_ztop(i,:) > 0._kind_phys)) mf_freq_nadv(i) = mf_freq_nadv(i) + 1._kind_phys

         mf_thlforcup_nadv(i,:nzt_clubb) = mf_thlforcup_nadv(i,:nzt_clubb) + mf_thlforcup(i,:nzt_clubb)
         mf_qtforcup_nadv(i,:nzt_clubb)  = mf_qtforcup_nadv(i,:nzt_clubb) + mf_qtforcup(i,:nzt_clubb)
         mf_thlforcdn_nadv(i,:nzt_clubb) = mf_thlforcdn_nadv(i,:nzt_clubb) + mf_thlforcdn(i,:nzt_clubb)
         mf_qtforcdn_nadv(i,:nzt_clubb)  = mf_qtforcdn_nadv(i,:nzt_clubb) + mf_qtforcdn(i,:nzt_clubb)
         mf_thlforc_nadv(i,:nzt_clubb) = mf_thlforc_nadv(i,:nzt_clubb) + thlm_forcing(i,:nzt_clubb)
         mf_qtforc_nadv(i,:nzt_clubb)  = mf_qtforc_nadv(i,:nzt_clubb) + rtm_forcing(i,:nzt_clubb)
         mf_ent_nadv(i,:nzm_clubb)     = mf_ent_nadv(i,:nzm_clubb) + s_awu(i,:nzm_clubb)
         max_cfl_nadv(i)               = MAX(max_cfl(i),max_cfl_nadv(i))
      end do

    end subroutine compute_mf_plumes

  end subroutine clubb1_run

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

  subroutine clubb2_run(ncol, pver, ixcldliq, ixcldice, ixnumliq, ixnumice, & ! in
                        clubb_detliq_rad, clubb_detice_rad, clubb_detphase_lowtemp, &! in
                        meltpt_temp, latice, rga, & ! in
                        dlf, state_t, state_pdel, state_pdeldry, & ! in
                        ptend_q, s, det_s, det_ice, & ! inout
                        dlf_liq_out, dlf_ice_out, & ! out
                        errmsg, errflg ) ! out

    use ccpp_kinds, only: kind_phys

    ! Input variables, intent(in)
    integer, intent(in) :: ncol, pver, ixcldliq, ixcldice, ixnumliq, ixnumice
    real(kind_phys), intent(in) :: clubb_detliq_rad, clubb_detice_rad, clubb_detphase_lowtemp
    real(kind_phys), intent(in) :: meltpt_temp, latice, rga
    real(kind_phys), intent(in) :: dlf(:,:), state_t(:,:), state_pdel(:,:), state_pdeldry(:,:)

    ! Input variables, intent(inout)
    real(kind_phys), intent(inout) :: ptend_q(:,:,:)
    real(kind_phys), intent(inout) :: s(:,:)
    real(kind_phys), intent(inout) :: det_s(:), det_ice(:)

    ! Input variables, intent(out)
    real(kind_phys), intent(out) :: dlf_liq_out(:,:), dlf_ice_out(:,:)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Local variables
    integer :: i, k
    real(kind_phys) :: dlf2, dum1, dl_rad, di_rad, dt_low

    errmsg = ''
    errflg = 0

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

        if( state_t(i,k) > meltpt_temp ) then
          dum1 = 0.0_kind_phys
        elseif ( state_t(i,k) < dt_low ) then
          dum1 = 1.0_kind_phys
        else
          dum1 = ( meltpt_temp - state_t(i,k) ) / ( meltpt_temp - dt_low )
        endif

        ptend_q(i,k,ixcldliq) = dlf(i,k) * ( 1._kind_phys - dum1 )
        ptend_q(i,k,ixcldice) = dlf(i,k) * dum1
        ptend_q(i,k,ixnumliq) = 3._kind_phys * ( max(0._kind_phys, ( dlf(i,k) - dlf2 )) * ( 1._kind_phys - dum1 ) ) &
                                   / (4._kind_phys*3.14_kind_phys*dl_rad**3*997._kind_phys) + & ! Deep    Convection
                                   3._kind_phys * (                         dlf2    * ( 1._kind_phys - dum1 ) ) &
                                   / (4._kind_phys*3.14_kind_phys*10.e-6_kind_phys**3*997._kind_phys)     ! Shallow Convection
        ptend_q(i,k,ixnumice) = 3._kind_phys * ( max(0._kind_phys, ( dlf(i,k) - dlf2 )) *  dum1 ) &
                                   / (4._kind_phys*3.14_kind_phys*di_rad**3*500._kind_phys) + & ! Deep    Convection
                                   3._kind_phys * (                         dlf2    *  dum1 ) &
                                   / (4._kind_phys*3.14_kind_phys*50.e-6_kind_phys**3*500._kind_phys)     ! Shallow Convection
        s(i,k)          = dlf(i,k) * dum1 * latice

        dlf_liq_out(i,k) = dlf(i,k) * ( 1._kind_phys - dum1 )
        dlf_ice_out(i,k) = dlf(i,k) * dum1

        ! convert moist dlf tendencies to dry
        ptend_q(i,k,ixcldliq) = ptend_q(i,k,ixcldliq)*state_pdel(i,k)/state_pdeldry(i,k)
        ptend_q(i,k,ixcldice) = ptend_q(i,k,ixcldice)*state_pdel(i,k)/state_pdeldry(i,k)

        ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
        !   track of the integrals of ice and static energy that is effected from conversion to ice
        !   so that the energy checker doesn't complain.
        det_s(i)                  = det_s(i)   + s(i,k)          * state_pdel(i,k)    * rga
        det_ice(i)                = det_ice(i) - ptend_q(i,k,ixcldice) * state_pdeldry(i,k) * rga
      enddo
    enddo

    det_ice(:ncol) = det_ice(:ncol) / 1000._kind_phys  ! divide by density of water


  end subroutine clubb2_run

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

  subroutine clubb3_run(ncol, pver, pverp, pcnst, top_lev, & ! in
                        ixq, ixcldice, ixcldliq, ixnumice, & ! in
                        rhminis_const, rhmaxis_const, rhmini_const, rhmaxi_const, & ! in
                        rhminl_const, rhminl_adj_land_const, rhminh_const, & ! in
                        dp1, dp2, zvir, rair, cpair, gravit, karman, & ! in
                        state_t, state_q, & ! in
                        state_pmid, state_pdel, state_pdeldry, & ! in
                        state_zm, state_zi, state_u, state_v, state_exner, & ! in
                        landfrac, snowhland, wsx, & ! in
                        wsy, shf, cflx, ptend_q, & ! in
                        rcm, cloud_frac, cmfmc, & ! in
                        single_column, scm_cambfb_mode, lq, & ! in
                        cnst_type, scm_clubb_iop_name, subcol_scheme, & ! in
                        mf_cloudfrac_output, mf_qc_output, & ! in
                        pblh_pbuf, alst_pbuf, qlst_pbuf, deepcu_pbuf, shalcu_pbuf, & ! inout
                        cmfmc_sh_pbuf, dp_icwmr_pbuf, concld_pbuf, aist_pbuf, & ! inout
                        qsatfac_pbuf, ast_pbuf, qist_pbuf, cld_pbuf, troplev, & ! inout
                        errmsg, errflg ) ! out

    use ccpp_kinds, only: kind_phys
    use clubb_mf,              only: do_clubb_mf, do_clubb_mf_rad, clubb_mf_cldfrac_fac
    use holtslag_boville_diff, only: hb_pbl_dependent_coefficients_run
    use compute_cloud_fraction_two_moment, only: aist_vector
    use atmos_phys_pbl_utils,  only: calc_friction_velocity, calc_obukhov_length, calc_ideal_gas_rrho, &
                                     calc_kinematic_heat_flux, calc_kinematic_water_vapor_flux, &
                                     calc_kinematic_buoyancy_flux

    ! Input variables, intent(in)
    integer, intent(in) :: ncol, pver, pverp, pcnst, top_lev
    integer, intent(in) :: ixq, ixcldice, ixcldliq, ixnumice
    real(kind_phys), intent(in) :: rhminis_const, rhmaxis_const, rhmini_const, rhmaxi_const
    real(kind_phys), intent(in) :: rhminl_const, rhminl_adj_land_const, rhminh_const
    real(kind_phys), intent(in) :: dp1, dp2, zvir, rair, cpair, gravit, karman
    real(kind_phys), intent(in) :: landfrac(:), snowhland(:), &
                                   wsx(:), wsy(:), shf(:)
    real(kind_phys), intent(in) :: state_pmid(:,:), state_pdel(:,:), state_pdeldry(:,:), rcm(:,:), &
                                   cloud_frac(:,:), state_t(:,:), state_exner(:,:), &
                                   state_zm(:,:), state_zi(:,:), state_u(:,:), state_v(:,:), cmfmc(:,:), cflx(:,:), &
                                   mf_cloudfrac_output(:,:), mf_qc_output(:,:)
    real(kind_phys), intent(in) :: state_q(:,:,:)
   ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    !real(kind_phys), intent(in)         :: tropp_p_loc(:,:)
    logical, intent(in) :: single_column, scm_cambfb_mode
    logical, intent(in) :: lq(:)
    character(len=3), intent(in) :: cnst_type(:)
    character(len=20), intent(in) :: scm_clubb_iop_name
    character(len=16), intent(in) :: subcol_scheme

    ! Input variables, intent(inout)
    real(kind_phys), intent(inout) :: pblh_pbuf(:)
    real(kind_phys), intent(inout) :: alst_pbuf(:,:), qlst_pbuf(:,:), deepcu_pbuf(:,:), shalcu_pbuf(:,:), &
                                      cmfmc_sh_pbuf(:,:), dp_icwmr_pbuf(:,:), concld_pbuf(:,:), &
                                      aist_pbuf(:,:), qsatfac_pbuf(:,:), ast_pbuf(:,:), qist_pbuf(:,:), &
                                      cld_pbuf(:,:)
    real(kind_phys), intent(inout) :: ptend_q(:,:,:)

    ! Input variables, intent(out)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Local variables
    integer :: i, k, ixind, k_clubb
    integer, intent(in) :: troplev(:)
    real(kind_phys) :: rhmini(ncol), rhmaxi(ncol)
    real(kind_phys) :: frac_limit, ic_limit
    real(kind_phys) :: rrho(ncol), ustar2(ncol), kinheat(ncol), kinwat(ncol), kbfs(ncol), obklen(ncol), &
                       dummy2(ncol), dummy3(ncol)
    real(kind_phys) :: th(ncol,pver), thv(ncol,pver)
    real(kind_phys) :: rhmini_default(ncol)
    real(kind_phys) :: rhmaxi_default(ncol)
    real(kind_phys) :: rhminl_arr(ncol)
    real(kind_phys) :: rhminl_adj_land_arr(ncol)
    real(kind_phys) :: rhminh_arr(ncol)

    ! ---------------------------------------------------------

    errmsg = ''
    errflg = 0

    ! ptend_all now has all accumulated tendencies.  Convert the tendencies for the
    ! wet constituents to wet air basis.
    do ixind = 1, pcnst
      if (lq(ixind) .and. cnst_type(ixind) == 'wet') then
        do k = 1, pver
          do i = 1, ncol
            ptend_q(i,k,ixind) = ptend_q(i,k,ixind)*state_pdeldry(i,k)/state_pdel(i,k)
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
    !  SH_FRAC has no other setter with shallow_scheme='CLUBB_SGS' (cam7
    !  convect_diagnostics zeroes ICWMRSH etc., but not SH_FRAC), and it is
    !  read downstream (e.g. conv_water), so it must be zeroed here.
    shalcu_pbuf(:,:) = 0.0_kind_phys

    do k = 1, pver-1
      do i = 1, ncol
        if (do_clubb_mf .and. do_clubb_mf_rad) then
          ! for the deep hookup the MF plume cloud is attributed to the DEEP fields
          ! (DP_FRAC/ICWMRDP feed aero_convproc and wetdep's conicw/cldcu the same
          ! way ZM's do); shalcu/sh_icwmr stay zero so nothing is double counted.
          deepcu_pbuf(i,k)   = clubb_mf_cldfrac_fac*mf_cloudfrac_output(i,k)
          dp_icwmr_pbuf(i,k) = mf_qc_output(i,k)
        else
          ! diagnose the deep convective cloud fraction, as done in macrophysics based on the
          ! deep convective mass flux, read in from pbuf.  Since shallow convection is never
          ! called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
          ! fraction is purely from deep convection scheme.
          deepcu_pbuf(i,k) = max(0.0_kind_phys,min(dp1*log(1.0_kind_phys+dp2*(cmfmc(i,k+1)-cmfmc_sh_pbuf(i,k+1))),0.6_kind_phys))
        end if

        if (deepcu_pbuf(i,k) <= frac_limit .or. dp_icwmr_pbuf(i,k) < ic_limit) then
          deepcu_pbuf(i,k) = 0._kind_phys
        endif

        ! NOTE: concld_pbuf used to be calculated in the commented-out version below, but since we
        ! set alst_pbuf=cloud_frac_pbuf, this simplifies to only using deepcu_pbuf.
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

    !       NOTE on macmic sub-cycling: with CLUBB_MF the fields DP_FRAC,
    !       ICWMRDP, CONCLD (and the deep part of CLD) computed above evolve
    !       every macmic iteration so in-loop consumers (PUMAS reads
    !       CLD/CONCLD) stay tightly coupled to the current sub-step.  Their
    !       post-loop consumers (aero_convproc/wetdep, conv_water, radiation)
    !       must instead see the sub-cycle MEAN; that averaging is done in the
    !       cam7 physpkg driver after the macmic loop, alongside
    !       prec_dp/snow_dp.

    ! --------------------------------------------------------------------------------- !
    !  COMPUTE THE ICE CLOUD FRACTION PORTION                                           !
    !  use the aist_vector function to compute the ice cloud fraction                   !
    ! --------------------------------------------------------------------------------- !

    aist_pbuf(:,:top_lev-1) = 0._kind_phys
    qsatfac_pbuf(:, :) = 0._kind_phys ! Zero out entire profile in case qsatfac is left undefined in aist_vector below

    rhmini_default(:) = rhmini_const
    rhmaxi_default(:) = rhmaxi_const
    rhminl_arr(:) = rhminl_const
    rhminl_adj_land_arr(:) = rhminl_adj_land_const
    rhminh_arr(:) = rhminh_const

    do k = top_lev, pver

      ! For Type II PSC and for thin cirrus, the clouds can be thin, but
      ! extensive and they should start forming when the gridbox mean saturation
      ! reaches 1.0.
      !
      ! For now, use the tropopause diagnostic to determine where the Type II
      ! PSC should be, but in the future wold like a better metric that can also
      ! identify the level for thin cirrus. Include the tropopause level so that
      ! the cold point tropopause will use the stratospheric values.
      where (k <= troplev)
        rhmini = rhminis_const
        rhmaxi = rhmaxis_const
      elsewhere
        rhmini = rhmini_const
        rhmaxi = rhmaxi_const
      end where

      !REMOVECAM: this is no longer needed when CAM is retired and pcols no longer exists
      aist_pbuf(:,k) = 0._kind_phys
      !REMOVECAM_END
      if ( trim(subcol_scheme) == 'SILHS' ) then
        call aist_vector(state_q(:ncol,k,ixq), state_t(:ncol,k), &
             state_pmid(:ncol,k), state_q(:ncol,k,ixcldice), &
             state_q(:ncol,k,ixnumice), landfrac(:ncol), &
             snowhland(:ncol), aist_pbuf(:ncol,k), ncol, &
             rhmaxi_in=rhmaxi_default(:ncol), &
             rhmini_in=rhmini_default(:ncol), &
             rhminl_in=rhminl_arr(:ncol), &
             rhminl_adj_land_in=rhminl_adj_land_arr(:ncol), &
             rhminh_in=rhminh_arr(:ncol))
      else
        call aist_vector(state_q(:ncol,k,ixq), state_t(:ncol,k), &
             state_pmid(:ncol,k), state_q(:ncol,k,ixcldice), &
             state_q(:ncol,k,ixnumice), landfrac(:ncol), &
             snowhland(:ncol), aist_pbuf(:ncol,k), ncol, &
             rhmaxi_in=rhmaxi(:ncol), &
             rhmini_in=rhmini(:ncol), &
             rhminl_in=rhminl_arr(:ncol), &
             rhminl_adj_land_in=rhminl_adj_land_arr(:ncol), &
             rhminh_in=rhminh_arr(:ncol), &
             qsatfac_out=qsatfac_pbuf(:ncol,k))
      endif
    enddo

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
         th(i,k) = state_t(i,k)*state_exner(i,k)
         !thv should have condensate loading to be consistent with earlier def's in this module
         thv(i,k) = th(i,k)*(1.0_kind_phys+zvir*state_q(i,k,ixq) - state_q(i,k,ixcldliq))
      enddo
    enddo

    ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)
    rrho   (1:ncol) = calc_ideal_gas_rrho(rair, state_t(1:ncol,pver), state_pmid(1:ncol,pver))
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
      z         = state_zm(:ncol,:pver),                    &
      zi        = state_zi(:ncol,:pverp),                   &
      u         = state_u(:ncol,:pver),                     &
      v         = state_v(:ncol,:pver),                     &
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

  end subroutine clubb3_run

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

  subroutine stats_init_clubb( pcols, l_stats_in, stats_tsamp_in, stats_tout_in, &
                               nnzp, nnrad_zt,nnrad_zm, delt, &
                               edsclr_dim, sclr_dim, hydromet_dim, &
                               hm_metadata, stats_metadata, &
                               clubb_vars_zt, clubb_vars_zm, clubb_vars_sfc, &
                               clubb_vars_rad_zt, clubb_vars_rad_zm, &
                               stats_zt, stats_zm, stats_sfc, &
                               stats_rad_zt, stats_rad_zm, &
                               errmsg, errflg )
    !
    ! Description: Initializes the statistics saving functionality of
    !   the CLUBB model.  This is for purpose of CAM-CLUBB interface.  Here
    !   the traditional stats_init of CLUBB is not called, as it is not compatible
    !   with CAM output.

    !-----------------------------------------------------------------------

    use ccpp_kinds, only: kind_phys
    use clubb_api_module, only: time_precision, stats, stats_metadata_type, hm_metadata_type, &   !
                                nvarmax_zm, stats_init_zm_api, & !
                                nvarmax_zt, stats_init_zt_api, & !
                                nvarmax_rad_zt, stats_init_rad_zt_api, & !
                                nvarmax_rad_zm, stats_init_rad_zm_api, & !
                                nvarmax_sfc, stats_init_sfc_api, & !
                                fstderr, var_length !

    implicit none

    !----------------------- Input Variables -----------------------

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    integer, intent(in) :: pcols

    integer, intent(in) :: sclr_dim, edsclr_dim, hydromet_dim

    integer, intent(inout) :: errflg
    character(len=*), intent(inout) :: errmsg

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

    integer :: i, ntot, j
    integer :: ierr
    character(len=256) :: alloc_errmsg

    !----------------------- Begin Code -----------------------

    !  Initialize
    errmsg = ''
    errflg = 0
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
         errmsg = "stats_init_clubb:  CLUBB stats_tsamp must be an even multiple of the timestep"
         errflg = 1
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
         errmsg = "stats_init_clubb:  number of zt statistical variables exceeds limit"
         errflg = 1
         return
      endif

      stats_zt(j)%num_output_fields = ntot
      stats_zt(j)%kk = nnzp - 1

      allocate( stats_zt(j)%z( stats_zt(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%z: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if

      allocate( stats_zt(j)%accum_field_values( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%accum_field_values: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_zt(j)%accum_num_samples( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%accum_num_samples: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_zt(j)%l_in_update( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%l_in_update: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      call stats_zero( stats_zt(j)%kk, stats_zt(j)%num_output_fields, stats_zt(j)%accum_field_values, &
                       stats_zt(j)%accum_num_samples, stats_zt(j)%l_in_update )

      allocate( stats_zt(j)%file%grid_avg_var( stats_zt(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%file%grid_avg_var: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_zt(j)%file%z( stats_zt(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%file%z: "//trim(alloc_errmsg)
        errflg = 1
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
         errmsg = "stats_init_clubb:  number of zm statistical variables exceeds limit"
         errflg = 1
         return
      endif

      stats_zm(j)%num_output_fields = ntot
      stats_zm(j)%kk = nnzp

      allocate( stats_zm(j)%z( stats_zm(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zm%z: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if

      allocate( stats_zm(j)%accum_field_values( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zm%accum_field_values: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_zm(j)%accum_num_samples( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zm%accum_num_samples: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_zm(j)%l_in_update( 1, 1, stats_zm(j)%kk, stats_zm(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zm%l_in_update: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      call stats_zero( stats_zm(j)%kk, stats_zm(j)%num_output_fields, stats_zm(j)%accum_field_values, &
                       stats_zm(j)%accum_num_samples, stats_zm(j)%l_in_update )

      allocate( stats_zm(j)%file%grid_avg_var( stats_zm(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zm%file%grid_avg_var: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_zm(j)%file%z( stats_zm(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zm%file%z: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if

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
            errmsg = "stats_init_clubb:  number of rad_zt statistical variables exceeds limit"
            errflg = 1
            return
         endif

        stats_rad_zt(j)%num_output_fields = ntot
        stats_rad_zt(j)%kk = nnrad_zt

        allocate( stats_rad_zt(j)%z( stats_rad_zt(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
        if( ierr /= 0 ) then
          errmsg = "stats_init_clubb: Failed to allocate stats_rad_zt%z: "//trim(alloc_errmsg)
          errflg = 1
          return
        end if

        allocate( stats_rad_zt(j)%accum_field_values( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ), &
                  stat=ierr, errmsg=alloc_errmsg )
        if( ierr /= 0 ) then
          errmsg = "stats_init_clubb: Failed to allocate stats_rad_zt%accum_field_values: "//trim(alloc_errmsg)
          errflg = 1
          return
        end if
        allocate( stats_rad_zt(j)%accum_num_samples( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ), &
                  stat=ierr, errmsg=alloc_errmsg )
        if( ierr /= 0 ) then
          errmsg = "stats_init_clubb: Failed to allocate stats_rad_zt%accum_num_samples: "//trim(alloc_errmsg)
          errflg = 1
          return
        end if
        allocate( stats_rad_zt(j)%l_in_update( 1, 1, stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields ), &
                  stat=ierr, errmsg=alloc_errmsg )
        if( ierr /= 0 ) then
          errmsg = "stats_init_clubb: Failed to allocate stats_rad_zt%l_in_update: "//trim(alloc_errmsg)
          errflg = 1
          return
        end if

        call stats_zero( stats_rad_zt(j)%kk, stats_rad_zt(j)%num_output_fields, stats_rad_zt(j)%accum_field_values, &
                       stats_rad_zt(j)%accum_num_samples, stats_rad_zt(j)%l_in_update )

        allocate( stats_rad_zt(j)%file%grid_avg_var( stats_rad_zt(j)%num_output_fields ), &
                  stat=ierr, errmsg=alloc_errmsg )
        if( ierr /= 0 ) then
          errmsg = "stats_init_clubb: Failed to allocate stats_rad_zt%file%grid_avg_var: "//trim(alloc_errmsg)
          errflg = 1
          return
        end if
        allocate( stats_rad_zt(j)%file%z( stats_rad_zt(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
        if( ierr /= 0 ) then
          errmsg = "stats_init_clubb: Failed to allocate stats_rad_zt%file%z: "//trim(alloc_errmsg)
          errflg = 1
          return
        end if

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
            errmsg = "stats_init_clubb:  number of rad_zm statistical variables exceeds limit"
            errflg = 1
            return
         endif

         stats_rad_zm(j)%num_output_fields = ntot
         stats_rad_zm(j)%kk = nnrad_zm

         allocate( stats_rad_zm(j)%z( stats_rad_zm(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
         if( ierr /= 0 ) then
           errmsg = "stats_init_clubb: Failed to allocate stats_rad_zm%z: "//trim(alloc_errmsg)
           errflg = 1
           return
         end if

         allocate( stats_rad_zm(j)%accum_field_values( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ), &
                   stat=ierr, errmsg=alloc_errmsg )
         if( ierr /= 0 ) then
           errmsg = "stats_init_clubb: Failed to allocate stats_rad_zm%accum_field_values: "//trim(alloc_errmsg)
           errflg = 1
           return
         end if
         allocate( stats_rad_zm(j)%accum_num_samples( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ), &
                   stat=ierr, errmsg=alloc_errmsg )
         if( ierr /= 0 ) then
           errmsg = "stats_init_clubb: Failed to allocate stats_rad_zm%accum_num_samples: "//trim(alloc_errmsg)
           errflg = 1
           return
         end if
         allocate( stats_rad_zm(j)%l_in_update( 1, 1, stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields ), &
                   stat=ierr, errmsg=alloc_errmsg )
         if( ierr /= 0 ) then
           errmsg = "stats_init_clubb: Failed to allocate stats_rad_zm%l_in_update: "//trim(alloc_errmsg)
           errflg = 1
           return
         end if

         call stats_zero( stats_rad_zm(j)%kk, stats_rad_zm(j)%num_output_fields, stats_rad_zm(j)%accum_field_values, &
                       stats_rad_zm(j)%accum_num_samples, stats_rad_zm(j)%l_in_update )

         allocate( stats_rad_zm(j)%file%grid_avg_var( stats_rad_zm(j)%num_output_fields ), &
                   stat=ierr, errmsg=alloc_errmsg )
         if( ierr /= 0 ) then
           errmsg = "stats_init_clubb: Failed to allocate stats_rad_zm%file%grid_avg_var: "//trim(alloc_errmsg)
           errflg = 1
           return
         end if
         allocate( stats_rad_zm(j)%file%z( stats_rad_zm(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
         if( ierr /= 0 ) then
           errmsg = "stats_init_clubb: Failed to allocate stats_rad_zm%file%z: "//trim(alloc_errmsg)
           errflg = 1
           return
         end if

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
         errmsg = "stats_init_clubb:  number of sfc statistical variables exceeds limit"
         errflg = 1
         return
      endif

      stats_sfc(j)%num_output_fields = ntot
      stats_sfc(j)%kk = 1

      allocate( stats_sfc(j)%z( stats_sfc(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_sfc%z: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if

      allocate( stats_sfc(j)%accum_field_values( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_sfc%accum_field_values: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_sfc(j)%accum_num_samples( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_sfc%accum_num_samples: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_sfc(j)%l_in_update( 1, 1, stats_sfc(j)%kk, stats_sfc(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_sfc%l_in_update: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if

      call stats_zero( stats_sfc(j)%kk, stats_sfc(j)%num_output_fields, stats_sfc(j)%accum_field_values, &
                       stats_sfc(j)%accum_num_samples, stats_sfc(j)%l_in_update )

      allocate( stats_sfc(j)%file%grid_avg_var( stats_sfc(j)%num_output_fields ), &
                stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_sfc%file%grid_avg_var: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if
      allocate( stats_sfc(j)%file%z( stats_sfc(j)%kk ), stat=ierr, errmsg=alloc_errmsg )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_sfc%file%z: "//trim(alloc_errmsg)
        errflg = 1
        return
      end if

      call stats_init_sfc_api( clubb_vars_sfc, &
                               l_error, &
                               stats_metadata, stats_sfc(j) )
    end do

    ! Check for errors

    if ( l_error ) then
       errmsg = 'stats_init:  errors found'
       errflg = 1
       return
    endif

  end subroutine stats_init_clubb

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

    !-----------------------------------------------------------------------
  subroutine stats_zero( kk, num_output_fields, x, n, l_in_update )

    !     Description:
    !     Initialize stats to zero
    !-----------------------------------------------------------------------

    use ccpp_kinds, only: kind_phys
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

  end subroutine stats_zero

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine stats_end_timestep_clubb(thecol, stats_zt, stats_zm, stats_rad_zt, stats_rad_zm, stats_sfc, &
                                      out_zt, out_zm, out_radzt, out_radzm, out_sfc, &
                                      stats_metadata, clubb_l_ascending_grid, pver, pverp, top_lev, &
                                      errmsg, errflg )
    !-----------------------------------------------------------------------
    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
    !-----------------------------------------------------------------------

    use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan

    use ccpp_kinds, only: kind_phys

    use clubb_api_module, only: &
        stats, &
        fstderr, & ! Constant(s)
        clubb_at_least_debug_level_api, & ! Procedure(s)
        stats_metadata_type

    implicit none

    integer :: thecol

    integer, intent(in) :: pver, pverp, top_lev

    ! Input Variables
    type (stats), intent(inout) :: stats_zt,      & ! stats_zt grid
                                   stats_zm,      & ! stats_zm grid
                                   stats_rad_zt,  & ! stats_rad_zt grid
                                   stats_rad_zm,  & ! stats_rad_zm grid
                                   stats_sfc        ! stats_sfc

    logical, intent(in) :: clubb_l_ascending_grid

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    ! Inout variables
    real(kind_phys), intent(inout) :: out_zt(:,:,:)
    real(kind_phys), intent(inout) :: out_zm(:,:,:)
    real(kind_phys), intent(inout) :: out_radzt(:,:,:)
    real(kind_phys), intent(inout) :: out_radzm(:,:,:)
    real(kind_phys), intent(inout) :: out_sfc(:,:,:)

    ! Out variables
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Local Variables

    integer :: i, k
    logical :: l_error

    errmsg = ''
    errflg = 0

    !  Check if it is time to write to file

    if ( .not. stats_metadata%l_stats_last ) return

    !  Initialize
    l_error = .false.

    !  Compute averages
    call stats_avg( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, stats_zt%accum_num_samples )
    call stats_avg( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, stats_zm%accum_num_samples )
    if (stats_metadata%l_output_rad_files) then
      call stats_avg( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                      stats_rad_zt%accum_num_samples )
      call stats_avg( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                      stats_rad_zm%accum_num_samples )
    end if
    call stats_avg( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, stats_sfc%accum_num_samples )

   !  Here we are not outputting the data, rather reading the stats into
   !  arrays which are conformable to CAM output.  Also, the data is "flipped"
   !  in the vertical level to be the same as CAM output.
    do i = 1, stats_zt%num_output_fields
      do k = 1, stats_zt%kk

        ! The data stored in stats types are ascending if clubb_l_ascending_grid = .true.
        if ( clubb_l_ascending_grid ) then
          out_zt(thecol,pver+1-k,i) = stats_zt%accum_field_values(1,1,k,i)
        else
          out_zt(thecol,top_lev-1+k,i) = stats_zt%accum_field_values(1,1,k,i)
        end if

        if(is_nan(out_zt(thecol,k,i))) out_zt(thecol,k,i) = 0.0_kind_phys

      enddo
    enddo

    do i = 1, stats_zm%num_output_fields
      do k = 1, stats_zm%kk

        ! The data stored in stats types are ascending if clubb_l_ascending_grid = .true.
        if ( clubb_l_ascending_grid ) then
          out_zm(thecol,pverp+1-k,i) = stats_zm%accum_field_values(1,1,k,i)
        else
          out_zm(thecol,top_lev-1+k,i) = stats_zm%accum_field_values(1,1,k,i)
        end if

        if(is_nan(out_zm(thecol,k,i))) out_zm(thecol,k,i) = 0.0_kind_phys

      enddo
    enddo

    if (stats_metadata%l_output_rad_files) then
      do i = 1, stats_rad_zt%num_output_fields
        do k = 1, stats_rad_zt%kk

          ! The data stored in stats types are ascending if clubb_l_ascending_grid = .true.
          if ( clubb_l_ascending_grid ) then
            out_radzt(thecol,pver+1-k,i) = stats_rad_zt%accum_field_values(1,1,k,i)
          else
            out_radzt(thecol,top_lev-1+k,i) = stats_rad_zt%accum_field_values(1,1,k,i)
          end if

          if(is_nan(out_radzt(thecol,k,i))) out_radzt(thecol,k,i) = 0.0_kind_phys

        enddo
      enddo

      do i = 1, stats_rad_zm%num_output_fields
        do k = 1, stats_rad_zm%kk

          ! The data stored in stats types are ascending if clubb_l_ascending_grid = .true.
          if ( clubb_l_ascending_grid ) then
            out_radzm(thecol,pverp+1-k,i) = stats_rad_zm%accum_field_values(1,1,k,i)
          else
            out_radzm(thecol,top_lev-1+k,i) = stats_rad_zm%accum_field_values(1,1,k,i)
          end if

          if(is_nan(out_radzm(thecol,k,i))) out_radzm(thecol,k,i) = 0.0_kind_phys

        enddo
      enddo

      ! Fill in values above the CLUBB top.
      out_zt(thecol,:top_lev-1,:) = 0.0_kind_phys
      out_zm(thecol,:top_lev-1,:) = 0.0_kind_phys
      out_radzt(thecol,:top_lev-1,:) = 0.0_kind_phys
      out_radzm(thecol,:top_lev-1,:) = 0.0_kind_phys

    endif ! l_output_rad_files

    do i = 1, stats_sfc%num_output_fields
      out_sfc(thecol,1,i) = stats_sfc%accum_field_values(1,1,1,i)
      if(is_nan(out_sfc(thecol,1,i))) out_sfc(thecol,1,i) = 0.0_kind_phys
    enddo

    !  Reset sample fields
    call stats_zero( stats_zt%kk, stats_zt%num_output_fields, stats_zt%accum_field_values, &
                     stats_zt%accum_num_samples, stats_zt%l_in_update )
    call stats_zero( stats_zm%kk, stats_zm%num_output_fields, stats_zm%accum_field_values, &
                     stats_zm%accum_num_samples, stats_zm%l_in_update )
    if (stats_metadata%l_output_rad_files) then
      call stats_zero( stats_rad_zt%kk, stats_rad_zt%num_output_fields, stats_rad_zt%accum_field_values, &
                       stats_rad_zt%accum_num_samples, stats_rad_zt%l_in_update )
      call stats_zero( stats_rad_zm%kk, stats_rad_zm%num_output_fields, stats_rad_zm%accum_field_values, &
                       stats_rad_zm%accum_num_samples, stats_rad_zm%l_in_update )
    end if
    call stats_zero( stats_sfc%kk, stats_sfc%num_output_fields, stats_sfc%accum_field_values, &
                     stats_sfc%accum_num_samples, stats_sfc%l_in_update )

  end subroutine stats_end_timestep_clubb

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

    !-----------------------------------------------------------------------
  subroutine stats_avg( kk, num_output_fields, x, n )

    !     Description:
    !     Compute the average of stats fields
    !-----------------------------------------------------------------------
    use clubb_api_module, only: &
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    !  Input
    integer, intent(in) :: num_output_fields, kk
    integer(kind=stat_nknd), dimension(1,1,kk,num_output_fields), intent(in) :: n

    !  Output
    real(kind=stat_rknd), dimension(1,1,kk,num_output_fields), intent(inout)  :: x

    !  Internal

    integer k,m

    !  Compute averages

    do m = 1, num_output_fields
       do k = 1, kk

          if ( n(1,1,k,m) > 0 ) then
             x(1,1,k,m) = x(1,1,k,m) / real( n(1,1,k,m) )
          end if

       end do
    end do

  end subroutine stats_avg

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  !
  ! DISCLAIMER : this code appears to be correct but has not been
  !              very thouroughly tested. If you do notice any
  !              anomalous behaviour then please contact Andy and/or
  !              Bjorn
  !
  ! Function diag_ustar:  returns value of ustar using the below
  ! similarity functions and a specified buoyancy flux (bflx) given in
  ! kinematic units
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !

  function diag_ustar( z, bflx, wnd, z0, shr_const_karman, &
                       shr_const_pi, shr_const_g )

    use ccpp_kinds, only: kind_phys
!    use shr_const_mod, only : shr_const_karman, shr_const_pi, shr_const_g

    implicit none

    real(kind_phys) :: diag_ustar

    real(kind_phys), parameter      :: am   =  4.8_kind_phys   !   "          "         "
    real(kind_phys), parameter      :: bm   = 19.3_kind_phys  !   "          "         "

    real(kind_phys) :: grav
    real(kind_phys) :: vonk
    real(kind_phys) :: pi

    real(kind_phys), intent (in)    :: z             ! height where u locates
    real(kind_phys), intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
    real(kind_phys), intent (in)    :: wnd           ! wind speed at z
    real(kind_phys), intent (in)    :: z0            ! momentum roughness height

    real(kind_phys), intent (in)    :: shr_const_karman, shr_const_pi, shr_const_g

    integer :: iterate
    real(kind_phys)    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

    grav = shr_const_g
    vonk = shr_const_karman
    pi = shr_const_pi

    lnz   = log( z / z0 )
    klnz  = vonk/lnz
    c1    = pi / 2.0_kind_phys - 3.0_kind_phys*log( 2.0_kind_phys )

    ustar =  wnd*klnz
    if (abs(bflx) > 1.e-6_kind_phys) then
      do iterate = 1, 4

          if (ustar > 1.e-6_kind_phys) then
            lmo   = -ustar**3 / ( vonk * bflx )
            zeta  = z/lmo
            if (zeta > 0._kind_phys) then
                ustar =  vonk*wnd  /(lnz + am*zeta)
            else
                x     = sqrt( sqrt( 1.0_kind_phys - bm*zeta ) )
                psi1  = 2._kind_phys*log( 1.0_kind_phys+x ) + log( 1.0_kind_phys+x*x ) - 2._kind_phys*atan( x ) + c1
                ustar = wnd*vonk/(lnz - psi1)
            end if

          endif

      end do
    end if

    diag_ustar = ustar

  end function diag_ustar

  ! clamped linear-in-pressure interpolation of a midpoint
  ! field to a target pressure level, mirroring the CLUBB-core pvertinterp
  ! (advance_helper_module): below the lowest midpoint or above the highest
  ! midpoint the boundary value is used.  Kept verbatim from the first
  ! implementation for bit-for-bit reproducibility of the LTS trigger.
  subroutine mf_pinterp(ncol, pver, pmid, pout, fld, outv)
    use ccpp_kinds, only: kind_phys

    integer,  intent(in)  :: ncol, pver
    real(kind_phys), intent(in)  :: pmid(:,:)   ! (ncol,pver) midpoint pressure (Pa), increasing downward
    real(kind_phys), intent(in)  :: pout        ! target pressure (Pa)
    real(kind_phys), intent(in)  :: fld(:,:)    ! (ncol,pver) field on midpoints
    real(kind_phys), intent(out) :: outv(ncol)

    integer  :: i, k
    real(kind_phys) :: dpu, dpl

    do i = 1, ncol
      if (pout >= pmid(i,pver)) then
        outv(i) = fld(i,pver)
      else if (pout <= pmid(i,1)) then
        outv(i) = fld(i,1)
      else
        do k = 1, pver-1
          if (pout > pmid(i,k) .and. pout <= pmid(i,k+1)) then
            dpu = pout - pmid(i,k)
            dpl = pmid(i,k+1) - pout
            outv(i) = (fld(i,k)*dpl + fld(i,k+1)*dpu) / (dpl + dpu)
            exit
          end if
        end do
      end if
    end do

  end subroutine mf_pinterp

end module clubb

