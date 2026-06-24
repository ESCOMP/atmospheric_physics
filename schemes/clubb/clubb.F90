module clubb

  use ccpp_kinds, only: kind_phys

  implicit none

  private

  save

  ! Subroutines to make public
  public :: clubb_init, clubb1_run, clubb2_run, clubb3_run, stats_zero

  contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine clubb_init(pcols, pver, pverp, pcnst, & ! in
                        masterproc, mpicom, mstrid, mpi_character, & ! in
                        iulog, max_fieldname_len, & ! in
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
    integer, intent(in) :: mpicom, mpi_character, mstrid
    integer, intent(in) :: iulog, pcnst, max_fieldname_len
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
    character(len=200), intent(out) :: errmsg


    ! Local variables
    real(kind=time_precision) :: dum1, dum2, dum3

    type(err_info_type) :: &
      err_info          ! err_info struct used in CLUBB containing err_code and err_header

    integer :: j, m
    integer :: ierr = 0


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
       errmsg = 'clubb_ini_cam: FATAL ERROR CALLING CHECK_CLUBB_SETTINGS_API'
       errflg = 1
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
         errmsg = 'clubb_ini_cam: Unable to allocate hm_metadata%hydromet_list'
         errflg = 1
         return
       end if
       allocate( hm_metadata%l_mix_rat_hm(1), stat=ierr)
       if( ierr /= 0 ) then
         errmsg = 'clubb_ini_cam: Unable to allocate hm_metadata%l_mix_rat_hm'
         errflg = 1
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
                             errmsg, errflg)

      if (errflg /= 0) return 

       allocate(out_zt(pcols,pver,stats_zt(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) then 
         errmsg = 'clubb_ini_cam: Unable to allocate out_zt'
         errflg = 1
         return
       end if
       allocate(out_zm(pcols,pverp,stats_zm(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) then 
         errmsg = 'clubb_ini_cam: Unable to allocate out_zm'
         errflg = 1
         return
       end if
       allocate(out_sfc(pcols,1,stats_sfc(1)%num_output_fields), stat=ierr)
       if( ierr /= 0 ) then
         errmsg = 'clubb_ini_cam: Unable to allocate out_sfc'
         errflg = 1
         return
       end if

       if ( stats_metadata%l_output_rad_files ) then
          allocate(out_radzt(pcols,pver,stats_rad_zt(1)%num_output_fields), stat=ierr)
          if( ierr /= 0 ) then
            errmsg = 'clubb_ini_cam: Unable to allocate out_radzt'
            errflg = 1
            return
          end if
          allocate(out_radzm(pcols,pverp,stats_rad_zm(1)%num_output_fields), stat=ierr)
          if( ierr /= 0 ) then
            errmsg = 'clubb_ini_cam: Unable to allocate out_radzm'
            errflg = 1
            return
          end if
       end if

    endif

    end subroutine clubb_init

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
                        rga, gravit, clubb_rnevap_effic, do_cldcool, do_rainturb, & ! in
                        do_clubb_mf, l_implemented, grid_type, lq, deep_scheme, & ! in
                        state_q, t, pmid, zm, & ! in
                        phis, pdel, pdeldry, & ! in
                        pint, zi, omega, wsx, & ! in
                        wsy, cflx, shf, landfrac, & ! in
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
                        s_ae, s_aw, s_awthl, s_awqt, s_awql, s_awqi, s_awu, s_awv, & ! inout
                        mf_thlflx, mf_qtflx, & ! inout
                        thlm, rtm, um, vm, wm_zt, rcm, rcm_in_layer, & ! inout
                        wp2_pbuf, wp3_pbuf, wpthlp_pbuf, wprtp_pbuf, & ! inout
                        rtpthlp_pbuf, rtp2_pbuf, thlp2_pbuf, rtp3_pbuf, & ! inout
                        thlp3_pbuf, up2_pbuf, vp2_pbuf, up3_pbuf, vp3_pbuf, & ! inout
                        upwp_pbuf, vpwp_pbuf, wpthvp_pbuf, wp2thvp_pbuf, wp2up_pbuf, & ! inout
                        rtpthvp_pbuf, thlpthvp_pbuf, pdf_zm_w_1_pbuf, pdf_zm_w_2_pbuf, & ! inout
                        pdf_zm_varnce_w_1_pbuf, pdf_zm_varnce_w_2_pbuf, pdf_zm_mixt_frac_pbuf, & ! inout
                        wp2rtp_pbuf, wp2thlp_pbuf, uprcp_pbuf, vprcp_pbuf, rc_coef_zm_pbuf, & ! inout
                        wp4_pbuf, wpup2_pbuf, wpvp2_pbuf, wp2up2_pbuf, wp2vp2_pbuf, cld_pbuf, & ! inout
                        concld_pbuf, ast_pbuf, alst_pbuf, aist_pbuf, qlst_pbuf, qist_pbuf, & ! inout
                        deepcu_pbuf, shalcu_pbuf, khzm_pbuf, pblh_pbuf, tke_pbuf, dp_icwmr_pbuf, & ! inout
                        ice_supersat_frac_pbuf, relvar_pbuf, naai_pbuf, cmeliq_pbuf, & ! inout
                        cmfmc_sh_pbuf, qsatfac_pbuf, npccn_pbuf, prer_evap_pbuf, qrl_pbuf, & ! inout
                        rtp2_mc_zt_pbuf, thlp2_mc_zt_pbuf, wprtp_mc_zt_pbuf, & ! inout
                        wpthlp_mc_zt_pbuf, rtpthlp_mc_zt_pbuf, ttend_clubb_pbuf, & ! inout
                        upwp_clubb_gw_pbuf, vpwp_clubb_gw_pbuf, thlp2_clubb_gw_pbuf, & ! inout
                        wpthlp_clubb_gw_pbuf, ttend_clubb_mc_pbuf, upwp_clubb_gw_mc_pbuf, & ! inout
                        vpwp_clubb_gw_mc_pbuf, thlp2_clubb_gw_mc_pbuf, wpthlp_clubb_gw_mc_pbuf, & ! inout
                        stats_zt, stats_zm, stats_sfc, stats_rad_zt, stats_rad_zm, & ! inout
                        out_zt, out_zm, out_sfc, out_radzt, out_radzm, & ! inout
                        invrs_cpairv, clubbtop_pbuf, & ! inout
                        errmsg, errflg ) ! out

    use clubb_mf,                  only: integrate_mf
    use atmos_phys_pbl_utils,  only: calc_friction_velocity, calc_ideal_gas_rrho
    use clubb_api_module, only: &
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
      implicit_coefs_terms


    ! Input variables, intent(in)
    integer, intent(in) :: ncol, iam, ixq, ixcldliq, ixcldice, ixrtpthlp, ixwpthlp, &
                           ixwprtp, ixwp3, ixwp2, ixthlp2, ixrtp2, ixup2, ixvp2, sclr_dim, edsclr_dim, &
                           macmic_it, top_lev, cld_macmic_num_steps, grid_type, hydromet_dim, nstep, &
                           nzt_clubb, nzm_clubb, pver, pverp, pcnst
    real(kind_phys), intent(in) :: hdtime, clubb_timestep
    real(kind_phys), intent(in) :: sclr_tol(:)
    real(kind_phys), intent(in) :: shr_const_karman, shr_const_pi, shr_const_g, omega_const, cpair, &
                                   rga, inv_p0_clubb, rair, zvir, latvap, latice, gravit, theta0, &
                                   ts_nudge, rtm_min, rtm_nudge_max_altitude, clubb_rnevap_effic
    real(kind_phys), intent(in) :: cpairv(:,:), rairv(:,:)
    real(kind_phys), intent(in) :: rtpthlp_const, wpthlp_const, wprtp_const, wp3_const
    real(kind_phys), intent(in) :: state_q(:,:,:)
    real(kind_phys), intent(in) :: wsx(:), wsy(:), shf(:)
    real(kind_phys), intent(in) :: cflx(:,:)
    real(kind_phys), intent(in) :: clubb_params_single_col(:,:)
    real(kind_phys), intent(in) :: lat(:), lon(:), phis(:)
    real(kind_phys), intent(in) :: pint(:,:)
    real(kind_phys), intent(in) :: pmid(:,:)
    real(kind_phys), intent(in) :: landfrac(:)
    real(kind_phys), intent(in) :: pdel(:,:), pdeldry(:,:), omega(:,:), &
                                   t(:,:), zm(:,:), zi(:,:)

    logical, intent(in) :: clubb_do_adv, first_step, first_restart_step, l_implemented, &
                           do_clubb_mf, single_column, scm_cambfb_mode, clubb_l_intr_sfc_flux_smooth, &
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
    real(kind_phys), intent(inout) :: concld_pbuf(:,:)
    real(kind_phys), intent(inout) :: ast_pbuf(:,:)
    real(kind_phys), intent(inout) :: alst_pbuf(:,:)
    real(kind_phys), intent(inout) :: aist_pbuf(:,:)
    real(kind_phys), intent(inout) :: qlst_pbuf(:,:)
    real(kind_phys), intent(inout) :: qist_pbuf(:,:)
    real(kind_phys), intent(inout) :: deepcu_pbuf(:,:)
    real(kind_phys), intent(inout) :: shalcu_pbuf(:,:)
    real(kind_phys), intent(inout) :: khzm_pbuf(:,:)
    real(kind_phys), intent(inout) :: pblh_pbuf(:)
    real(kind_phys), intent(inout) :: tke_pbuf(:,:)
    real(kind_phys), intent(inout) :: dp_icwmr_pbuf(:,:)
    real(kind_phys), intent(inout) :: ice_supersat_frac_pbuf(:,:)
    real(kind_phys), intent(inout) :: relvar_pbuf(:,:)
    real(kind_phys), intent(inout) :: naai_pbuf(:,:)
    real(kind_phys), intent(inout) :: cmeliq_pbuf(:,:)
    real(kind_phys), intent(inout) :: cmfmc_sh_pbuf(:,:)

    real(kind_phys), intent(inout) :: qsatfac_pbuf(:,:)
    real(kind_phys), intent(inout) :: npccn_pbuf(:,:)
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
      s_awql(:,:),     s_awqi(:,:),        &
      s_awu(:,:),      s_awv(:,:),         &
      mf_thlflx(:,:),  mf_qtflx(:,:)

    real(kind_phys), intent(inout) :: apply_const

    integer, intent(inout) :: clubbtop_pbuf(:)

    real(kind_phys), intent(inout) :: ztodtptr

    type(grid), intent(inout) :: &
      gr

    ! Input variables, intent(out)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! ------------- Local variables --------------
    real(kind_phys) :: rad2deg

    real(kind_phys) :: &
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
      ! MF local thermodynamic vars
      invrs_exner_zt,                 & ! thermodynamic grid
      kappa_zt                          ! thermodynamic grid

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
      rtpthlp_mc

    real(kind_phys), dimension(ncol,nzm_clubb) :: &
      ! MF local momentum vars
      rtm_zm,     thlm_zm,       & ! momentum grid
      kappa_zm,   p_in_Pa_zm,    & ! momentum grid
                  invrs_exner_zm   ! momentum grid

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
      rrho_tmp,                 &
      ! Variables below are needed to compute energy integrals for conservation
      te_a, se_a, ke_a, wv_a, wl_a, &
      te_b, se_b, ke_b, wv_b, wl_b

    intrinsic :: max

    real(kind_phys), dimension(ncol,nparams) :: &
      clubb_params

    character(len=200) :: temp1, sub             ! Strings needed for CLUBB output

    integer :: &
      i, j, k, tt, ixind, nadv, n,      & ! Loop variables
      k_cam, k_clubb, sclr, iedsclr, & ! Loop variables
      icnt, &
      stats_nsamp, stats_nout         ! Stats sampling and output intervals for CLUBB [timestep]

    ! -------------------------------------------------
    ! ----------------- Start -------------------------

    rad2deg = 180.0_kind_phys / shr_const_pi

    ! Initialize err_info with parallelization and geographical info
    ! BAS putting a -1 here in place of "lchnk"
    call init_err_info_api(ncol, -1, iam, lat*rad2deg, lon*rad2deg, err_info)

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
    !$acc              state_q, u, v, t, pmid, &
    !$acc              zm, phis, pdel, pdeldry, s, &
    !$acc              pint, zi, omega, lat, &
    !$acc              wsx, wsy, cflx, shf, &
    !$acc              err_info, err_info%err_header, &
    !$acc              cpairv, rairv, se_dis, eleak, cld_pbuf, clubb_params_single_col, grid_dx, grid_dy ) &
    !$acc     copyout( clubb_s, clubbtop_pbuf, &
    !$acc              qclvar, wprcp, rcm_in_layer, rcm, cloud_frac, thlm, rtm, &
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
    !call t_stopf('clubb_tend_cam:acc_copyin')
    !call t_startf('clubb_tend_cam:acc_region')

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

    !----------------------------------------- Initializing arrays -----------------------------------------

    if ( clubb_do_adv ) then

      if (macmic_it  ==  1) then
        
        !  Note that some of the moments below can be positive or negative.
        !    Remove a constant that was added to prevent dynamics from clipping
        !    them to prevent dynamics from making them positive.
        do k = 1, nzm_clubb
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
        zt_g(i,k) = zm(i,k_cam) - zi(i,pverp)

        invrs_dz_g(i,k) = 1._kind_phys / ( zi(i,k_cam) - zi(i,k_cam+1) )  ! compute thickness

        rho_zt(i,k)          = rga * pdel(i,k_cam)    * invrs_dz_g(i,k)

        rho_ds_zt(i,k)       = rga * pdeldry(i,k_cam) * invrs_dz_g(i,k)

        invrs_rho_ds_zt(i,k) = 1._kind_phys / rho_ds_zt(i,k)

      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol

        k_cam = top_lev - 1 + k  

        p_in_Pa(i,k) = pmid(i,k_cam)
        
        !  Compute inverse exner function consistent with CLUBB's definition, which uses a constant
        !  surface pressure.  CAM's exner (in state) does not.  Therefore, for consistent
        !  treatment with CLUBB code, anytime exner is needed to treat CLUBB variables
        !  (such as thlm), use "invrs_exner_zt" otherwise use the exner in state
        exner(i,k) = ( p_in_Pa(i,k) * inv_p0_clubb )**( rairv(i,k_cam) * invrs_cpairv(i,k_cam) )

        invrs_exner_zt(i,k) = 1._kind_phys / exner(i,k)

        ! exception - setting this to moist thv_ds_zt
        thv_ds_zt(i,k) = t(i,k_cam) * invrs_exner_zt(i,k)  &
                         * (1._kind_phys + zvir * state_q(i,k_cam,ixq) - state_q(i,k_cam,ixcldliq))

        rcm(i,k)    = state_q(i,k_cam,ixcldliq)
        rtm(i,k)    = state_q(i,k_cam,ixq) + state_q(i,k_cam,ixcldliq)

        thlm(i,k)   = ( t(i,k_cam) - ( latvap * invrs_cpairv(i,k_cam) ) &
                                               * state_q(i,k_cam,ixcldliq) ) * invrs_exner_zt(i,k)
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol

        k_cam = top_lev - 1 + k  

        !  Compute mean w wind on thermo grid, convert from omega to w
        wm_zt(i,k) = -1._kind_phys * ( omega(i,k_cam) - omega(i,pver) ) / ( rho_zt(i,k) * gravit )

        cloud_frac(i,k)       = cld_pbuf(i,k_cam)

        um(i,k) = u(i,k_cam)
        vm(i,k) = v(i,k_cam)

        rfrzm(i,k)  = state_q(i,k_cam,ixcldice)

      end do
    end do

    !$acc parallel loop gang vector default(present)
    do i = 1, ncol

      deltaz(i)         = zi(i,pverp-1) - zi(i,pverp)

      !  Set the surface pressure      
      p_sfc(i)          = pint(i,pverp)

      !  Set the elevation of the surface
      sfc_elevation(i)  = zi(i,pverp)

    end do

    ! Define the CLUBB momentum grid (in height, units of m)
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm_clubb
      do i = 1, ncol
        k_cam = top_lev - 1 + k  
        zi_g(i,k) = zi(i,k_cam) - zi(i,pverp)
      end do
    end do

    if (do_clubb_mf) then

      do k = 1, nzt_clubb
        do i = 1, ncol
          k_cam = top_lev - 1 + k
          kappa_zt(i,k)   = rairv(i,k_cam) * invrs_cpairv(i,k_cam)
          dz_g(i,k)       = zi(i,k_cam) - zi(i,k_cam+1)  ! compute thickness
        end do
      end do

      ! pressure on momentum grid needed for mass flux calc.
      do k = 1, nzm_clubb
        do i = 1, ncol
          k_cam = top_lev - 1 + k
          p_in_Pa_zm(i,k)     = pint(i,k_cam)
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
    !call t_stopf('clubb_tend_cam:acc_region')
    !call t_startf('clubb_tend_cam:non_acc_region')
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
       errflg = 0
       return
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
       errflg = 0
       return
    end if

    ! CLUBB's grid data structure (gr) and nu_vert_res_dep contain arrays that need to
    ! be copied to the GPU
!    !call t_stopf('clubb_tend_cam:non_acc_region')
!    !call t_startf('clubb_tend_cam:acc_copyin')
    !$acc data copyin( gr, gr%zm, gr%zt, gr%dzm, gr%dzt, gr%invrs_dzt, gr%invrs_dzm, &
    !$acc              gr%weights_zt2zm, gr%weights_zm2zt, &
    !$acc              nu_vert_res_dep, nu_vert_res_dep%nu2, nu_vert_res_dep%nu9, &
    !$acc              nu_vert_res_dep%nu1, nu_vert_res_dep%nu8, nu_vert_res_dep%nu10, &
    !$acc              nu_vert_res_dep%nu6)
    !call t_stopf('clubb_tend_cam:acc_copyin')
    !call t_startf('clubb_tend_cam:acc_region')
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

      !call t_stopf('clubb_tend_cam:acc_region')
      !call t_startf('clubb_tend_cam:non_acc_region')
      !$acc update host( u, v, t, pmid, wsx, wsy )

      ! Adjust surface stresses using winds from the prior macmic iteration
      do i = 1, ncol
        ubar = sqrt(u(i,pver)**2+v(i,pver)**2)
        if (ubar <  0.25_kind_phys) ubar = 0.25_kind_phys

        rrho_tmp = calc_ideal_gas_rrho(rair, t(i,pver), pmid(i,pver))
        ustar    = calc_friction_velocity(wsx(i), wsy(i), rrho_tmp)

        upwp_sfc(i) = -u(i,pver)*ustar**2/ubar
        vpwp_sfc(i) = -v(i,pver)*ustar**2/ubar
      end do

      !$acc update device( upwp_sfc, vpwp_sfc )
      !call t_stopf('clubb_tend_cam:non_acc_region')
      !call t_startf('clubb_tend_cam:acc_region')

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

      !#######################################################################
      !###################### CALL MF DIAGNOSTIC PLUMES ######################
      !#######################################################################
      if (do_clubb_mf) then
        !call t_startf('clubb_tend_cam:do_clubb_mf')

        rtm_zm     = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr,  rtm(:ncol,:) )
        thlm_zm    = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, thlm(:ncol,:) )

        ! exner on momentum grid needed for mass flux calc.
        kappa_zm = zt2zm_api( nzm_clubb, nzt_clubb, ncol, gr, kappa_zt )

        do k = 1, nzm_clubb
          do i = 1, ncol
            k_cam = top_lev - 1 + k
            invrs_exner_zm(i,k) = 1._kind_phys / ( (p_in_Pa_zm(i,k) * inv_p0_clubb)**kappa_zm(i,k) )
          end do
        end do

        !--------------------------------------- integrate_mf call ---------------------------------------
        ! integrate_mf expects arguments of individual columns.
        ! If the column loop gets pushed into it, we can also avoid the array slicing.

        do i = 1, ncol
          call integrate_mf( nzm_clubb, nzt_clubb, dz_g(i,:), zi_g(i,:), p_in_Pa_zm(i,:), invrs_exner_zm(i,:),  & ! input
                                                                            p_in_Pa(i,:), invrs_exner_zt(i,:),  & ! input
                            um(i,:), vm(i,:), thlm(i,:),        rtm(i,:), thv_ds_zt(i,:),                       & ! input
                                                            thlm_zm(i,:),    rtm_zm(i,:),                       & ! input
                                                            wpthlp_sfc(i),  wprtp_sfc(i),  pblh_pbuf(i),        & ! input
                            mf_dry_a(i,:),    mf_moist_a(i,:),                                        & ! output - plume diagnostics
                            mf_dry_w(i,:),    mf_moist_w(i,:),                                        & ! output - plume diagnostics
                            mf_dry_qt(i,:),   mf_moist_qt(i,:),                                       & ! output - plume diagnostics
                            mf_dry_thl(i,:),  mf_moist_thl(i,:),                                      & ! output - plume diagnostics
                            mf_dry_u(i,:),    mf_moist_u(i,:),                                        & ! output - plume diagnostics
                            mf_dry_v(i,:),    mf_moist_v(i,:),                                        & ! output - plume diagnostics
                                              mf_moist_qc(i,:),                                       & ! output - plume diagnostics
                            s_ae(i,:),        s_aw(i,:),                                              & ! output - plume diagnostics
                            s_awthl(i,:),     s_awqt(i,:),                                            & ! output - plume diagnostics
                            s_awql(i,:),      s_awqi(i,:),                                            & ! output - plume diagnostics
                            s_awu(i,:),       s_awv(i,:),                                             & ! output - plume diagnostics
                            mf_thlflx(i,:),   mf_qtflx(i,:) )                                 ! output - variables needed for solver
        end do

        !--------------------------------------- END integrate_mf call ---------------------------------------

        ! pass MF turbulent advection term as CLUBB explicit forcing term
        do k = 1, nzt_clubb
          do i = 1, ncol
            rtm_forcing(i,k)  = rtm_forcing(i,k) - invrs_rho_ds_zt(i,k) * invrs_dz_g(i,k) * &
                              ((rho_ds_zm(i,k) * mf_qtflx(i,k)) - (rho_ds_zm(i,k+1) * mf_qtflx(i,k+1)))

            thlm_forcing(i,k) = thlm_forcing(i,k) - invrs_rho_ds_zt(i,k) * invrs_dz_g(i,k) * &
                               ((rho_ds_zm(i,k) * mf_thlflx(i,k)) - (rho_ds_zm(i,k+1) * mf_thlflx(i,k+1)))
          end do
        end do
        !call t_stopf('clubb_tend_cam:do_clubb_mf')

      end if

      
      if ( clubb_l_ascending_grid ) then

        ! CLUBB is to be run in ascending mode, which has the surface at k=1, which is 
        ! the opposite of the cam grid that the rest of clubb_intr uses, so
        ! we need to flip the fields (in the vertical dimensions) before calling advance_clubb_core
        !
        ! NOTE: We do not neccesarily flip all arrays, only ones that are used within this
        !       subroutine (advance_clubb_core). For example, only the pdf_params fields that 
        !       are used within this subroutine (or used in a subroutine we call) need to
        !       be flipped. 
        
        !call t_startf('clubb_tend_cam:ascending_grid_flip')

        thlm_forcing              =              thlm_forcing(:,nzt_clubb:1:-1)
        rtm_forcing               =               rtm_forcing(:,nzt_clubb:1:-1)
        um_forcing                =                um_forcing(:,nzt_clubb:1:-1)
        vm_forcing                =                vm_forcing(:,nzt_clubb:1:-1)
        wm_zt                     =                     wm_zt(:,nzt_clubb:1:-1)
        rho_zt                    =                    rho_zt(:,nzt_clubb:1:-1)
        rho_ds_zt                 =                 rho_ds_zt(:,nzt_clubb:1:-1)
        invrs_rho_ds_zt           =           invrs_rho_ds_zt(:,nzt_clubb:1:-1)
        thv_ds_zt                 =                 thv_ds_zt(:,nzt_clubb:1:-1)
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
                             zi_g(:,nzm_clubb:1:-1), zt_g(:,nzt_clubb:1:-1),   & ! intent(in)
                             gr, err_info )                                      ! intent(inout)
  
        !call t_stopf('clubb_tend_cam:ascending_grid_flip')

      end if

      !  Advance CLUBB CORE one timestep in the future
      !call t_startf('clubb_tend_cam:advance_clubb_core_api')

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

      !call t_stopf('clubb_tend_cam:advance_clubb_core_api')
          

      if ( clubb_l_ascending_grid ) then

        !call t_startf('clubb_tend_cam:ascending_grid_flip')

        ! If running in ascending mode, we flip the arrays before calling advance_clubb_core
        ! so we need to flip them back. This section should flip every array that was flipped 
        ! before the advance_clubb_core call.

        thlm_forcing               =               thlm_forcing(:,nzt_clubb:1:-1)
        rtm_forcing                =                rtm_forcing(:,nzt_clubb:1:-1)
        um_forcing                 =                 um_forcing(:,nzt_clubb:1:-1)
        vm_forcing                 =                 vm_forcing(:,nzt_clubb:1:-1)
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

        !call t_stopf('clubb_tend_cam:ascending_grid_flip')

      end if

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) "Fatal error in CLUBB advance_clubb_core: at timestep ", nstep
        errmsg = 'clubb1_run: '//err_info%err_header_global//NEW_LINE('a')//'Fatal error in CLUBB advance_clubb_core'
        errflg = 0
        return
      end if

      if ( do_rainturb ) then

!        !call t_startf('clubb_tend_cam:do_rainturb')

        do k = 1, nzt_clubb
          do i = 1, ncol
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

        !call t_stopf('clubb_tend_cam:do_rainturb')

      end if

      if (do_cldcool) then

        !call t_startf('clubb_tend_cam:do_cldcool')

        thlp2_rad(:,:) = 0._kind_phys
        
        do k = 1, nzt_clubb
          do i = 1, ncol
            k_cam = top_lev - 1 + k
            qrl_clubb(i,k) = qrl_pbuf(i,k_cam) / ( cpairv(i,k_cam) * pdeldry(i,k_cam) )
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

        !call t_stopf('clubb_tend_cam:do_cldcool')

      end if

      !  Check to see if stats should be output, here stats are read into
      !  output arrays to make them conformable to CAM output
      if (stats_metadata%l_stats) then
        !call t_startf('clubb_tend_cam:stats_end_timestep_clubb')
        do i = 1, ncol
          call stats_end_timestep_clubb(i, stats_zt(i), stats_zm(i), stats_rad_zt(i), stats_rad_zm(i), stats_sfc(i), &
                                        out_zt, out_zm, out_radzt, out_radzm, out_sfc, &
                                        stats_metadata, clubb_l_ascending_grid, pver, pverp, top_lev, &
                                        errmsg, errflg )
        end do
        !call t_stopf('clubb_tend_cam:stats_end_timestep_clubb')
      end if

    end do  ! end time loop
    !----------------------------------------- END substepping loop -----------------------------------------


    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt_clubb
      do i = 1, ncol
        k_cam = top_lev - 1 + k
        qclvar(i,k)        = min( 1._kind_phys, qclvar(i,k) ) ! We should move this clipping inside clubb
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm_clubb
      do i = 1, ncol
        k_cam = top_lev - 1 + k
        khzm_pbuf(i,k_cam)         = khzm(i,k)
      end do
    end do

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
                       + gravit * zm(i,k) + phis(i)
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
        se_a = se_a + clubb_s(i,k_clubb)*pdel(i,k)*rga
        ke_a = ke_a + 0.5_kind_phys*(um(i,k_clubb)**2+vm(i,k_clubb)**2)*pdel(i,k)*rga
        wv_a = wv_a + (rtm(i,k_clubb)-rcm(i,k_clubb))*pdeldry(i,k)*rga
        wl_a = wl_a + (rcm(i,k_clubb))*pdeldry(i,k)*rga
      end do

      ! Based on these integrals, compute the total energy after CLUBB call
      te_a = se_a + ke_a + (latvap+latice) * wv_a + latice * wl_a

      do k = top_lev, pver
        ! Do the same as above, but for before CLUBB was called.
        se_b = se_b + s(i,k)*pdel(i,k)*rga
        ke_b = ke_b + 0.5_kind_phys*(u(i,k)**2+v(i,k)**2)*pdel(i,k)*rga
        wv_b = wv_b + state_q(i,k,ixq)*pdeldry(i,k)*rga
        wl_b = wl_b + state_q(i,k,ixcldliq)*pdeldry(i,k)*rga
      end do

      ! Based on these integrals, compute the total energy before CLUBB call
      te_b = se_b + ke_b + (latvap+latice) * wv_b + latice * wl_b

      ! Take into account the surface fluxes of heat and moisture
      !  Use correct qflux from cam_in, not lhf/latvap as was done previously
      te_b = te_b + (shf(i)+cflx(i,1)*(latvap+latice)) * hdtime

      ! Compute the disbalance of total energy, over depth where CLUBB is active
      se_dis(i) = ( te_a - te_b ) / ( pint(i,pverp) - pint(i,clubbtop_pbuf(i)) )

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

    !call t_stopf('clubb_tend_cam:acc_region')

    !call t_startf('clubb_tend_cam:acc_copyout')
    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data
    !call t_stopf('clubb_tend_cam:acc_copyout')

    !call t_startf('clubb_tend_cam:non_acc_region')

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

  end subroutine clubb1_run

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

  subroutine clubb2_run(ncol, pver, ixcldliq, ixcldice, ixnumliq, ixnumice, & ! in
                        clubb_detliq_rad, clubb_detice_rad, clubb_detphase_lowtemp, &! in
                        meltpt_temp, latice, rga, & ! in
                        dlf, t, pdel, pdeldry, & ! in
                        q, s, det_s, det_ice, & ! inout
                        dlf_liq_out, dlf_ice_out ) ! out

    ! Input variables, intent(in)
    integer, intent(in) :: ncol, pver, ixcldliq, ixcldice, ixnumliq, ixnumice
    real(kind_phys), intent(in) :: clubb_detliq_rad, clubb_detice_rad, clubb_detphase_lowtemp
    real(kind_phys), intent(in) :: meltpt_temp, latice, rga
    real(kind_phys), intent(in) :: dlf(:,:), t(:,:), pdel(:,:), pdeldry(:,:)

    ! Input variables, intent(inout)
    real(kind_phys), intent(inout) :: q(:,:,:)
    real(kind_phys), intent(inout) :: s(:,:)
    real(kind_phys), intent(inout) :: det_s(:), det_ice(:)
 
    ! Input variables, intent(out)
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

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

  subroutine clubb3_run(ncol, pver, pverp, pcnst, top_lev, & ! in
                        ixq, ixcldice, ixcldliq, ixnumice, & ! in
                        rhminis_const, rhmaxis_const, rhmini_const, rhmaxi_const, & ! in
                        dp1, dp2, zvir, rair, cpair, gravit, karman, & ! in
                        calday, tropp_days, & ! in
                        lat, phis, landfrac, snowhland, & ! in
                        wsx, wsy, shf, & ! in
                        pint, pmid, pdel, pdeldry, & ! in
                        rcm, cloud_frac, t, exner, & ! in
                        state_exner, zm, zi, u, & ! in
                        v, cmfmc, cflx, state_q, & ! in
                        single_column, scm_cambfb_mode, lq, & ! in
                        cnst_type, scm_clubb_iop_name, subcol_scheme, & ! in
                        pblh_pbuf, alst_pbuf, qlst_pbuf, deepcu_pbuf, shalcu_pbuf, & ! inout
                        cmfmc_sh_pbuf, dp_icwmr_pbuf, concld_pbuf, aist_pbuf,  & ! inout
                        qsatfac_pbuf, ast_pbuf, qist_pbuf, cld_pbuf, ptend_q, troplev, & ! inout
                        errmsg, errflg ) ! out

!BAS prob need to switch to:
!    use tropopause_find,       only: tropopause_find_run
!but I need to calculate tropp_p_loc
    use holtslag_boville_diff, only: hb_pbl_dependent_coefficients_run
!BAS I need to figure out what to do about aist_vector
!    use cldfrc2m,              only: aist_vector
    use atmos_phys_pbl_utils,  only: calc_friction_velocity, calc_obukhov_length, calc_ideal_gas_rrho, &
                                     calc_kinematic_heat_flux, calc_kinematic_water_vapor_flux, &
                                     calc_kinematic_buoyancy_flux

    ! Input variables, intent(in)
    integer, intent(in) :: ncol, pver, pverp, pcnst, top_lev 
    integer, intent(in) :: ixq, ixcldice, ixcldliq, ixnumice
    real(kind_phys), intent(in) :: rhminis_const, rhmaxis_const, rhmini_const, rhmaxi_const
    real(kind_phys), intent(in) :: dp1, dp2, zvir, rair, cpair, gravit, karman
    real(kind_phys), intent(in) :: calday
    real(kind_phys), intent(in) :: tropp_days(:)
    real(kind_phys), intent(in) :: lat(:), phis(:), landfrac(:), snowhland(:), &
                                   wsx(:), wsy(:), shf(:)
    real(kind_phys), intent(in) :: pint(:,:), pmid(:,:), pdel(:,:), pdeldry(:,:), rcm(:,:), &
                                   cloud_frac(:,:), t(:,:), exner(:,:), state_exner(:,:), &
                                   zm(:,:), zi(:,:), u(:,:), v(:,:), cmfmc(:,:), cflx(:,:)
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
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Local variables
    integer :: i, k, ixind, k_clubb
    !BAS troplev is inout right now, but ultimately could be local
    integer, intent(inout) :: troplev(:)
    real(kind_phys) :: rhmini(ncol), rhmaxi(ncol)
    real(kind_phys) :: frac_limit, ic_limit
    real(kind_phys) :: rrho(ncol), ustar2(ncol), kinheat(ncol), kinwat(ncol), kbfs(ncol), obklen(ncol), &
                       dummy2(ncol), dummy3(ncol)
    real(kind_phys) :: th(ncol,pver), thv(ncol,pver)

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

!    if (errflg /= 0) return

!    aist_pbuf(:,:top_lev-1) = 0._kind_phys
!    qsatfac_pbuf(:, :) = 0._kind_phys ! Zero out entire profile in case qsatfac is left undefined in aist_vector below

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
         th(i,k) = t(i,k)*state_exner(i,k)
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

  end subroutine clubb3_run

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

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
                               errmsg, errflg )
    !
    ! Description: Initializes the statistics saving functionality of
    !   the CLUBB model.  This is for purpose of CAM-CLUBB interface.  Here
    !   the traditional stats_init of CLUBB is not called, as it is not compatible
    !   with CAM output.

    !-----------------------------------------------------------------------

    use clubb_api_module, only: time_precision, stats, stats_metadata_type, hm_metadata_type, &   !
                                nvarmax_zm, stats_init_zm_api, & !
                                nvarmax_zt, stats_init_zt_api, & !
                                nvarmax_rad_zt, stats_init_rad_zt_api, & !
                                nvarmax_rad_zm, stats_init_rad_zm_api, & !
                                nvarmax_sfc, stats_init_sfc_api, & !
                                fstderr, var_length !
!BAS addfld here? 
    use cam_history,      only: addfld, horiz_only

    implicit none

    !----------------------- Input Variables -----------------------

    logical, intent(in) :: l_stats_in ! Stats on? T/F
  
    logical, intent(in) :: masterproc

    integer, intent(in) :: mpicom, mpi_character, mstrid, pcols, max_fieldname_len

    integer, intent(in) :: sclr_dim, edsclr_dim, hydromet_dim

    integer, intent(inout) :: errflg
    character(len=200), intent(inout) :: errmsg

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

      allocate( stats_zt(j)%z( stats_zt(j)%kk ), stat=ierr )
      if( ierr /= 0 ) then
      errmsg = "stats_init_clubb: Failed to allocate stats_zt%z"
      errflg = 1
      return
    end if

      allocate( stats_zt(j)%accum_field_values( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%accum_field_values"
        errflg = 1
        return
      end if
      allocate( stats_zt(j)%accum_num_samples( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%accum_num_samples"
        errflg = 1
        return
      end if
      allocate( stats_zt(j)%l_in_update( 1, 1, stats_zt(j)%kk, stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%l_in_update"
        errflg = 1
        return
      end if
      call stats_zero( stats_zt(j)%kk, stats_zt(j)%num_output_fields, stats_zt(j)%accum_field_values, &
                       stats_zt(j)%accum_num_samples, stats_zt(j)%l_in_update )

      allocate( stats_zt(j)%file%grid_avg_var( stats_zt(j)%num_output_fields ), stat=ierr )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%file%grid_avg_var"
        errflg = 1
        return
      end if
      allocate( stats_zt(j)%file%z( stats_zt(j)%kk ), stat=ierr )
      if( ierr /= 0 ) then
        errmsg = "stats_init_clubb: Failed to allocate stats_zt%file%z"
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
            errmsg = "stats_init_clubb:  number of rad_zt statistical variables exceeds limit"
            errflg = 1
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
            errmsg = "stats_init_clubb:  number of rad_zm statistical variables exceeds limit"
            errflg = 1
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
         errmsg = "stats_init_clubb:  number of sfc statistical variables exceeds limit"
         errflg = 1
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
       errmsg = 'stats_init:  errors found'
       errflg = 1
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

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

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

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS
  subroutine stats_end_timestep_clubb(thecol, stats_zt, stats_zm, stats_rad_zt, stats_rad_zm, stats_sfc, &
                                      out_zt, out_zm, out_radzt, out_radzm, out_sfc, &
                                      stats_metadata, clubb_l_ascending_grid, pver, pverp, top_lev, &
                                      errmsg, errflg )
    !-----------------------------------------------------------------------
    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
    !-----------------------------------------------------------------------

    use shr_infnan_mod, only: is_nan => shr_infnan_isnan

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
    character(len=512), intent(out) :: errmsg
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

    return

  end subroutine stats_end_timestep_clubb
#endif

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

#ifdef CLUBB_SGS
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

    return

  end subroutine stats_avg
#endif

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------

#ifdef CLUBB_SGS
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

  real(kind_phys) function diag_ustar( z, bflx, wnd, z0, shr_const_karman, &
                                       shr_const_pi, shr_const_g )

!    use shr_const_mod, only : shr_const_karman, shr_const_pi, shr_const_g

    implicit none

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

    return

  end function diag_ustar
#endif


end module clubb

