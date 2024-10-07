module vertical_diffusion_solver

    implicit none
    private
    save
    
    public :: fin_vol_solve
    
    contains
    
    ! Designed to solve the equation:
    !
    ! w * dq/dt = d/dp (D q' - v q) + c q
    !
    ! where q is a grid-cell average, and p is the vertical coordinate
    ! (presumably pressure).
    !
    ! In this function, coef_q_weight == w, coef_q_diff == D,
    ! coef_q_adv == v, and coef_q == c. All these are optional; omitting a
    ! coefficient is equivalent to setting the entire array to 0.
    !
    ! coef_q_diff and coef_q_adv are defined at the level interfaces, while
    ! coef_q and coef_q_weight are grid-cell averages.

    function fin_vol_solve(dt, p, toSolve, ncols, pver, coef_q, coef_q_diff, coef_q_adv, &
        coef_q_weight, upper_bndry, lower_bndry, l_cond, r_cond)  result(solution)
    
        use linear_1d_operators, only: &
        zero_operator,               &
        diagonal_operator,           &
        diffusion_operator,          &
        advection_operator,          &
        BoundaryType,                &
        TriDiagDecomp,               &
        TriDiagOp,                   &
        BoundaryCond,                &
        operator(+)
        use shr_kind_mod,        only: r8 => shr_kind_r8
        use coords_1d,           only: Coords1D
    
        ! ---------------------- !
        ! Input-Output Arguments !
        ! ---------------------- !
    
        ! Time step.
        real(r8), intent(in) :: dt
        ! Grid spacings.
        type(Coords1D), intent(in) :: p
    
        ! Matrix to decomp from.
        ! real(r8), intent(in) :: u(ncols,pver)
        integer,  intent(in)    :: ncols
        integer,  intent(in)    :: pver
        real(r8), intent(in)    :: toSolve(ncols,pver)
    
        ! Coefficients for diffusion and advection.
        !
        ! The sizes must be consistent among all the coefficients that are
        ! actually present, i.e. coef_q_diff and coef_q_adv should be one level
        ! bigger than coef_q and coef_q_weight, and have the same column number.
        real(r8), contiguous, intent(in), optional :: coef_q(:,:), &
        coef_q_diff(:,:), coef_q_adv(:,:), coef_q_weight(:,:)
    
        ! Boundary conditions (optional, default to 0 flux through boundary).
        class(BoundaryType), target, intent(in), optional :: &
        upper_bndry, lower_bndry
    
        ! Objects representing boundary conditions.
        class(BoundaryCond), intent(in), optional :: l_cond, r_cond
    
        real(r8) :: solution(ncols,pver)
    
        ! decomposition.
        type(TriDiagDecomp) :: decomp
    
        ! --------------- !
        ! Local Variables !
        ! --------------- !
    
        ! Operator objects.
        type(TriDiagOp) :: add_term
        type(TriDiagOp) :: net_operator
    
        ! ----------------------- !
        ! Main Computation Begins !
        ! ----------------------- !
    
        ! A diffusion term is probably present, so start with that. Otherwise
        ! start with an operator of all 0s.
    
        if (present(coef_q_diff)) then
        net_operator = diffusion_operator(p, coef_q_diff, &
            upper_bndry, lower_bndry)
        else
        net_operator = zero_operator(p%n, p%d)
        end if
    
        ! Constant term (damping).
        if (present(coef_q)) then
        add_term = diagonal_operator(coef_q)
        call net_operator%add(add_term)
        end if
    
        ! Effective advection.
        if (present(coef_q_adv)) then
        add_term = advection_operator(p, coef_q_adv, &
            upper_bndry, lower_bndry)
        call net_operator%add(add_term)
        end if
    
        ! We want I-dt*(w^-1)*A for a single time step, implicit method, where
        ! A is the right-hand-side operator (i.e. what net_operator is now).
        if (present(coef_q_weight)) then
        call net_operator%lmult_as_diag(-dt/coef_q_weight)
        else
        call net_operator%lmult_as_diag(-dt)
        end if
        call net_operator%add_to_diag(1._r8)
    
        ! Decompose
        decomp = TriDiagDecomp(net_operator)
        solution = toSolve
    
        call net_operator%finalize()
        call add_term%finalize()
    
        call decomp%left_div(solution(:ncols, :), l_cond=l_cond, r_cond=r_cond)
        !tendency = tendency - toSolve
    
        ! Ensure local objects are deallocated.
        call decomp%finalize()
    
    end function fin_vol_solve
    
end module vertical_diffusion_solver
    