#include "Sloan.hpp"

void Sloan::solve(
    Cauchy &sigma_prime_ep, 
    State &state_ep, 
    Voigt Delta_epsilon_tilde_p, 
    double FTOL,
    double STOL,
    double DT_MIN,
    double EPS,
    int MAXITS_YSC,
    YieldFunction compute_f,
    ConstitutiveMatrixFunction compute_D_e,
    DerivativeFunction compute_derivatives,
    StateIncrementFunction compute_plastic_state_variable_increment) {
    // Sloan et al. (2001) substepping with automatic error control.
    int substeps = 0;
    int corrections = 0;
    double dT = 1.0;
    double T = 0.0;
    while (T < 1.0) {
        PLOG_DEBUG << "Plastic increment, dT = " << dT;
        Voigt Delta_epsilon_tilde_p_dT = Delta_epsilon_tilde_p*dT;
        double Delta_epsilon_vol_p_dT = to_cauchy(Delta_epsilon_tilde_p_dT).trace();

        // Compute stress and state variable increment estimates.
        Voigt Delta_sigma_prime_1, Delta_sigma_prime_2;
        State delta_state_1, delta_state_2;
        Sloan::compute_plastic_increment(sigma_prime_ep, state_ep, Delta_epsilon_tilde_p_dT, Delta_sigma_prime_1, delta_state_1, compute_D_e, compute_derivatives, compute_plastic_state_variable_increment);
        Cauchy sigma_prime_1 = sigma_prime_ep + to_cauchy(Delta_sigma_prime_1);
        State state_1 = state_ep + delta_state_1;
        Sloan::compute_plastic_increment(sigma_prime_1, state_ep, Delta_epsilon_tilde_p_dT, Delta_sigma_prime_2, delta_state_2, compute_D_e, compute_derivatives, compute_plastic_state_variable_increment);
        PLOG_DEBUG << "State variable increment 1, delta_state_1 = " << delta_state_1;
        PLOG_DEBUG << "State variable increment 2, delta_state_2 = " << delta_state_2;

        // Calculate modified Euler stresses and state variables.
        Cauchy sigma_prime_ini = sigma_prime_ep + to_cauchy(1.0/2.0*(Delta_sigma_prime_1 + Delta_sigma_prime_2));
        State state_ini = state_ep + 1.0/2.0*(delta_state_1+delta_state_2);
        PLOG_DEBUG << "Initial stress estimate, sigma_prime_ini_tilde = " << to_voigt(sigma_prime_ini);
        PLOG_DEBUG << "State variable estimate, state_ini = " << state_ini;

        // Compute error estimate.
        bool accepted = false;
        double R_n = Sloan::compute_error_estimate(sigma_prime_ini, state_ini, Delta_sigma_prime_1, Delta_sigma_prime_2, delta_state_1, delta_state_2, EPS);
        if (R_n < STOL) {
            // Accept increment and correct stresses and state variables back to yield surface.
            accepted = true;
            Cauchy sigma_prime_u, sigma_prime_c;
            sigma_prime_u = sigma_prime_c = sigma_prime_ini;
            State state_u, state_c;
            state_u = state_c = state_ini;

            // Correct stresses and state variables back to yield surface. 
            double f_u, f_c;
            f_u = f_c = compute_f(sigma_prime_u, state_u);
            int ITS_YSC = 0;
            while (ITS_YSC < MAXITS_YSC && std::abs(f_c) > FTOL) {
                // If yield surface drift correction unsuccessful, log fault.
                if (ITS_YSC >= MAXITS_YSC && std::abs(f_c) > FTOL) {
                    PLOG_FATAL << "Maximum number of yield surface correction iterations performed and f = " << f_c << " > FTOL = " << FTOL << ".";
                    assert(false);
                } 

                // Calculate uncorrected elastic constitutive matrix using tangent moduli and elastic stress increment.
                Constitutive D_e_u = compute_D_e(sigma_prime_u, Cauchy::Zero());

                // Calculate uncorrected derivatives.
                Cauchy df_dsigma_prime_u, dg_dsigma_prime_u;
                Voigt a_u, b_u;
                double H_u;
                compute_derivatives(sigma_prime_u, state_u, df_dsigma_prime_u, a_u, dg_dsigma_prime_u, b_u, H_u);

                // Compute consistent correction to yield surface.
                Sloan::compute_yield_surface_correction(sigma_prime_u, state_u, f_u, H_u, a_u, b_u, D_e_u, sigma_prime_c, state_c, compute_plastic_state_variable_increment);

                // Check yield surface function.
                f_c = compute_f(sigma_prime_c, state_c);
                f_u = compute_f(sigma_prime_u, state_u);
                if (std::abs(f_c) > std::abs(f_u)) {
                    // Apply normal correction instead.
                    Sloan::compute_normal_yield_surface_correction(sigma_prime_u, state_u, f_u, a_u, sigma_prime_c, state_c);
                }

                // Correct stress and state variables.
                sigma_prime_u = sigma_prime_c;
                state_u = state_c;
                ITS_YSC += 1;
            }

            // Update stress state and state variables.
            corrections += ITS_YSC;
            sigma_prime_ep = sigma_prime_c;
            state_ep = state_c;

            // Increment substep and pseudotime T.
            substeps += 1;
            T += dT;
        }
        dT = Sloan::compute_new_substep_size(accepted, dT, T, R_n, DT_MIN, STOL);
    }
    // Strain increment solved.
    PLOG_INFO << "Elastoplastic stress increment integrated to within a tolerance FTOL = " << FTOL << " via " << substeps << " substeps and " << corrections << " drift corrections.";
}

double Sloan::compute_new_substep_size(
    bool accepted, 
    double dT, 
    double T, 
    double R_n, 
    double DT_MIN, 
    double STOL) {
    if (dT == DT_MIN) {
        PLOG_FATAL << "Minimum step size DT_MIN = " << DT_MIN << " failed to generated an accepted increment.";
        assert(false);
    }
    double q_step;
    if (accepted) {
        // Step size accepted: allow substep size to grow.
        q_step = std::min(0.9*std::sqrt(STOL/R_n),1.1);
    } else {
        // Step size rejected: calculate reduced substep size factor and limit substep size growth factor.
        q_step = std::max(0.9*std::sqrt(STOL/R_n),0.1);
        q_step = std::min(q_step, 1.0);
    }
    dT *= q_step;

    // Delimit substep size to pseudomtime remaining and minimum substep size.
    dT = std::max(dT, DT_MIN);
    dT = std::min(dT, 1.0-T);
    return dT;
}

void Sloan::compute_yield_surface_correction(
    Cauchy sigma_prime_u, 
    State state_u, 
    double f_u, 
    double H_u, 
    Voigt a_u, 
    Voigt b_u, 
    Constitutive D_e_u, 
    Cauchy &sigma_prime_c, 
    State &state_c, 
    StateIncrementFunction compute_plastic_state_variable_increment) {
    // Compute correction factor.
    double delta_lambda_c = f_u/(H_u + a_u.transpose()*D_e_u*b_u);

    // Update stress and state variables using corrections.
    Voigt Delta_sigma_prime_c = -delta_lambda_c*D_e_u*b_u;
    Cauchy df_dsigma_prime_u = to_cauchy(a_u);
    State delta_state_c = compute_plastic_state_variable_increment(delta_lambda_c, df_dsigma_prime_u, H_u, Voigt::Zero());
    sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
    state_c = state_u + delta_state_c;
}

void Sloan::compute_normal_yield_surface_correction(
    Cauchy sigma_prime_u, 
    State state_u, 
    double f_u, 
    Voigt a_u, 
    Cauchy &sigma_prime_c, 
    State &state_c) {
    // Compute corretion factor.
    double delta_lambda_c = f_u/(a_u.transpose()*a_u);
    
    //Update stress and state variables using correction.
    Voigt Delta_sigma_prime_c = -delta_lambda_c*a_u;
    sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
    state_c = state_u; /* i.e. no correction to state variables. */
}

void Sloan::compute_plastic_increment(
    Cauchy sigma_prime, 
    State state, 
    Voigt Delta_epsilon_tilde_p_dT, 
    Voigt &Delta_sigma_prime, 
    State &delta_state, 
    ConstitutiveMatrixFunction compute_D_e, 
    DerivativeFunction compute_derivatives, 
    StateIncrementFunction compute_plastic_state_variable_increment) {   
    // Calculate elastic constitutive matrix using tangent moduli and elastic stress increment.
    Constitutive D_e = compute_D_e(sigma_prime, Cauchy::Zero());

    // Compute elastoplastic constitutive matrix and elastoplastic multiplier. 
    Cauchy df_dsigma_prime, dg_dsigma_prime;
    Voigt a, b;
    double H;
    compute_derivatives(sigma_prime, state, df_dsigma_prime, a, dg_dsigma_prime, b, H);
    Constitutive D_ep = D_e-(D_e*b*a.transpose()*D_e)/(H+a.transpose()*D_e*b);
    Voigt Delta_sigma_prime_e = D_e*Delta_epsilon_tilde_p_dT;
    double delta_lambda = (double)(a.transpose()*Delta_sigma_prime_e)/(double)(H + a.transpose()*D_e*b);
    PLOG_DEBUG << "Plastic multiplier, delta_lambda = " << delta_lambda;

    // Update stress and state variable increment by reference.
    Delta_sigma_prime = D_ep*Delta_epsilon_tilde_p_dT;
    delta_state = compute_plastic_state_variable_increment(delta_lambda, df_dsigma_prime, H, Delta_epsilon_tilde_p_dT);
}

double Sloan::compute_error_estimate(
    Cauchy sigma_prime_ini, 
    State state_ini, 
    Voigt Delta_sigma_prime_1, 
    Voigt Delta_sigma_prime_2, 
    State delta_state_1, 
    State delta_state_2, 
    double EPS) {
    int size_state = state_ini.size();
    State error(size_state);
    error[0] = (to_cauchy(Delta_sigma_prime_2 - Delta_sigma_prime_1)).norm()/sigma_prime_ini.norm();
    for (int i=1; i<size_state; i++) {
        error[i] = std::abs((delta_state_2[i] - delta_state_1[i]))/state_ini[i];
    }
    double R_n = 1.0/2.0*error.maxCoeff();

    // Check against machine tolerance.
    R_n = std::max(R_n, EPS);
    return R_n;
}