#include "Elastoplastic.hpp"

void Elastoplastic::solve(void) {
    if (!solved) {
        
        PLOG_INFO << "Solving current strain increment.";
        State state = get_state_variables();
        alpha = compute_alpha(sigma_prime, state, Delta_epsilon_tilde, FTOL);
        substeps = 0;
        corrections = 0;

        // Update stress and state variables based on elastic portion of strain increment.
        if (alpha > 0.0) {
            Delta_epsilon_tilde_e = alpha*Delta_epsilon_tilde;
            Delta_epsilon_e = to_cauchy(Delta_epsilon_tilde_e);
            Delta_epsilon_vol_e = compute_Delta_epsilon_vol(Delta_epsilon_e);
            Delta_sigma_prime_e = compute_elastic_stress_increment(D_e, Delta_epsilon_tilde_e);
            sigma_prime_e = compute_isotropic_linear_elastic_stress(sigma_prime, alpha, Delta_epsilon_tilde_e);
            state_e = compute_elastic_state_variable(Delta_epsilon_tilde_e);
            PLOG_INFO << "Elastic portion of strain increment alpha = " << alpha << " applied.";
        }
        // If fully plastic increment, use current stress and state variables.
        if (alpha == 0.0) {
            sigma_prime_e = sigma_prime;
            state_e = get_state_variables();
        }
        PLOG_DEBUG << "sigma_prime_e = \n" << sigma_prime_e;
        PLOG_DEBUG << "state_e = \n" << state_e;

        if (alpha < 1.0) {
            // Perform elastoplastic stress integration on remaining plastic portion of strain increment.
            Delta_epsilon_tilde_p = (1.0-alpha)*Delta_epsilon_tilde;
            sigma_prime_ep = sigma_prime_e;
            state_ep = state_e;

            // Substepping with automatic error control.
            dT = 1.0;
            T = 0.0;
            while (T < 1.0) {
                Delta_epsilon_tilde_p_dT = Delta_epsilon_tilde_p*dT;
                Delta_epsilon_vol_p_dT = to_cauchy(Delta_epsilon_tilde_p_dT).trace();

                // Compute stress and state variable increment estimates.
                compute_plastic_increment(sigma_prime_ep, state_ep, Delta_epsilon_tilde_p_dT, Delta_sigma_prime_1, delta_state_1);
                sigma_prime_1 = sigma_prime_ep + to_cauchy(Delta_sigma_prime_1);
                state_1 = state_ep + delta_state_1;
                compute_plastic_increment(sigma_prime_1, state_ep, Delta_epsilon_tilde_p_dT, Delta_sigma_prime_2, delta_state_2);

                // Debug output.
                PLOG_DEBUG << "sigma_prime_ep = \n" << sigma_prime_ep;
                PLOG_DEBUG << "state_ep = \n" << state_ep;
                PLOG_DEBUG << "sigma_prime_1 = \n" << sigma_prime_1;
                PLOG_DEBUG << "state_1 = \n" << state_1;
                PLOG_DEBUG << "Delta_sigma_prime_1 = " << Delta_sigma_prime_1;
                PLOG_DEBUG << "delta_state_1 = " << delta_state_1;
                PLOG_DEBUG << "Delta_sigma_prime_2 = " << Delta_sigma_prime_2;
                PLOG_DEBUG << "delta_state_2 = " << delta_state_2;

                // Calculate modified Euler stresses and state variables.
                sigma_prime_ini = sigma_prime_ep + to_cauchy(1.0/2.0*(Delta_sigma_prime_1 + Delta_sigma_prime_2));
                state_ini = state_ep + 1.0/2.0*(delta_state_1+delta_state_2);

                // Compute error estimate.
                accepted = false;
                R_n = compute_error_estimate();
                if (R_n < STOL) {
                    // Accept increment and correct stresses and state variables back to yield surface.a
                    accepted = true;
                    sigma_prime_u = sigma_prime_c = sigma_prime_ini;
                    state_u = state_c = state_ini;

                    // Debug output.
                    PLOG_DEBUG << "sigma_prime_u =\n" << sigma_prime_u;
                    PLOG_DEBUG << "state_u = " << state_u;

                    // Correct stresses and state variables back to yield surface. 
                    f_0 = f_u = f_c = compute_f(sigma_prime_u, state_u);

                    ITS_YSC = 0;
                    while (ITS_YSC < MAXITS_YSC) {
                        if (std::abs(f_c) > FTOL) {
                            // Compute consistent correction to yield surface.
                            compute_yield_surface_correction(f_0, f_u, sigma_prime_u, state_u, f_c, sigma_prime_c, state_c);

                            // Debug output.
                            PLOG_DEBUG << "sigma_prime_c =\n" << sigma_prime_c;
                            PLOG_DEBUG << "state_c =\n" << state_c;
                            PLOG_DEBUG << "f_c = " << f_c;

                            // Correct stress and state variables.
                            sigma_prime_u = sigma_prime_c;
                            state_u = state_c;
                            ITS_YSC += 1;
                        } else {
                            // Stresses and state variables corrected; breakout of correction while loop.
                            sigma_prime_c = sigma_prime_u;
                            state_c = state_u;
                            goto CORRECTED;
                        }
                        // If yield surface drift correction unsuccessful, log fault.
                        if (ITS_YSC >= MAXITS_YSC && std::abs(f_c) > FTOL) {
                            PLOG_FATAL << "Maximum number of yield surface correction iterations performed and f = " << f_c << " > FTOL = " << FTOL << ".";
                            assert(false);
                        } 
                    }

                    // Update stress state and state variables.
                    CORRECTED: corrections += ITS_YSC;
                    sigma_prime_ep = sigma_prime_c;
                    state_ep = state_c;

                    // Increment substep and pseudotime T.
                    substeps += 1;
                    T += dT;
                }
                dT = compute_new_substep_size();
            }
            // Set stress and state variables as final elastoplastic values.
            sigma_prime = sigma_prime_ep;
            set_state_variables(state_ep);
            solved = true;
            PLOG_INFO << "Elastoplastic strain increment integrated to within a tolerance FTOL = " << FTOL << " via " << substeps << " substeps and " << corrections << " drift corrections.";
        } else {
            // Fully elastic increment. Update stress and state variables.
            sigma_prime = sigma_prime_e;
            set_state_variables(state_e);
            PLOG_INFO << "Fully elastic stress increment integrated directly.";
        }

        // Compute final stress variables.
        compute_stress_variables();
    }
}

double Elastoplastic::compute_new_substep_size(void) {
    if (dT == dT_min) {
        PLOG_FATAL << "Minimum step size dT_min = " << dT_min << " failed to generated an accepted increment.";
        assert(false);
    }
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
    dT = std::max(dT, dT_min);
    dT = std::min(dT, 1.0-T);

    // Debug output.
    PLOG_DEBUG << "dT = " << dT;

    return dT;
}

void Elastoplastic::compute_yield_surface_correction(double f_0, double f_u, Cauchy sigma_prime_u, State state_u, double &f_c, Cauchy &sigma_prime_c, State &state_c) {
    // Calculate elastic constitutive matrix using tangent moduli and elastic stress increment.
    double p_prime_u = compute_p_prime(sigma_prime_u);
    double K_tan_u = compute_K(0, p_prime_u);
    double G_tan_u = compute_G(K_tan_u);
    Constitutive D_e_u = compute_isotropic_linear_elastic_matrix(K_tan_u, G_tan_u);

    // Calculate derivatives.
    Cauchy df_dsigma_prime_u, dg_dsigma_prime_u;
    Voigt a_u, b_u;
    double H_u;
    compute_derivatives(sigma_prime_u, state_u, df_dsigma_prime_u, a_u, dg_dsigma_prime_u, b_u, H_u);

    // Compute correction factor.
    double delta_lambda_c = f_u/(H_u + a_u.transpose()*D_e_u*b_u);

    // Update stress and state variables using corrections.
    Voigt Delta_sigma_prime_c = -delta_lambda_c*D_e_u*b_u;
    State delta_state_c = compute_plastic_state_variable_increment(sigma_prime_c, delta_lambda_c, H_u);
    sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
    state_c = state_u + delta_state_c;

    // Check yield surface function.
    f_c = compute_f(sigma_prime_c, state_c);

    // Apply normal correction instead.
    if (std::abs(f_c) > std::abs(f_0)) {
        // Compute corretion factor.
        double delta_lambda_c = f_u/(a_u.transpose()*a_u);

        //Update stress and state variables using correction.
        Voigt Delta_sigma_prime_c = -delta_lambda_c*a_u;
        sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
        state_c = state_u; /* i.e. no correction to state variables. */

        // Recalculate yield surface function.
        f_c = compute_f(sigma_prime_c, state_c);
    }
    std::cout << "delta_lambda = " << delta_lambda_c <<std::endl;
    std::cout << "Delta_sigma_prime_c = \n" << Delta_sigma_prime_c << std::endl;
    std::cout << "a_u = " << a_u << std::endl;
    std::cout << "b_u = " << b_u << std::endl;
    std::cout << "D_e_u = \n" << D_e_u << std::endl;
    std::cout << "H_u = " << H_u << std::endl;
    std::cout << "f_c = " << f_c << std::endl;
    std::cout << "sigma_prime_c = \n" << sigma_prime_c << std::endl;
}

void Elastoplastic::compute_plastic_increment(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde_p_dT, Voigt &Delta_sigma_prime, State &delta_state) {   
    // Calculate elastic constitutive matrix using tangent moduli and elastic stress increment.
    double p_prime = compute_p_prime(sigma_prime);
    double K_tan = compute_K(0, p_prime);
    double G_tan = compute_G(K_tan);
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K_tan, G_tan);
    Voigt Delta_sigma_prime_e = D_e*Delta_epsilon_tilde_p_dT;

    // Compute elastoplastic constitutive matrix and elastoplastic multiplier. 
    Cauchy df_dsigma_prime, dg_dsigma_prime;
    Voigt a, b;
    double H;
    compute_derivatives(sigma_prime, state, df_dsigma_prime, a, dg_dsigma_prime, b, H);
    Constitutive D_ep = compute_elastoplastic_matrix(D_e, a, b, H);
    double delta_lambda = compute_elastoplastic_multiplier(Delta_sigma_prime_e, D_e, a, b, H);

    // Update stress and state variable increment by reference.
    Delta_sigma_prime = D_ep*Delta_epsilon_tilde_p_dT;
    delta_state = compute_plastic_state_variable_increment(sigma_prime, Delta_epsilon_tilde_p_dT, delta_lambda, H);
}

double Elastoplastic::compute_error_estimate(void) {
    int size_state = get_state_variables().size();
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

double Elastoplastic::compute_elastoplastic_multiplier(Voigt Delta_sigma_prime_e, Constitutive D_e, Voigt a, Voigt b, double H) {
    double lambda = (double)(a.transpose()*Delta_sigma_prime_e)/(double)(H + a.transpose()*D_e*b);

    // Debug output.
    PLOG_DEBUG << "lambda = " << lambda;

    return lambda;
}

Constitutive Elastoplastic::compute_elastoplastic_matrix(Constitutive D_e, Voigt a, Voigt b, double H) {
    Constitutive D_ep = D_e-(D_e*b*a.transpose()*D_e)/(H+a.transpose()*D_e*b);

    // Debug output.
    PLOG_DEBUG << "D_ep = \n" << D_ep;

    return D_ep;
}

// void Elastoplastic::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H) {
//     // Current state variables.
//     double e = state[0];
//     double p_c = state[1];

//     // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
//     double q = compute_q(sigma_prime);
//     double p_prime = compute_p_prime(sigma_prime);
//     Cauchy sigma = compute_sigma(sigma_prime, u);
//     Cauchy s = compute_s(sigma, p);
//     double p = compute_p(sigma);
//     double I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar;
//     compute_stress_invariants(sigma, p, s, I_1, I_2, I_3, J_1, J_2, J_3);
//     compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
//     Cauchy dq_dsigma_prime = compute_dq_dsigma_prime(sigma_prime, s, q);
//     Cauchy dtheta_dsigma_prime = compute_dtheta_dsigma_prime(sigma_prime);
    
//     // Compute derivatives of the yield surface.
//     double df_dq = compute_df_dq();
//     double df_dp_prime = compute_df_dp_prime();
//     double df_dtheta = compute_df_dtheta();
//     df_dsigma_prime = df_dq*dq_dsigma_prime + df_dp_prime*dp_prime_dsigma_prime + df_dtheta*dtheta_dsigma_prime;

//     // Compute derivatives of the plastic potential function.
//     double dg_dq = compute_dg_dq();
//     double dg_dp_prime = compute_dg_dp_prime();
//     double dg_dtheta = compute_dg_dtheta();
//     dg_dsigma_prime = dg_dq*dq_dsigma_prime + dg_dp_prime*dp_prime_dsigma_prime + dg_dtheta*dtheta_dsigma_prime;
    
//     // Convert to Voigt notation.
//     a = to_voigt(df_dsigma_prime);
//     b = to_voigt(dg_dsigma_prime);

//     // Compute hardening modulus.
//     H = compute_H();
// }

double Elastoplastic::compute_alpha(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde, double FTOL){
    PLOG_DEBUG << "Computing alpha.";
    double alpha;

    // Current stress state (i.e. alpha = 0.0).
    double f_0 = compute_f_alpha(sigma_prime, state, 0.0, Delta_epsilon_tilde);

    // Fully elastic increment trial stress state.
    double f_1 = compute_f_alpha(sigma_prime, state, 1.0, Delta_epsilon_tilde);

    // Check increment type by finding alpha.
    if (f_1 <= FTOL) {
        PLOG_INFO << "Fully elastic increment; alpha = 1.0.";
        alpha = 1.0;
    } else if (std::abs(f_0) <= FTOL && f_1> FTOL) {
        PLOG_INFO << "Checking potential for elastoplastic unloading-reloading.";
        if (!check_unload_reload(sigma_prime, state)) {
            // Compute bounds on alpha.
            double alpha_0 = 0.0;
            double alpha_1 = 1.0;
            compute_alpha_bounds(sigma_prime, state, alpha_0, alpha_1);

            // Trial stress state for alpha_0.
            double f_0 = compute_f_alpha(sigma_prime, state, alpha_0, Delta_epsilon_tilde);

            // Trial stress state for alpha_1.
            double f_1 = compute_f_alpha(sigma_prime, state, alpha_1, Delta_epsilon_tilde);

            // Perform Pegasus intersection within the bounds of alpha_0 and alpha_1.
            PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = " << alpha_0 << " and alpha_1 = " << alpha_1 <<".";
            alpha = pegasus_regula_falsi(sigma_prime, state, alpha_0, alpha_1, f_0, f_1);
        } else {
            PLOG_INFO << "Fully plastic increment; alpha = 0.0.";
            alpha = 0.0;
        } 
    } else if (f_0 < -FTOL && f_1 > FTOL) {
        PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = 0.0 and alpha_1 = 1.0.";
        alpha = pegasus_regula_falsi(sigma_prime, state, 0.0, 1.0, f_0, f_1);
    } else {
        PLOG_FATAL << "Illegal stress state.";
        assert(false);
    }
    return alpha;
}

bool Elastoplastic::check_unload_reload(Cauchy sigma_prime, State state) {
    // Compute required derivatives for the given stress state.
    Cauchy df_dsigma_prime_check, dg_dsigma_prime_check;
    Voigt a_check, b_check;
    double H_check;
    compute_derivatives(sigma_prime, state, df_dsigma_prime_check, a_check, dg_dsigma_prime_check, b_check, H_check);

    // Compute the elastic matrix using tangent moduli.
    double K_tan = compute_K(0, p_prime);
    double G_tan = compute_G(K_tan);
    Constitutive D_e_tan = compute_isotropic_linear_elastic_matrix(K_tan, G_tan);

    // Compute elastic stress increment.
    Voigt Delta_sigma_e = D_e_tan*Delta_epsilon_tilde;

    // Check unloading-reloading criterion.
    double cos_theta = (double)(a.transpose()*Delta_sigma_e)/(double)(a.squaredNorm()*Delta_sigma_e.squaredNorm());

    // Check against tolerance.
    if (cos_theta >= -LTOL) {
        return false;
    } else {
        return true;
    }
}

void Elastoplastic::compute_alpha_bounds(Cauchy sigma_prime, State state, double &alpha_0, double &alpha_1) {
    PLOG_DEBUG << "Computing alpha bounds.";
    double Delta_epsilon_vol_e = Delta_epsilon.trace();
    double alpha_n, d_alpha;
    int i = 0;
    int j = 0;
    while (i < MAXITS_YSI) {
        d_alpha = (alpha_1-alpha_0)/NSUB;
        alpha_n = alpha_0+d_alpha;
        while (j < NSUB) {
            // Compute the elastic matrix.
            double K_trial = compute_K(alpha_n*Delta_epsilon_vol_e, p_prime);
            double G_trial = compute_G(K_trial);
            Constitutive D_e_trial = compute_isotropic_linear_elastic_matrix(K_trial, G_trial);

            // Compute elastic stress increment.
            Voigt Delta_sigma_e_trial = alpha_n*D_e_trial*Delta_epsilon_tilde;

            // Check yield function.
            State state_trial = state;
            Cauchy sigma_prime_trial = sigma_prime + to_cauchy(Delta_sigma_e_trial);
            double f_trial = compute_f(sigma_prime_trial, state_trial);

            // Check criterion.
            if (f_trial > FTOL) {
                // Break from loops.
                alpha_1 = alpha_n;
                goto BREAK;
            } else {
                // Continue iterating.
                alpha_0 = alpha_n;
                alpha_n = alpha_n+d_alpha;
                j += 1; 
            }
        }
        i += 1;
    }
    BREAK: return;
}

double Elastoplastic::pegasus_regula_falsi(Cauchy sigma_prime, State state, double alpha_0, double alpha_1, double f_0, double f_1) {
    // Pegasus algorithm.
    double alpha, alpha_n;
    double f_n = FTOL;
    Cauchy sigma_prime_n = sigma_prime;
    State state_n = state;

    // Iterate to find optimal alpha if a plastic increment.
    ITS_YSI = 0;
    while (ITS_YSI < MAXITS_YSI && std::abs(f_n) >= FTOL) {
        alpha_n = alpha_1 - f_1*(alpha_1-alpha_0)/(f_1-f_0);  
            
        sigma_prime_n = compute_isotropic_linear_elastic_stress(sigma_prime, alpha_n, Delta_epsilon_tilde);
        f_n = compute_f(sigma_prime_n, state_n);

        // Update trial using Pegasus method rules.
        if (f_n*f_1 < 0) {
            alpha_1 = alpha_0;
            f_1 = f_0;
        } else {
            f_1 = f_1*f_0/(f_0+f_n);
        }
        alpha_0 = alpha_n;
        f_0 = f_n;
        ITS_YSI += 1;
    }
    alpha = alpha_n;
    alpha = std::max(alpha, 0.0);
    alpha = std::min(alpha, 1.0);
    double f = f_n;
    if (std::abs(f) >= FTOL) {
        PLOG_FATAL << "Performed " << MAXITS_YSI << " Pegasus iteration(s): alpha = " << alpha << "; |f| = " << std::abs(f) << " > tolerance = " << FTOL << ".";
        assert(false);
    } else{
        PLOG_INFO << "Performed " << ITS_YSI << " Pegasus iteration(s): " << "alpha = " << alpha << "; " << "|f| = " << std::abs(f) << " < tolerance = " << FTOL << ".";
        assert(alpha >= 0.0 && alpha <= 1.0);
    }
    return alpha;
}

double Elastoplastic::compute_f_alpha(Cauchy sigma_prime, State state, double alpha, Voigt Delta_epsilon_tilde) {
    Cauchy sigma_prime_alpha = compute_isotropic_linear_elastic_stress(sigma_prime, alpha, Delta_epsilon_tilde);
    return compute_f(sigma_prime_alpha, state);
}