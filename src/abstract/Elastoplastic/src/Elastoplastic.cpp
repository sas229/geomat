#include "Elastoplastic.hpp"

void Elastoplastic::solve(void) {
    if (!solved) {
        compute_alpha();
        substeps = 0;
        corrections = 0;

        if (alpha == 0.0) {
            // Fully plastic increment. Update stress and state variables.
            sigma_prime_e = sigma_prime;
            state_e = get_state_variables();
        }
        if (alpha > 0.0) {
            // Update stress and state variables based on elastic portion of strain increment.
            Delta_epsilon_tilde_e = alpha*Delta_epsilon_tilde;
            Delta_epsilon_e = to_cauchy(Delta_epsilon_tilde_e);
            Delta_epsilon_vol_e = compute_Delta_epsilon_vol(Delta_epsilon_e);
            Delta_sigma_prime_e = compute_elastic_stress_increment(D_e, Delta_epsilon_tilde_e);
            sigma_prime_e = compute_isotropic_linear_elastic_stress(sigma_prime, alpha, Delta_epsilon_tilde_e);
            state_e = compute_elastic_state_variable(Delta_epsilon_tilde_e);
        } 
        if (alpha < 1.0) {
            // Perform elastoplastic stress integration on plastic portion of strain increment.
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

                // Calculate modified Euler stresses and state variables.
                sigma_prime_ini = sigma_prime_ep + to_cauchy(1.0/2.0*(Delta_sigma_prime_1 + Delta_sigma_prime_2));
                state_ini = state_ep + 1.0/2.0*(delta_state_1+delta_state_2);

                // Compute error estimate.
                accepted = false;
                compute_error_estimate();
                if (R_n < STOL) {
                    // Accept increment and correct stresses and state variables back to yield surface.a
                    accepted = true;
                    sigma_prime_u = sigma_prime_c = sigma_prime_ini;
                    state_u = state_c = state_ini;

                    // Correct stresses and state variables back to yield surface. 
                    f_u = f_c = compute_f(sigma_prime_u, state_u);
                    ITS_YSC = 0;
                    while (ITS_YSC < MAXITS_YSC) {
                        if (std::abs(f_c) > FTOL) {
                            // Compute consistent correction to yield surface.
                            compute_yield_surface_correction();

                            // Check yield surface function.
                            f_c = compute_f(sigma_prime_c, state_c);
                            f_u = compute_f(sigma_prime_u, state_u);
                            if (std::abs(f_c) > std::abs(f_u)) {
                                // Apply normal correction instead.
                                compute_normal_yield_surface_correction();
                            }

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
            PLOG_INFO << "Elastoplastic stress increment integrated to within a tolerance FTOL = " << FTOL << " via " << substeps << " substeps and " << corrections << " drift corrections.";
        } else {
            // Fully elastic increment. Update stress and state variables.
            sigma_prime = sigma_prime_e;
            set_state_variables(state_e);
            PLOG_INFO << "Fully elastic stress increment integrated directly.";
        }

        // Compute final stress invariants.
        p_prime = compute_p_prime(sigma_prime);
        q = compute_q(sigma_prime);
        compute_principal_stresses(sigma_prime, sigma_1, sigma_2, sigma_3, R, S);
        Eigen::VectorXd p_prime_surface, q_surface;
        compute_yield_surface(get_state_variables(), 1000, p_prime_surface, q_surface);
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
    return dT;
}

void Elastoplastic::compute_yield_surface_correction(void) {
    // Calculate elastic constitutive matrix using tangent moduli and elastic stress increment.
    p_prime_u = compute_p_prime(sigma_prime_u);
    K_tan_u = compute_K(0, p_prime_u);
    G_tan_u = compute_G(K_tan_u);
    D_e_u = compute_isotropic_linear_elastic_matrix(K_tan_u, G_tan_u);

    // Calculate derivatives.
    compute_derivatives(sigma_prime_u, state_u, df_dsigma_prime_u, a_u, dg_dsigma_prime_u, b_u, dg_dp_prime_u, H_u);

    // Compute correction factor.
    delta_lambda_c = f_u/(H_u + a_u.transpose()*D_e_u*b_u);

    // Update stress and state variables using corrections.
    Delta_sigma_prime_c = -delta_lambda_c*D_e_u*b_u;
    delta_state_c = compute_plastic_state_variable_increment(delta_lambda_c, H_u);
    sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
    state_c = state_u + delta_state_c;
}

void Elastoplastic::compute_normal_yield_surface_correction(void) {
    // Compute corretion factor.
    delta_lambda_c = f_u/(a_u.transpose()*a_u);
    
    //Update stress and state variables using correction.
    Delta_sigma_prime_c = -delta_lambda_c*a_u;
    sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
    state_c = state_u; /* i.e. no correction to state variables. */
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
    double dg_dp_prime, H;
    compute_derivatives(sigma_prime, state, df_dsigma_prime, a, dg_dsigma_prime, b, dg_dp_prime, H);
    Constitutive D_ep = compute_elastoplastic_matrix(D_e, a, b, H);
    double delta_lambda = compute_elastoplastic_multiplier(Delta_sigma_prime_e, D_e, a, b, H);

    // Update stress and state variable increment by reference.
    Delta_sigma_prime = D_ep*Delta_epsilon_tilde_p_dT;
    delta_state = compute_plastic_state_variable_increment(Delta_epsilon_tilde_p_dT, delta_lambda, H);
}

void Elastoplastic::compute_error_estimate(void) {
    int size_state = get_state_variables().size();
    State error(size_state);
    error[0] = (to_cauchy(Delta_sigma_prime_2 - Delta_sigma_prime_1)).norm()/sigma_prime_ini.norm();
    for (int i=1; i<size_state; i++) {
        error[i] = std::abs((delta_state_2[i] - delta_state_1[i]))/state_ini[i];
    }
    R_n = 1.0/2.0*error.maxCoeff();

    // Check against machine tolerance.
    R_n = std::max(R_n, EPS);
}

double Elastoplastic::compute_elastoplastic_multiplier(Voigt Delta_sigma_prime_e, Constitutive D_e, Voigt a, Voigt b, double H) {
    return (double)(a.transpose()*Delta_sigma_prime_e)/(double)(H + a.transpose()*D_e*b);
}

Constitutive Elastoplastic::compute_elastoplastic_matrix(Constitutive D_e, Voigt a, Voigt b, double H) {
    return D_e-(D_e*b*a.transpose()*D_e)/(H+a.transpose()*D_e*b);
}

void Elastoplastic::compute_alpha_bounds(double &alpha_0, double &alpha_1) {
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
            State state_trial = get_state_variables();
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

bool Elastoplastic::check_unload_reload(Cauchy sigma_prime) {
    // Compute required derivatives for the given stress state.
    State state = get_state_variables();
    Cauchy df_dsigma_prime_check, dg_dsigma_prime_check;
    Voigt a_check, b_check;
    double dg_dp_prime_check, H_check;
    compute_derivatives(sigma_prime, state, df_dsigma_prime_check, a_check, dg_dsigma_prime_check, b_check, dg_dp_prime_check, H_check);

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

void Elastoplastic::compute_alpha(void) {
    // Current stress state (i.e. alpha = 0.0).
    State state = get_state_variables();
    double f_0 = compute_f(sigma_prime, state);

    // Fully elastic increment trial stress state.
    Cauchy sigma_prime_1 = compute_isotropic_linear_elastic_stress(sigma_prime, 1.0, Delta_epsilon_tilde);
    double f_1 = compute_f(sigma_prime_1, state);

    // Check increment type by finding alpha.
    if (f_1 <= FTOL) {
        PLOG_INFO << "Fully elastic increment; alpha = 1.0.";
        alpha = 1.0;
    } else if (std::abs(f_0) <= FTOL && f_1> FTOL) {
        PLOG_INFO << "Check potential for elastoplastic unloading-reloading.";
        if (!check_unload_reload(sigma_prime)) {
            // Compute bounds on alpha.
            double alpha_0 = 0.0;
            double alpha_1 = 1.0;
            compute_alpha_bounds(alpha_0, alpha_1);

            // Trial stress state for alpha_0.
            Cauchy sigma_prime_0 = compute_isotropic_linear_elastic_stress(sigma_prime, alpha_0, Delta_epsilon_tilde);
            double f_0 = compute_f(sigma_prime_0, state);

            // Trial stress state for alpha_1.
            Cauchy sigma_prime_1 = compute_isotropic_linear_elastic_stress(sigma_prime, alpha_1, Delta_epsilon_tilde);
            double f_1 = compute_f(sigma_prime_1, state);

            // Perform Pegasus intersection within the bounds of alpha_0 and alpha_1.
            PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = " << alpha_0 << " and alpha_1 = " << alpha_1 <<".";
            alpha = pegasus_regula_falsi(alpha_0, alpha_1, f_0, f_1);
        } else {
            PLOG_INFO << "Fully plastic increment; alpha = 0.0.";
            alpha = 0.0;
        } 
    } else if (f_0 < -FTOL && f_1 > FTOL) {
        PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = 0.0 and alpha_1 = 1.0.";
        alpha = pegasus_regula_falsi(0.0, 1.0, f_0, f_1);
    } else {
        PLOG_FATAL << "Illegal stress state.";
        assert(false);
    }
}

double Elastoplastic::pegasus_regula_falsi(double alpha_0, double alpha_1, double f_0, double f_1) {
    // Pegasus algorithm.
    double alpha_n;
    double f_n = FTOL;
    Cauchy sigma_prime_n = sigma_prime;
    State state_n = get_state_variables();

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

void Elastoplastic::compute_yield_surface(State state, int points, Eigen::VectorXd &p_prime_surface, Eigen::VectorXd &q_surface) {
    // Calculate minimum and maximum p_prime.

    // For the calculated range of p_prime, solve for q.  
}