#include "Elastoplastic.hpp"

void Elastoplastic::solve(void) { 
    compute_alpha();
    int substeps = 0;
    int drift_corrections = 0;
    double R_n;

    if (alpha > 0.0) {
        // Update stress and state variables based on elastic portion of strain increment.
        delta_epsilon_tilde_e = alpha*delta_epsilon_tilde;
        delta_sigma_prime_e = compute_elastic_stress_increment(D_e, delta_epsilon_tilde_e);
        sigma_prime_e = compute_isotropic_linear_elastic_stress(sigma_prime, alpha, delta_epsilon_tilde_e);
        delta_epsilon_vol_e = compute_delta_epsilon_vol(delta_epsilon_tilde_e.cauchy());
        compute_elastic_state_variable();
    } 
    if (alpha < 1.0) {
        // Perform elastoplastic stress integration on plastic portion of strain increment.
        delta_epsilon_tilde_p = (1.0-alpha)*delta_epsilon_tilde;
        Cauchy sigma_prime_ep = sigma_prime_e;
        Eigen::VectorXd state_ep = get_state_variables();

        // Substepping with automatic error control.
        dT = 1.0;
        solved = false;
        T = 0.0;
        double q_step;
        while (T < 1.0) {
            delta_epsilon_tilde_p_dT = delta_epsilon_tilde_p*dT;
            delta_epsilon_vol_p_dT = (delta_epsilon_tilde_p_dT.cauchy()).trace();

            // Compute stress and state variable increment estimates.
            Voigt delta_sigma_prime_1, delta_sigma_prime_2;
            Eigen::VectorXd delta_state_1(2), delta_state_2(2);
            compute_plastic_increment(sigma_prime_ep, state_ep, delta_epsilon_tilde_p_dT, delta_sigma_prime_1, delta_state_1);
            sigma_prime_1 = sigma_prime_ep + delta_sigma_prime_1.cauchy();
            compute_plastic_increment(sigma_prime_1, state_ep, delta_epsilon_tilde_p_dT, delta_sigma_prime_2, delta_state_2);

            // Calculate modified Euler stresses and state variables.
            Cauchy sigma_prime_ini = sigma_prime_ep + 1.0/2.0*(delta_sigma_prime_1.cauchy() + delta_sigma_prime_2.cauchy());
            Eigen::VectorXd state_ini = state_ep + 1.0/2.0*(delta_state_1+delta_state_2);

            // Compute error estimate.
            R_n = compute_error_estimate(sigma_prime_ini, delta_sigma_prime_1, delta_sigma_prime_2, state_ini, delta_state_1, delta_state_2);
            if (R_n < STOL) {
                // Accept increment and correct stresses and state variables back to yield surface.
                Cauchy sigma_prime_u = sigma_prime_ini;
                Cauchy sigma_prime_c = sigma_prime_u; 
                Eigen::VectorXd state_u = state_ini;
                Eigen::VectorXd state_c = state_u;

                // Correct stresses and state variables back to yield surface. 
                double f_u = compute_f(sigma_prime_u, state_u);
                double f_c = f_u;
                int ITS_YSC = 0;
                while (ITS_YSC < MAXITS_YSC) {
                    if (std::abs(f_c) > FTOL) {
                        // Calculate elastic constitutive matrix using tangent moduli and elastic stress increment.
                        double p_prime_u = compute_p_prime(sigma_prime_u);
                        double K_tan_u = compute_K(0, p_prime_u);
                        double G_tan_u = compute_G(K_tan_u);
                        Constitutive D_e_u = compute_isotropic_linear_elastic_matrix(K_tan_u, G_tan_u);

                        // Calculate derivatives.
                        Cauchy df_dsigma_prime_u, dg_dsigma_prime_u;
                        Voigt a_u, b_u;
                        double dg_dp_prime_u, H_u;
                        compute_derivatives(sigma_prime_u, state_u, df_dsigma_prime_u, a_u, dg_dsigma_prime_u, b_u, dg_dp_prime_u, H_u);

                        // Compute correction factor.
                        double delta_lambda_c = f_u/(H_u + a_u.transpose()*D_e_u*b_u);
                        Voigt delta_sigma_prime_c = -delta_lambda_c*D_e_u*b_u;
                        Eigen::VectorXd delta_state_c = compute_plastic_state_variable_correction(delta_lambda_c, H_u);

                        // Update stress and state variables using corrections.
                        sigma_prime_c = sigma_prime_u + delta_sigma_prime_c.cauchy();
                        state_c = state_u + delta_state_c;

                        // Check yield surface function.
                        f_c = compute_f(sigma_prime_c, state_c);
                        f_u = compute_f(sigma_prime_u, state_u);
                        if (std::abs(f_c) > std::abs(f_u)) {
                            // Apply normal correction instead.
                            delta_lambda_c = f_u/(a.transpose()*a);
                            delta_sigma_prime_c = -delta_lambda_c*a;
                            sigma_prime_c = sigma_prime_u + delta_sigma_prime_c.cauchy();
                            state_c = state_u;
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
                CORRECTED: drift_corrections += ITS_YSC;
                solved = true;
                sigma_prime_ep = sigma_prime_c;
                state_ep = state_c;

                // Increment substep and pseudotime T; reset solved boolean.
                T += dT;
                substeps += 1;
                solved = false;

                // Compute new step size factor.
                q_step = std::min(0.9*std::sqrt(STOL/R_n),1.1);
            } else {
                // Step size rejected: calculate reduced step size factor and limit step size growth factor..
                q_step = std::max(0.9*std::sqrt(STOL/R_n),0.1);
                q_step = std::min(q_step, 1.0);
            }
            // Update pseudo-time T and check (i) that the next step size is not smaller than the minimum size specified and (ii) does not exceed T = 1.
            dT = q_step*dT;
            dT = std::max(dT, dT_min);
            dT = std::min(dT, 1.0-T);
        }
        // Set stress state as final elastoplastic stress state.
        sigma_prime = sigma_prime_ep;
        PLOG_INFO << "Elastoplastic stress increment integrated to within a tolerance FTOL = " << FTOL << " via " << substeps << " substeps and " << drift_corrections << " drift corrections.";
    } else {
        // Fully elastic increment. Update stress state.
        sigma_prime = sigma_prime_e;
        PLOG_INFO << "Fully elastic tress increment integrated directly.";
    }
    
    // Compute final stress invariants.
    compute_principal_stresses(sigma_prime, sigma_1, sigma_2, sigma_3, R, S);
    std::cout << "Principal stresses: sigma_1 = " << sigma_1 << "; sigma_2 = " << sigma_2 << "; sigma_3 = " << sigma_3 << "\n";
}

void Elastoplastic::compute_plastic_increment(Cauchy sigma_prime, Eigen::VectorXd state, Voigt delta_epsilon_tilde_p_dT, Voigt &delta_sigma_prime, Eigen::VectorXd &delta_state) {   
    // Calculate elastic constitutive matrix using tangent moduli and elastic stress increment.
    double p_prime = compute_p_prime(sigma_prime);
    double K_tan = compute_K(0, p_prime);
    double G_tan = compute_G(K_tan);
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K_tan, G_tan);
    Voigt delta_sigma_prime_e = D_e*delta_epsilon_tilde_p_dT;

    // Compute elastoplastic constitutive matrix and elastoplastic multiplier. 
    Cauchy df_dsigma_prime, dg_dsigma_prime;
    Voigt a, b;
    double dg_dp_prime, H;
    compute_derivatives(sigma_prime, state, df_dsigma_prime, a, dg_dsigma_prime, b, dg_dp_prime, H);
    Constitutive D_ep = compute_elastoplastic_matrix(D_e, a, b, H);
    double delta_lambda = compute_elastoplastic_multiplier(delta_sigma_prime_e, D_e, a, b, H);

    // Update stress and state variable increment by reference.
    delta_sigma_prime = D_ep*delta_epsilon_tilde_p_dT;
    delta_state = compute_plastic_state_variable_increment(delta_lambda, H);
}

double Elastoplastic::compute_error_estimate(Cauchy sigma_prime_ini, Voigt delta_sigma_prime_1, Voigt delta_sigma_prime_2, Eigen::VectorXd state_ini, Eigen::VectorXd delta_state_1, Eigen::VectorXd delta_state_2) {
    Eigen::VectorXd error(state_ini.size()+1);
    error[0] = ((delta_sigma_prime_2.cauchy() - delta_sigma_prime_1.cauchy())).norm()/sigma_prime_ini.norm();
    for (int i=1; i<state_ini.size(); i++) {
        error[i] = std::abs((delta_state_2[i] - delta_state_1[i]))/state_ini[i];
    }
    double R_n = 1.0/2.0*error.maxCoeff();

    // Check against machine tolerance.
    R_n = std::max(R_n, EPS);
    return R_n;
}

double Elastoplastic::compute_elastoplastic_multiplier(Voigt delta_sigma_prime_e, Constitutive D_e, Voigt a, Voigt b, double H) {
    double numerator = (a.transpose()*delta_sigma_prime_e);
    double denominator = H + a.transpose()*D_e*b;
    return numerator/denominator;
}

Constitutive Elastoplastic::compute_elastoplastic_matrix(Constitutive D_e, Voigt a, Voigt b, double H) {
    return D_e-(D_e*b*a.transpose()*D_e)/(H+a.transpose()*D_e*b);
}

void Elastoplastic::compute_alpha_bounds(double &alpha_0, double &alpha_1) {
    double delta_epsilon_vol_e = delta_epsilon.trace();
    double alpha_n, d_alpha;
    int i = 0;
    int j = 0;
    while (i < MAXITS_YSI) {
        d_alpha = (alpha_1-alpha_0)/NSUB;
        alpha_n = alpha_0+d_alpha;
        while (j < NSUB) {
            // Compute the elastic matrix.
            double K_trial = compute_K(alpha_n*delta_epsilon_vol_e, p_prime);
            double G_trial = compute_G(K_trial);
            Constitutive D_e_trial = compute_isotropic_linear_elastic_matrix(K_trial, G_trial);

            // Compute elastic stress increment.
            Voigt delta_sigma_e_trial = alpha_n*D_e_trial*delta_epsilon_tilde;

            // Check yield function.
            Eigen::VectorXd state_trial = get_state_variables();
            Cauchy sigma_prime_trial = sigma_prime + delta_sigma_e_trial.cauchy();
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
    Eigen::VectorXd state = get_state_variables();
    Cauchy df_dsigma_prime_check, dg_dsigma_prime_check;
    Voigt a_check, b_check;
    double dg_dp_prime_check, H_check;
    compute_derivatives(sigma_prime, state, df_dsigma_prime_check, a_check, dg_dsigma_prime_check, b_check, dg_dp_prime_check, H_check);

    // Compute the elastic matrix using tangent moduli.
    double K_tan = compute_K(0, p_prime);
    double G_tan = compute_G(K_tan);
    Constitutive D_e_tan = compute_isotropic_linear_elastic_matrix(K_tan, G_tan);

    // Compute elastic stress increment.
    Voigt delta_sigma_e = D_e_tan*delta_epsilon_tilde;

    // Check unloading-reloading criterion.
    double norm_a = a.squaredNorm();
    double norm_delta_sigma_e = delta_sigma_e.squaredNorm();
    double numerator = a.transpose()*delta_sigma_e;       
    double denominator = norm_a*norm_delta_sigma_e;
    double cos_theta = numerator/denominator;

    // Check against tolerance.
    if (cos_theta >= -LTOL) {
        return false;
    } else {
        return true;
    }
}

void Elastoplastic::compute_alpha(void) {
    // Current stress state (i.e. alpha = 0.0).
    Eigen::VectorXd state = get_state_variables();
    double f_0 = compute_f(sigma_prime, state);

    // Fully elastic increment trial stress state.
    Cauchy sigma_prime_1 = compute_isotropic_linear_elastic_stress(sigma_prime, 1.0, delta_epsilon_tilde);
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
            Cauchy sigma_prime_0 = compute_isotropic_linear_elastic_stress(sigma_prime, alpha_0, delta_epsilon_tilde);
            double f_0 = compute_f(sigma_prime_0, state);

            // Trial stress state for alpha_1.
            Cauchy sigma_prime_1 = compute_isotropic_linear_elastic_stress(sigma_prime, alpha_1, delta_epsilon_tilde);
            double f_1 = compute_f(sigma_prime_1, state);

            // Perform Pegasus intersection within the bounds of alpha_0 and alpha_1.
            PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = " << alpha_0 << " and alpha_1 = " << alpha_1 <<".";
            alpha = pegasus_regula_falsi(alpha_0, alpha_1, f_0, f_1);
        } else {
            PLOG_INFO << "Fully elastic increment; alpha = 0.0.";
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
    int iterations = 0;
    double alpha_n;
    double f_n = FTOL;
    Cauchy sigma_prime_n = sigma_prime;
    Eigen::VectorXd state_n = get_state_variables();

    // Iterate to find optimal alpha if a plastic increment.
    while (iterations < MAXITS_YSI && std::abs(f_n) >= FTOL) {
        alpha_n = alpha_1 - f_1*(alpha_1-alpha_0)/(f_1-f_0);  
            
        sigma_prime_n = compute_isotropic_linear_elastic_stress(sigma_prime, alpha_n, delta_epsilon_tilde);
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
        iterations += 1;
    }
    double alpha = alpha_n;
    double f = f_n;
    if (std::abs(f) >= FTOL) {
        PLOG_FATAL << "Performed " << MAXITS_YSI << " Pegasus iteration(s): alpha = " << alpha << "; |f| = " << std::abs(f) << " > tolerance = " << FTOL << ".";
        assert(false);
    } else{
        PLOG_INFO << "Performed " << iterations << " Pegasus iteration(s): " << "alpha = " << alpha << "; " << "|f| = " << std::abs(f) << " < tolerance = " << FTOL << ".";
        assert(alpha >= 0.0 && alpha <= 1.0);
    }
    return alpha;
}