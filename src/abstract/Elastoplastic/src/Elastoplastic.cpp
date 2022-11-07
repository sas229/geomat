#include "Elastoplastic.hpp"

void Elastoplastic::solve(void) { 
    compute_alpha();
    std::cout << "alpha = " << alpha << "\n";

    // Update stress and state variables based on elastic portion of strain increment.
    Voigt delta_epsilon_tilde_e = alpha*delta_epsilon_tilde;
    Cauchy sigma_prime_e = compute_isotropic_linear_elastic_stress(sigma_prime, alpha, delta_epsilon_tilde_e);
    delta_epsilon_vol_e = compute_delta_epsilon_vol(delta_epsilon_tilde_e.cauchy());
    compute_elastic_state_variable_update();

    // Perform elastoplastic stress integration on plastic portion of strain increment.
    if (alpha < 1.0) {
        delta_epsilon_tilde_p = (1.0-alpha)*delta_epsilon_tilde;
    }

    D_e = compute_isotropic_linear_elastic_matrix(50000, 25000);
    compute_derivatives(sigma_prime, df_dsigma_prime, a, dg_dsigma_prime, b, dg_dp_prime, H);
    D_ep = D_e-(D_e*b*a.transpose()*D_e)/(H+a.transpose()*D_e*b);
    std::cout << D_ep << "\n";
}

void Elastoplastic::compute_alpha_bounds(double &alpha_0, double &alpha_1) {
    double delta_epsilon_vol_e = delta_epsilon.trace();
    double alpha_n, d_alpha;
    int i = 0;
    int j = 0;
    while (i < MAXITS) {
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
            Cauchy sigma_prime_trial = sigma_prime + delta_sigma_e_trial.cauchy();
            double f_trial = compute_f(sigma_prime_trial);

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
    Cauchy df_dsigma_prime_check, dg_dsigma_prime_check;
    Voigt a_check, b_check;
    double dg_dp_prime_check, H_check;
    compute_derivatives(sigma_prime, df_dsigma_prime_check, a_check, dg_dsigma_prime_check, b_check, dg_dp_prime_check, H_check);

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
    double f_0 = compute_f(sigma_prime);

    // Fully elastic increment trial stress state.
    Cauchy sigma_prime_1 = compute_isotropic_linear_elastic_stress(sigma_prime, 1.0, delta_epsilon_tilde);
    double f_1 = compute_f(sigma_prime_1);

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
            double f_0 = compute_f(sigma_prime_0);

            // Trial stress state for alpha_1.
            Cauchy sigma_prime_1 = compute_isotropic_linear_elastic_stress(sigma_prime, alpha_1, delta_epsilon_tilde);
            double f_1 = compute_f(sigma_prime_1);

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

    // Iterate to find optimal alpha if a plastic increment.
    while (iterations < MAXITS && std::abs(f_n) >= FTOL) {
        alpha_n = alpha_1 - f_1*(alpha_1-alpha_0)/(f_1-f_0);  
            
        sigma_prime_n = compute_isotropic_linear_elastic_stress(sigma_prime, alpha_n, delta_epsilon_tilde);
        f_n = compute_f(sigma_prime_n);

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
        PLOG_FATAL << "Performed " << MAXITS << " Pegasus iterations: alpha = " << alpha << "; |f| = " << std::abs(f) << " > tolerance = " << FTOL << ".";
        assert(false);
    } else{
        PLOG_INFO << "Performed " << iterations << " Pegasus iterations: " << "alpha = " << alpha << "; " << "|f| = " << std::abs(f) << " < tolerance = " << FTOL << ".";
        assert(alpha >= 0.0 && alpha <= 1.0);
    }
    return alpha;
}