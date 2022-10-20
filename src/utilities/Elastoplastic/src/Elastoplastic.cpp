#include "Elastoplastic.hpp"

void Elastoplastic::solve(void) {
    // Check increment type.
    double alpha;

    // Current stress state.
    double f_0 = compute_f(sigma_prime);

    // Fully elastic increment trial stress state.
    Cauchy sigma_prime_1 = compute_isotropic_linear_elastic_trial_stress(sigma_prime, 1.0, delta_epsilon_tilde);
    double f_1 = compute_f(sigma_prime_1);

    // Check increment type.
    if (f_1 <= FTOL) {
        PLOG_INFO << "Fully elastic increment; alpha = 1.0.";
        alpha = 1.0;
    } else if (std::abs(f_0) <= FTOL and f_1> FTOL) {
        PLOG_INFO << "Check potential for elastoplastic unloading-reloading.";
    } else if (f_0 < -FTOL && f_1 > FTOL) {
        PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm.";
        alpha = compute_alpha(0.0, 1.0, f_0, f_1);
    } else {
        PLOG_FATAL << "Illegal stress state.";
        assert(false);
    }
    std::cout << "alpha = " << alpha << "\n";
}

double Elastoplastic::compute_alpha(double alpha_0, double alpha_1, double f_0, double f_1) {
    // Pegasus algorithm.
    int iterations = 0;
    double alpha_n;
    double f_n = FTOL;
    Cauchy sigma_prime_n = sigma_prime;

    // Iterate to find optimal alpha if a plastic increment.
    while (iterations < MAXITS && std::abs(f_n) >= FTOL) {
        alpha_n = alpha_1 - f_1*(alpha_1-alpha_0)/(f_1-f_0);  
            
        sigma_prime_n = compute_isotropic_linear_elastic_trial_stress(sigma_prime, alpha_n, delta_epsilon_tilde);
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