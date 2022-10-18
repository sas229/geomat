#include "Elastoplastic.hpp"

void Elastoplastic::solve(void) {
    std::cout << "Solving current strain increment.\n";
    alpha = compute_alpha();
}

double Elastoplastic::compute_alpha(void) {
    // Pegasus algorithm.
    int iterations = 0;
    int max_iterations = 100;
    double tolerance = 1e-10;
    double alpha;
    double alpha_n;
    double alpha_0 = 0.0;
    double alpha_1 = 1.0;
    Cauchy sigma_prime_0 = sigma_prime;
    Cauchy sigma_prime_1 = sigma_prime;
    Cauchy sigma_prime_n = sigma_prime;
    double f_0, f_1, f_n = tolerance;

    // First trial state: no strain increment.
    sigma_prime_0 = compute_isotropic_linear_elastic_trial_stress(sigma_prime, alpha_0, delta_epsilon_tilde);
    f_0 = compute_f(sigma_prime_0);

    // Second trial state: fully elastic strain increment.
    sigma_prime_1 = compute_isotropic_linear_elastic_trial_stress(sigma_prime, alpha_1, delta_epsilon_tilde);
    f_1 = compute_f(sigma_prime_1);

    // Check if plastic increment.
    if (f_1 > 0.0) {
        PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm.";
    } else {
        PLOG_INFO << "Fully elastic increment.";
        alpha = 1.0;
        return alpha;
    }

    // Iterate to find optimal alpha if a plastic increment.
    while (iterations < max_iterations && std::abs(f_n) >= tolerance) {
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
    alpha = alpha_n;
    PLOG_INFO << "Pegasus method iterations = " << iterations << "; " << "alpha = " << alpha << "; " << "f_n = " << f_n << ".";
    return alpha;
}