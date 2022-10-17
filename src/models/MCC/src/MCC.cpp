#include "MCC.hpp"

MCC::MCC(std::vector<double> parameters, std::vector<double> state) : parameters(parameters), state(state) {
    set_name("MCC");
    int parameters_required = 5;
    int state_required = 2;
    PLOG_FATAL_IF(parameters.size() != parameters_required) << parameters.size() << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state.size() != state_required) << state.size() << " parameters supplied when " << state_required << " expected.";
    assert(parameters.size() == parameters_required && state.size() == state_required);
    PLOG_INFO << name << " model instantiated with " << parameters.size() << " parameters and " << state.size() << " state variables.";    
}

double MCC::compute_f(Cauchy sigma_prime) {
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    return std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c);
}

void MCC::solve(void) {
    // Pegasus algorithm.
    int iterations = 0;
    int max_iterations = 100;
    double tolerance = 1e-10;
    double alpha_n;
    double alpha_0 = 0.0;
    double alpha_1 = 1.0;
    Cauchy sigma_prime_0 = sigma_prime;
    Cauchy sigma_prime_1 = sigma_prime;
    Cauchy sigma_prime_n = sigma_prime;
    double f_0, f_1, f_n = tolerance;

    // First trial state: no strain increment.
    sigma_prime_0 = compute_elastic_trial_stress(sigma_prime, alpha_0, delta_epsilon_tilde, kappa_star, nu);
    f_0 = compute_f(sigma_prime_0);

    // Second trial state: fully elastic strain increment.
    sigma_prime_1 = compute_elastic_trial_stress(sigma_prime, alpha_1, delta_epsilon_tilde, kappa_star, nu);
    f_1 = compute_f(sigma_prime_1);

    if (f_1 > 0.0) {
        std::cout << "Plastic increment: find alpha via the Pegasus algorithm...\n";
    }
    // Iterate to find optimal alpha.
    while (iterations < max_iterations && std::abs(f_n) >= tolerance) {
        alpha_n = alpha_1 - f_1*(alpha_1-alpha_0)/(f_1-f_0);  
             
        sigma_prime_n = compute_elastic_trial_stress(sigma_prime, alpha_n, delta_epsilon_tilde, kappa_star, nu);
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
    std::cout << "Pegasus method iterations = " << iterations << "\n";
    std::cout << "f_n = " << f_n << "\n";
    std::cout << "alpha = " << alpha_n << "\n";
}

