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
    double f = std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c);
    return f;
}

double MCC::compute_K(double delta_epsilon_vol, double p_prime) {
    if (delta_epsilon_vol != 0.0) {
        K = (p_prime/delta_epsilon_vol)*(std::exp(delta_epsilon_vol/kappa_star)-1);
    } else {
        K = p_prime/kappa_star;
    }
    return K;
}

double MCC::compute_G(double K) {
    G = (3.0*(1-2.0*nu)*K)/(2.0*(1.0+nu));
    return G;
}

