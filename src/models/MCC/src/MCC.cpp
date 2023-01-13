#include "MCC.hpp"

MCC::MCC(Parameters parameters, State state) : parameters(parameters), state(state) {   
    set_name("MCC");
    int parameters_required = 5;
    int state_required = 2;
    PLOG_FATAL_IF(parameters.size() != parameters_required) << parameters.size() << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state.size() != state_required) << state.size() << " parameters supplied when " << state_required << " expected.";
    assert(parameters.size() == parameters_required && state.size() == state_required);
    PLOG_INFO << name << " model instantiated with " << parameters.size() << " parameters and " << state.size() << " state variables.";  
}

State MCC::get_state_variables(void) {
    return state;
}

void MCC::set_state_variables(State new_state) {
    state = new_state;
}

State MCC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    double Delta_epsilon_vol_e = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_e));
    State elastic_state(state.size());
    elastic_state[0] = e-(e*Delta_epsilon_vol_e);
    elastic_state[1] = p_c;
    return elastic_state;
}

State MCC::compute_plastic_state_variable_increment(Voigt Delta_epsilon_tilde_p, double delta_lambda, double H) {
    double Delta_epsilon_vol_p = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_p));
    State delta_state(state.size());
    delta_state[0] = -(1+e)*Delta_epsilon_vol_p;
    delta_state[1] = delta_lambda*H/(std::pow(M,2)*p_prime);
    return delta_state;
}

State MCC::compute_plastic_state_variable_increment(double delta_lambda, double H) {
    // Note: only correct state variables that do not depend on the magnitude of the strain increment (hence strain increment is not passed in).
    Voigt Delta_epsilon_tilde_p = Voigt::Zero();
    return compute_plastic_state_variable_increment(Delta_epsilon_tilde_p, delta_lambda, H);
}

double MCC::compute_K(double Delta_epsilon_e_vol, double p_prime) {
    if (Delta_epsilon_e_vol != 0.0) {
        // Return secant bulk modulus.
        return (p_prime/Delta_epsilon_e_vol)*(std::exp(Delta_epsilon_e_vol/kappa_star)-1);
    } else {
        // Return tangent bulk modulus.
        return p_prime/kappa_star;
    }
}

double MCC::compute_G(double K) {
    double G = (3.0*(1-2.0*nu)*K)/(2.0*(1.0+nu));
    return G;
}

double MCC::compute_f(Cauchy sigma_prime, State state) {
    // State variables.
    double e = state[0];
    double p_c = state[1];

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    
    // Yield surface function.
    return std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c);
}

double MCC::compute_df_dq(void) {
    double df_dq = 2*q;
    return df_dq;
}

double MCC::compute_df_dp_prime(void) {
    double df_dp_prime = std::pow(M,2)*(2*p_prime-p_c);
    return df_dp_prime;
}

double MCC::compute_df_dtheta(void) {
    df_dtheta = 0.0;
    return df_dtheta;
}

double MCC::compute_dg_dq(void) {
    double dg_dq = 2*q;
    return dg_dq;
}

double MCC::compute_dg_dp_prime(void) {
    double dg_dp_prime = std::pow(M,2)*(2*p_prime-p_c);
    return dg_dp_prime;
}

double MCC::compute_dg_dtheta(void) {
    dg_dtheta = 0.0;
    return dg_dtheta;
}

double MCC::compute_H(void) {
    double H = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*df_dsigma_prime.trace();
    return H;
}
