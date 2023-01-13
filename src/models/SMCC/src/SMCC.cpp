#include "SMCC.hpp"

SMCC::SMCC(Parameters parameters, State state) : parameters(parameters), state(state) {   
    set_name("SMCC");
    int parameters_required = 8;
    int state_required = 3;
    PLOG_FATAL_IF(parameters.size() != parameters_required) << parameters.size() << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state.size() != state_required) << state.size() << " parameters supplied when " << state_required << " expected.";
    assert(parameters.size() == parameters_required && state.size() == state_required);
    PLOG_INFO << name << " model instantiated with " << parameters.size() << " parameters and " << state.size() << " state variables.";  
}

State SMCC::get_state_variables(void) {
    return state;
}

void SMCC::set_state_variables(State new_state) {
    state = new_state;
}

State SMCC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    double Delta_epsilon_vol_e = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_e));
    State elastic_state(state.size());
    elastic_state[0] = e-(e*Delta_epsilon_vol_e);
    elastic_state[1] = p_c;
    elastic_state[2] = s_ep;
    return elastic_state;
}

State SMCC::compute_plastic_state_variable_increment(Voigt Delta_epsilon_tilde_p, double delta_lambda, double H) {
    double Delta_epsilon_vol_p = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_p));
    State delta_state(state.size());
    delta_state[0] = -(1+e)*Delta_epsilon_vol_p;
    delta_state[1] = delta_lambda*(-((std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star))*s_ep*df_dsigma_prime.trace()/-(std::pow(M,2)*p_prime*s_ep));
    delta_state[2] = delta_lambda*(-((std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*-k*(s_ep-1.0)*std::sqrt((1-A)*std::pow(df_dsigma_prime.trace(),2) + (A*2.0/3.0*(df_dsigma_prime*(df_dsigma_prime.transpose())).trace())))/-(std::pow(M,2)*p_prime*p_c));
    return delta_state;
}

State SMCC::compute_plastic_state_variable_increment(double delta_lambda, double H) {
    // Note: only correct state variables that do not depend on the magnitude of the strain increment (hence strain increment is not passed in).
    Voigt Delta_epsilon_tilde_p = Voigt::Zero();
    return compute_plastic_state_variable_increment(Delta_epsilon_tilde_p, delta_lambda, H);
}

double SMCC::compute_K(double Delta_epsilon_e_vol, double p_prime) {
    if (Delta_epsilon_e_vol != 0.0) {
        // Return secant bulk modulus.
        return (p_prime/Delta_epsilon_e_vol)*(std::exp(Delta_epsilon_e_vol/kappa_star)-1);
    } else {
        // Return tangent bulk modulus.
        return p_prime/kappa_star;
    }
}

double SMCC::compute_G(double K) {
    double G = (3.0*(1-2.0*nu)*K)/(2.0*(1.0+nu));
    return G;
}

double SMCC::compute_f(Cauchy sigma_prime, State state) {
    // State variables.
    double e = state[0];
    double p_c = state[1];

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    
    // Yield surface function.
    return std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c*s_ep);
}

double SMCC::compute_df_dq(void) {
    double df_dq = 2*q;
    return df_dq;
}

double SMCC::compute_df_dp_prime(void) {
    double df_dp_prime = std::pow(M,2)*(2*p_prime-p_c*s_ep);
    return df_dp_prime;
}

double SMCC::compute_df_dtheta(void) {
    df_dtheta = 0.0;
    return df_dtheta;
}

double SMCC::compute_dg_dq(void) {
    double dg_dq = 2*q;
    return dg_dq;
}

double SMCC::compute_dg_dp_prime(void) {
    double dg_dp_prime = std::pow(M,2)*(2*p_prime-p_c*s_ep);
    return dg_dp_prime;
}

double SMCC::compute_dg_dtheta(void) {
    dg_dtheta = 0.0;
    return dg_dtheta;
}

double SMCC::compute_H(void) {
    double H_p_c = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*s_ep*df_dsigma_prime.trace();
    double H_s_ep = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*-k*(s_ep-1.0)*std::sqrt((1-A)*std::pow(df_dsigma_prime.trace(),2) + (A*2.0/3.0*(df_dsigma_prime*(df_dsigma_prime.transpose())).trace()));
    double H = H_p_c + H_s_ep;
    return H;
}
