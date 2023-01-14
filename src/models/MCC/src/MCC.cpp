#include "MCC.hpp"

MCC::MCC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state) {
    log_severity = log_severity;
    set_name("MCC");
    set_model_type("Elastoplastic");

    // Initialise logger.
    Logging::initialise_log(log_severity);

    // Check inputs.
    Checks::check_inputs(name, parameters.size(), state.size(), parameters_required, state_required);
}

State MCC::get_state_variables(void) {
    return state;
}

void MCC::set_state_variables(State new_state) {
    state = new_state;
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

State MCC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    double Delta_epsilon_vol_e = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_e));
    State elastic_state(state.size());
    elastic_state[0] = e-(e*Delta_epsilon_vol_e);
    elastic_state[1] = p_c;
    return elastic_state;
}

double MCC::compute_f(Cauchy sigma_prime, State state) {
    // State variables.
    double e = state[0];
    double p_c = state[1];

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    
    // Yield surface function.
    double f = std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c);

    // Debug output.
    PLOG_DEBUG << "f = " << f;

    return f;
}

void MCC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H) {
    // Current state variables.
    double e = state[0];
    double p_c = state[1];

    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    Cauchy sigma = compute_sigma(sigma_prime, u);
    Cauchy s = compute_s(sigma, p);
    double p = compute_p(sigma);
    double I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar;
    compute_stress_invariants(sigma, p, s, I_1, I_2, I_3, J_1, J_2, J_3);
    compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
    Cauchy dq_dsigma_prime = compute_dq_dsigma_prime(sigma_prime, s, q);
    Cauchy dtheta_dsigma_prime = compute_dtheta_dsigma_prime(sigma_prime);
    
    // Compute derivatives of the yield surface with respect to the stress state.
    double df_dq = 2*q;
    double df_dp_prime = std::pow(M,2)*(2*p_prime-p_c);
    double df_dtheta = 0.0;
    df_dsigma_prime = df_dq*dq_dsigma_prime + df_dp_prime*dp_prime_dsigma_prime + df_dtheta*dtheta_dsigma_prime;

    // Compute derivatives of the plastic potential function with respect to the stress state.
    double dg_dq = 2*q;
    double dg_dp_prime = std::pow(M,2)*(2*p_prime-p_c);
    double dg_dtheta = 0.0;
    dg_dsigma_prime = dg_dq*dq_dsigma_prime + dg_dp_prime*dp_prime_dsigma_prime + dg_dtheta*dtheta_dsigma_prime;
    
    // Convert to Voigt notation.
    a = to_voigt(df_dsigma_prime);
    b = to_voigt(dg_dsigma_prime);

    // Hardening modulus.
    H = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*df_dsigma_prime.trace();

    // Debug output.
    PLOG_DEBUG << "a = " << a;
    PLOG_DEBUG << "b = " << b;
    PLOG_DEBUG << "H = " << H;
}

State MCC::compute_plastic_state_variable_increment(Cauchy sigma_prime, Voigt Delta_epsilon_tilde_p, double delta_lambda, double H) {
    // Current state variables.
    double e = state[0];
    double p_c = state[1];

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);

    // Calculate increment in plastic state variables.
    double Delta_epsilon_vol_p = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_p));
    State delta_state(state.size());
    delta_state[0] = -(1+e)*Delta_epsilon_vol_p;
    delta_state[1] = delta_lambda*H/(std::pow(M,2)*p_prime);
    return delta_state;
}

State MCC::compute_plastic_state_variable_increment(Cauchy sigma_prime, double delta_lambda, double H) {
    // Note: only correct state variables that do not depend on the magnitude of the strain increment (hence strain increment is not passed in).
    Voigt Delta_epsilon_tilde_p = Voigt::Zero();
    return compute_plastic_state_variable_increment(sigma_prime, Delta_epsilon_tilde_p, delta_lambda, H);
}
