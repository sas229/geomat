#include "MCC.hpp"

MCC::MCC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state), Elastoplastic::Elastoplastic() {   
    set_model_name("MCC");
    set_model_type("Elastoplastic");

    // Initialise logger.
    Logging::initialise_log(log_severity);

    // Check inputs.
    Checks::check_inputs(name, (int)parameters.size(), (int)state.size(), parameters_required, state_required);
}

State MCC::get_state_variables(void) {
    return state;
}

void MCC::set_state_variables(State new_state) {
    state = new_state;
}

Constitutive MCC::compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon) {
    double Delta_epsilon_e_vol = Delta_epsilon.trace();
    double p_prime = compute_p_prime(sigma_prime);
    double K;
    if (Delta_epsilon_e_vol != 0.0) {
        // Secant bulk modulus.
        K = (p_prime/Delta_epsilon_e_vol)*(std::exp(Delta_epsilon_e_vol/kappa_star)-1);
    } else {
        // Tangent bulk modulus.
        K = p_prime/kappa_star;
    }
    double G = (3.0*(1-2.0*nu)*K)/(2.0*(1.0+nu));
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K, G);  
    return D_e;
}

double MCC::compute_f(Cauchy sigma_prime, State state) {
    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    
    // State variables.
    double e = state[0];
    double p_c = state[1];
    
    // Yield surface function.
    double f = std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c);

    return f;
}

void MCC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, HardeningModuli  &H_s) {
    // State variables.
    double e = state[0];
    double p_c = state[1];

    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double q, p_prime;
    Cauchy s, dq_dsigma_prime;
    q = compute_q(sigma_prime);
    p_prime = compute_p_prime(sigma_prime);
    s = compute_s(sigma_prime, p_prime);
    dq_dsigma_prime = compute_dq_dsigma_prime(sigma_prime, s, q);
    
    // Yield surface derivatives.
    double df_dq = 2*q;
    double df_dp_prime = std::pow(M,2)*(2*p_prime-p_c);
    double df_dtheta = 0;

    // Derivatives of yield surface and plastic potential function.
    df_dsigma_prime = (df_dp_prime*dp_prime_dsigma_prime) + (df_dq*dq_dsigma_prime);
    dg_dsigma_prime = df_dsigma_prime; // Associated flow.
    
    // Vectors of derivatives.
    a = to_voigt(df_dsigma_prime);
    b = to_voigt(dg_dsigma_prime);

    // Hardening modulus.
    H = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*df_dsigma_prime.trace();
}

State MCC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    double Delta_epsilon_vol_e = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_e));
    State elastic_state(state.size());
    elastic_state[0] = e-(e*Delta_epsilon_vol_e);
    elastic_state[1] = p_c;
    return elastic_state;
}

State MCC::compute_plastic_state_variable_increment(double delta_lambda, Cauchy df_dsigma_prime, HardeningModuli  H_s, Voigt Delta_epsilon_tilde_p) {
    double Delta_epsilon_vol_p = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_p));
    State delta_state(state.size());
    delta_state[0] = -(1+e)*Delta_epsilon_vol_p;
    delta_state[1] = delta_lambda*H/(pow(M,2)*p_prime);
    return delta_state;
}