#include "SMCC.hpp"
#include "SMCC_Definition.hpp"

SMCC::SMCC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state) {   
    set_model_name("SMCC");
    set_model_type("Elastoplastic");

    // Initialise logger.
    Logging::initialise_log(log_severity);

    // Check inputs.
    Checks::check_inputs(name, parameters.size(), state.size(), parameters_required, state_required);
}

double SMCC::compute_f(Cauchy sigma_prime, State state) {

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);

    // State variables.
    double e = state[0];
    double p_c = state[1];
    double s_ep = state[2];

    return SMCC_YIELD;
}

Cauchy SMCC::compute_elastic_stress(Cauchy sigma_prime, double alpha, Voigt Delta_epsilon_tilde) {
    return compute_isotropic_linear_elastic_stress(sigma_prime, alpha, Delta_epsilon_tilde);
}

Constitutive SMCC::compute_elastic_matrix(Cauchy sigma_prime, double Delta_epsilon_e_vol) {
    double p_prime = compute_p_prime(sigma_prime);
    double K;
    if (Delta_epsilon_e_vol != 0.0) {
        // Return secant bulk modulus.
        K = SMCC_SECANT_BULK_MODULUS;
    } else {
        // Return tangent bulk modulus.
        K = SMCC_TANGENT_BULK_MODULUS;
    }
    double G = SMCC_SHEAR_MODULUS;
    D_e = compute_isotropic_linear_elastic_matrix(K, G);
    return D_e;
}

double SMCC::compute_G(double K) {
    return SMCC_SHEAR_MODULUS;
}

double SMCC::compute_K(double Delta_epsilon_e_vol, double p_prime) {
    if (Delta_epsilon_e_vol != 0.0) {
        // Return secant bulk modulus.
        return SMCC_SECANT_BULK_MODULUS;
    } else {
        // Return tangent bulk modulus.
        return SMCC_TANGENT_BULK_MODULUS;
    }
}

State SMCC::get_state_variables(void) {
    return state;
}

void SMCC::set_state_variables(State new_state) {
    state = new_state;
}

void SMCC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H) {
    // State variables.
    double e = state[0];
    double p_c = state[1];
    double s_ep = state[2];

    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double q, p_prime, I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar;
    Cauchy s, dq_dsigma_prime, dJ_3_dsigma_prime, sigma;
    q = compute_q(sigma_prime);
    p_prime = compute_p_prime(sigma_prime);
    s = compute_s(sigma_prime, p_prime);
    dq_dsigma_prime = compute_dq_dsigma_prime(sigma_prime, s, q);
    dJ_3_dsigma_prime = compute_dJ_3_dsigma_prime(sigma_prime, s, q);
    sigma = compute_sigma(sigma_prime, u);
    compute_stress_invariants(sigma, I_1, I_2, I_3, J_1, J_2, J_3);
    compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
    
    // Compute derivatives.
    double df_dq = SMCC_DF_DQ;
    double df_dp_prime = SMCC_DF_DP_PRIME;
    double df_dtheta = SMCC_DF_DTHETA;
    double pi = 2*std::acos(0.0);
    Cauchy one = Cauchy::Constant(1.0); 
    if (q > 0.0 && df_dtheta != 0.0) {
        df_dsigma_prime = (df_dp_prime*dp_dsigma_prime) + ((df_dq - df_dtheta*tan(3*theta_s)/q)*dq_dsigma_prime) + (one*(sqrt(3)/(2.0*pow(q,3)*cos(3*theta_s)))*df_dtheta);
    } else { 
        df_dsigma_prime = (df_dp_prime*dp_dsigma_prime) + (df_dq*dq_dsigma_prime);
    }
    
    double dg_dq = SMCC_DG_DQ;
    double dg_dp_prime = SMCC_DG_DP_PRIME;
    double dg_dtheta = SMCC_DG_DTHETA;
    if (q > 0.0 && dg_dtheta != 0.0) {
        dg_dsigma_prime = (dg_dp_prime*dp_dsigma_prime) + ((dg_dq - dg_dtheta*tan(3*theta_s)/q)*dq_dsigma_prime) + (one*(sqrt(3)/(2.0*pow(q,3)*cos(3*theta_s)))*dg_dtheta);
    } else { 
        dg_dsigma_prime = (dg_dp_prime*dp_dsigma_prime) + (dg_dq*dq_dsigma_prime);
    }
    a = to_voigt(df_dsigma_prime);
    b = to_voigt(dg_dsigma_prime);
    H = SMCC_HARDENING_MODULUS;
}

State SMCC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    double Delta_epsilon_vol_e = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_e));
    State elastic_state(state.size());
    elastic_state[0] = SMCC_STATE_0_ELASTIC_UPDATE;
    elastic_state[1] = SMCC_STATE_1_ELASTIC_UPDATE;
    elastic_state[2] = SMCC_STATE_2_ELASTIC_UPDATE;
    return elastic_state;
}

State SMCC::compute_plastic_state_variable_increment(Voigt Delta_epsilon_tilde_p, double delta_lambda, Cauchy df_dsigma_prime, double H) {
    double Delta_epsilon_vol_p = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_p));
    State delta_state(state.size());
    delta_state[0] = SMCC_STATE_0_PLASTIC_INCREMENT;
    delta_state[1] = SMCC_STATE_1_PLASTIC_INCREMENT;
    delta_state[2] = SMCC_STATE_2_PLASTIC_INCREMENT;
    return delta_state;
}

State SMCC::compute_plastic_state_variable_increment(double delta_lambda, Cauchy df_dsigma_prime, double H) {
    // Note: only correct state variables that do not depend on the magnitude of the strain increment (hence strain increment is not passed in).
    Voigt Delta_epsilon_tilde_p = Voigt::Zero();
    return compute_plastic_state_variable_increment(Delta_epsilon_tilde_p, delta_lambda, df_dsigma_prime, H);
}