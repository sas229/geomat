#include "MCC.hpp"
#include "MCC_Definition.hpp"

MCC::MCC(Parameters parameters, State state) : parameters(parameters), state(state) {
    set_name("MCC");
    int parameters_required = 5;
    int state_required = 2;
    PLOG_FATAL_IF(parameters.size() != parameters_required) << parameters.size() << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state.size() != state_required) << state.size() << " parameters supplied when " << state_required << " expected.";
    assert(parameters.size() == parameters_required && state.size() == state_required);
    PLOG_INFO << name << " model instantiated with " << parameters.size() << " parameters and " << state.size() << " state variables.";  
}

double MCC::compute_f(Cauchy sigma_prime, State state) {
    // State variables.
    double e = state[0];
    double p_c = state[1];

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    
    return MCC_YIELD;
}

double MCC::compute_f(double p_prime, double q, State state) {
    // State variables.
    double e = state[0];
    double p_c = state[1];
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    return MCC_YIELD;
}


double MCC::compute_K(double Delta_epsilon_e_vol, double p_prime) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    if (Delta_epsilon_e_vol != 0.0) {
        return MCC_SECANT_BULK_MODULUS;
    } else {
        return MCC_TANGENT_BULK_MODULUS;
    }
}

double MCC::compute_G(double K) {
    return MCC_SHEAR_MODULUS;
}

void MCC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &dg_dp_prime, double &H) {
    // State variables.
    double e = state[0];
    double p_c = state[1];

    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double p_prime = compute_p_prime(sigma_prime);
    dq_dsigma_prime = compute_dq_dsigma_prime(sigma_prime);
    
    // Compute derivatives.
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    df_dq = MCC_DF_DQ;
    df_dp_prime = MCC_DF_DP_PRIME;
    df_dsigma_prime = (df_dq*dq_dsigma_prime) + (df_dp_prime*dp_dsigma_prime);
    dg_dq = MCC_DG_DQ;
    dg_dp_prime = MCC_DG_DP_PRIME;
    dg_dsigma_prime = (dg_dq*dq_dsigma_prime) + (dg_dp_prime*dp_dsigma_prime);
    a = to_voigt(df_dsigma_prime);
    b = to_voigt(dg_dsigma_prime);
    H = MCC_HARDENING_MODULUS;
}

State MCC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    double Delta_epsilon_vol_e = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_e));
    State elastic_state(state.size());
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    elastic_state[0] = MCC_STATE_0_ELASTIC_UPDATE;
    elastic_state[1] = MCC_STATE_1_ELASTIC_UPDATE;
    return elastic_state;
}

State MCC::get_state_variables(void) {
    return state;
}

void MCC::set_state_variables(State new_state) {
    state = new_state;
}

State MCC::compute_plastic_state_variable_increment(Voigt Delta_epsilon_tilde_p, double delta_lambda, double H) {
    double Delta_epsilon_vol_p = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_p));
    State delta_state(state.size());
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    delta_state[0] = MCC_STATE_0_PLASTIC_INCREMENT;
    delta_state[1] = MCC_STATE_1_PLASTIC_INCREMENT;
    return delta_state;
}

State MCC::compute_plastic_state_variable_increment(double delta_lambda, double H) {
    // Note: only correct state variables that do not depend on the magnitude of the strain increment (hence strain increment is not passed in).
    Voigt Delta_epsilon_tilde_p = Voigt::Zero();
    return compute_plastic_state_variable_increment(Delta_epsilon_tilde_p, delta_lambda, H);
}
