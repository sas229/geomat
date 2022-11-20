#include "MCC.hpp"
#include "MCC_Definition.hpp"

MCC::MCC(State parameters, State state) : parameters(parameters), state(state) {
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
    
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    return YIELD;
}

double MCC::compute_K(double delta_epsilon_e_vol, double p_prime) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    if (delta_epsilon_e_vol != 0.0) {
        return BULK_MODULUS_SECANT;
    } else {
        return BULK_MODULUS_TANGENT;
    }
}

double MCC::compute_G(double K) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    return SHEAR_MODULUS;
}

void MCC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &dg_dp_prime, double &H) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */

    // State variables.
    double e = state[0];
    double p_c = state[1];

    // Compute mean effective stress and deviatoric stress tensor for current stress state.
    double p_prime = compute_p_prime(sigma_prime);
    Cauchy s = compute_s(sigma_prime, p_prime);

    // Compute derivatives.
    df_dsigma_prime = DF_DSIGMA_PRIME;
    a = df_dsigma_prime.voigt();
    dg_dsigma_prime = DG_DSIGMA_PRIME;
    b = dg_dsigma_prime.voigt();
    dg_dp_prime = DG_DP_PRIME;
    H = HARDENING_MODULUS;    
}

State MCC::compute_elastic_state_variable(Voigt delta_epsilon_tilde_e) {
    double delta_epsilon_vol_e = compute_delta_epsilon_vol(delta_epsilon_tilde_e.cauchy());
    State elastic_state(state.size());
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    elastic_state[0] = STATE_0_ELASTIC_UPDATE;
    elastic_state[1] = STATE_1_ELASTIC_UPDATE;
    return elastic_state;
}

State MCC::get_state_variables(void) {
    return state;
}

State MCC::compute_plastic_state_variable_increment(double delta_lambda, double H) {
    State delta_state(state.size());
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    delta_state[0] = STATE_0_PLASTIC_INCREMENT;
    delta_state[1] = STATE_1_PLASTIC_INCREMENT;
    return delta_state;
}

State MCC::compute_plastic_state_variable_correction(double delta_lambda, double H) {
    // Note: only correct state variables that do not depend on the magnitude of the strain increment.
    State delta_state_c(state.size());
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    delta_state_c[0] = 0;
    delta_state_c[1] = STATE_1_PLASTIC_INCREMENT;
    return delta_state_c;
}

void MCC::compute_plastic_state_variable(void) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
}
