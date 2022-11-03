#include "MCC.hpp"
#include "MCC_Definition.hpp"

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
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
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

void MCC::compute_derivatives(Cauchy sigma_prime, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &dg_dp_prime, double &H) {
    df_dsigma_prime = DF_DSIGMA_PRIME;
    a = df_dsigma_prime.voigt();
    dg_dsigma_prime = DG_DSIGMA_PRIME;
    b = dg_dsigma_prime.voigt();
    dg_dp_prime = DG_DP_PRIME;
    H = HARDENING_MODULUS;    
}