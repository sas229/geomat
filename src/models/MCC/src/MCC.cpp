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
    using namespace std;
    double Delta_epsilon_e_vol = compute_Delta_epsilon_vol(Delta_epsilon);
    double K = compute_K_Butterfield(p_prime, Delta_epsilon_e_vol, kappa_star, settings.EPS);
    double G = compute_G_given_K_nu(K, nu);
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K, G);  
    return D_e;
}

double MCC::compute_f(Cauchy sigma_prime, State state) {
    using namespace std;
    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    
    // State variables.
    double p_c = state(0);
    
    // Yield surface function.
    double f = pow(q,2) + pow(M,2)*p_prime*(p_prime-p_c);

    return f;
}

void MCC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Cauchy &dg_dsigma_prime, HardeningModuli &H_s, StateFactors &B_s) {
    using namespace std;
    // State variables.
    double p_c = state(0);

    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double q, p_prime;
    Cauchy s, dq_dsigma_prime;
    q = compute_q(sigma_prime);
    p_prime = compute_p_prime(sigma_prime);
    s = compute_s(sigma_prime, p_prime);
    dq_dsigma_prime = compute_dq_dsigma_prime(sigma_prime, s, q);
    
    // Yield surface derivatives.
    double df_dq = 2*q;
    double df_dp_prime = pow(M,2)*(2*p_prime-p_c);
    double df_dtheta = 0;

    // Derivatives of yield surface and plastic potential function.
    df_dsigma_prime = (df_dp_prime*dp_prime_dsigma_prime) + (df_dq*dq_dsigma_prime);
    dg_dsigma_prime = df_dsigma_prime; // Associated flow.

    // Hardening moduli.
    H_s(0) = (pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*df_dsigma_prime.trace();

    // State variable increment factors.
    double df_dp_c = -pow(M,2)*p_prime;
    B_s(0) = H_s(0)/df_dp_c;
}