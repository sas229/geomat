#include "SMCC.hpp"

SMCC::SMCC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state), Elastoplastic::Elastoplastic() {   
    set_model_name("SMCC");
    set_model_type("Elastoplastic");

    // Initialise logger.
    Logging::initialise_log(log_severity);

    // Check inputs.
    Checks::check_inputs(name, (int)parameters.size(), (int)state.size(), parameters_required, state_required);
}

State SMCC::get_state_variables(void) {
    return state;
}

void SMCC::set_state_variables(State new_state) {
    state = new_state;
}

Constitutive SMCC::compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon) {
    using namespace std;
    double Delta_epsilon_e_vol = compute_Delta_epsilon_vol(Delta_epsilon);
    double K = compute_K_Butterfield(p_prime, Delta_epsilon_e_vol, kappa_star, settings.EPS);
    double G = compute_G_given_K_nu(K, nu);
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K, G);  
    return D_e;
}

double SMCC::compute_f(Cauchy sigma_prime, State state) {
    using namespace std;
    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);

    // State variables.
    double p_c = state(0);
    double s_ep = state(1);

    // Yield surface function.
    double f = pow(q,2) + pow(M,2)*p_prime*(p_prime-p_c*s_ep);
    return f;
}

void SMCC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Cauchy &dg_dsigma_prime, HardeningModuli &H_s, StateFactors &B_s) {
    using namespace std;
    // State variables.
    double p_c = state(0);
    double s_ep = state(1);

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
    
    // Yield surface derivatives.
    double df_dq = 2*q;
    double df_dp_prime = pow(M,2)*(2*p_prime-p_c*s_ep);
    double df_dtheta = 0;

    // Derivatives of yield surface and plastic potential function.
    df_dsigma_prime = (df_dp_prime*dp_prime_dsigma_prime) + (df_dq*dq_dsigma_prime);
    dg_dsigma_prime = df_dsigma_prime; // Associated flow.

    // Hardening moduli.
    H_s(0) = (pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*s_ep*df_dsigma_prime.trace();
    H_s(1) = (pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*-k*(s_ep-1.0)*sqrt((1-A)*pow(df_dsigma_prime.trace(),2)
        + (A*2.0/3.0*(double_dot_product(dev(df_dsigma_prime)))));

    // State variable increment factors.
    double df_dp_c = -pow(M,2)*p_prime*s_ep;
    double df_ds_ep = -pow(M,2)*p_prime*p_c;
    B_s(0) = H_s(0)/df_dp_c;
    B_s(1) = H_s(1)/df_ds_ep;
}