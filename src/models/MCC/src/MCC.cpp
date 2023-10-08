#include "MCC.hpp"

MCC::MCC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state), Elastoplastic::Elastoplastic() {   
    set_model_name("MCC");
    set_model_type("Elastoplastic");

    // Initialise logger.
    initialise_log(log_severity);

    // Check inputs.
    check_inputs(get_model_name(), (int)parameters.size(), (int)state.size(), parameters_required, state_required);
}

Constitutive MCC::compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon) {
    double Delta_epsilon_e_vol = compute_Delta_epsilon_vol(Delta_epsilon);
    double p_prime = compute_p_prime(sigma_prime);
    double K = compute_K_Butterfield(p_prime, Delta_epsilon_e_vol, kappa_star, settings.EPS);
    double G = compute_G_given_K_nu(K, nu);
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K, G);  
    return D_e;
}

double MCC::compute_f(Cauchy sigma_prime, State state) {
    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);
    
    // State variables.
    double p_c = state(0);
    
    // Yield surface function.
    using namespace std;
    double f = pow(q,2) + pow(M,2)*p_prime*(p_prime-p_c);

    return f;
}

Derivatives MCC::compute_derivatives(Cauchy sigma_prime, State state) {
    // State variables.
    double p_c = state(0);

    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double q, p_prime;
    Cauchy df_dsigma_prime, dg_dsigma_prime, s, dq_dsigma_prime;
    HardeningModuli H_s(state.size());
    StateFactors B_s(state.size());
    q = compute_q(sigma_prime);
    p_prime = compute_p_prime(sigma_prime);
    s = compute_s(sigma_prime, p_prime);
    dq_dsigma_prime = compute_dq_dsigma_prime(sigma_prime, s, q);
    
    // Yield surface derivatives.
    using namespace std;
    double df_dq = 2*q;
    double df_dp_prime = pow(M,2)*(2*p_prime-p_c);

    // Derivatives of yield surface and plastic potential function.
    df_dsigma_prime = (df_dp_prime*dp_prime_dsigma_prime) + (df_dq*dq_dsigma_prime);
    dg_dsigma_prime = df_dsigma_prime; // Associated flow.

    // Hardening moduli.
    H_s(0) = (pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*tr(dg_dsigma_prime);

    // State variable increment factors.
    double dg_dp_c = -pow(M,2)*p_prime;
    B_s(0) = H_s(0)/dg_dp_c;

    // Return Derivatives object.
    Derivatives derivatives;
    derivatives.df_dsigma_prime = df_dsigma_prime;
    derivatives.dg_dsigma_prime = dg_dsigma_prime;
    derivatives.H_s = H_s;
    derivatives.B_s = B_s;
    return derivatives;
}