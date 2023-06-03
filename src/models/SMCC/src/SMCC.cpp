#include "SMCC.hpp"

SMCC::SMCC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state) {   
    set_model_name("SMCC");
    set_model_type("Elastoplastic");

    // Initialise logger.
    Logging::initialise_log(log_severity);

    // Check inputs.
    Checks::check_inputs(name, parameters.size(), state.size(), parameters_required, state_required);
}

State SMCC::get_state_variables(void) {
    return state;
}

void SMCC::set_state_variables(State new_state) {
    state = new_state;
}

Constitutive SMCC::compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon) {

    /* USER DEFINED CODE STARTS HERE */
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
    Constitutive D_e =  compute_isotropic_linear_elastic_matrix(K, G);
    /* USER DEFINED CODE ENDS HERE */

    return D_e;
}

double SMCC::compute_f(Cauchy sigma_prime, State state) {

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);

    // State variables.
    double e = state[0];
    double p_c = state[1];
    double s_ep = state[2];

    /* USER DEFINED CODE STARTS HERE */
    double f = std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c*s_ep);
    /* USER DEFINED CODE ENDS HERE */

    return f;
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
    
    /* USER DEFINED CODE STARTS HERE */
    double df_dq = 2*q;
    double df_dp_prime = std::pow(M,2)*(2*p_prime-p_c*s_ep);
    double df_dtheta = 0;
    /* USER DEFINED CODE ENDS HERE */

    double pi = 2*std::acos(0.0);
    Cauchy one = Cauchy::Constant(1.0); 
    if (q > 0.0 && df_dtheta != 0.0) {
        df_dsigma_prime = (df_dp_prime*dp_dsigma_prime) + ((df_dq - df_dtheta*tan(3*theta_s)/q)*dq_dsigma_prime) 
            + (one*(sqrt(3)/(2.0*pow(q,3)*cos(3*theta_s)))*df_dtheta);
    } else { 
        df_dsigma_prime = (df_dp_prime*dp_dsigma_prime) + (df_dq*dq_dsigma_prime);
    }
    
    /* USER DEFINED CODE STARTS HERE */
    double dg_dq = 2*q;
    double dg_dp_prime = std::pow(M,2)*(2*p_prime-p_c*s_ep);
    double dg_dtheta = 0;
    /* USER DEFINED CODE ENDS HERE */

    if (q > 0.0 && dg_dtheta != 0.0) {
        dg_dsigma_prime = (dg_dp_prime*dp_dsigma_prime) + ((dg_dq - dg_dtheta*tan(3*theta_s)/q)*dq_dsigma_prime)
            + (one*(sqrt(3)/(2.0*pow(q,3)*cos(3*theta_s)))*dg_dtheta);
    } else { 
        dg_dsigma_prime = (dg_dp_prime*dp_dsigma_prime) + (dg_dq*dq_dsigma_prime);
    }
    a = to_voigt(df_dsigma_prime);
    b = to_voigt(dg_dsigma_prime);

    /* USER DEFINED CODE STARTS HERE */
    double H_p_c = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*s_ep*df_dsigma_prime.trace();
    double H_s_ep = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*-k*(s_ep-1.0)*std::sqrt((1-A)*std::pow(df_dsigma_prime.trace(),2)
        + (A*2.0/3.0*(df_dsigma_prime*(df_dsigma_prime.transpose())).trace()));
    H = H_p_c + H_s_ep;
    /* USER DEFINED CODE ENDS HERE */

}

State SMCC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    double Delta_epsilon_vol_e = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_e));
    State elastic_state(state.size());

    /* USER DEFINED CODE STARTS HERE */
    elastic_state[0] = e-(e*Delta_epsilon_vol_e);
    elastic_state[1] = p_c;
    elastic_state[2] = s_ep;
    /* USER DEFINED CODE ENDS HERE */

    return elastic_state;
}

State SMCC::compute_plastic_state_variable_increment(double delta_lambda, Cauchy df_dsigma_prime, double H, Voigt Delta_epsilon_tilde_p) {
    double Delta_epsilon_vol_p = compute_Delta_epsilon_vol(to_cauchy(Delta_epsilon_tilde_p));
    State delta_state(state.size());

    /* USER DEFINED CODE STARTS HERE */
    delta_state[0] = -(1+e)*Delta_epsilon_vol_p;
    delta_state[1] = delta_lambda*((((std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star))*s_ep*df_dsigma_prime.trace())/(std::pow(M,2)*p_prime*s_ep));
    delta_state[2] = delta_lambda*((std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*-k*(s_ep-1.0)*std::sqrt((1-A)*std::pow(df_dsigma_prime.trace(),2)
        + (A*2.0/3.0*(df_dsigma_prime*(df_dsigma_prime.transpose())).trace())))/(std::pow(M,2)*p_prime*p_c);
    /* USER DEFINED CODE ENDS HERE */

    return delta_state;
}