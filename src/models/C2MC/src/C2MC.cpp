#include "C2MC.hpp"

C2MC::C2MC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state) {   
    set_model_name("C2MC");
    set_model_type("Elastoplastic");

    // Initialise logger.
    Logging::initialise_log(log_severity);

    // Check inputs.
    Checks::check_inputs(name, (int)parameters.size(), (int)state.size(), parameters_required, state_required);
}

State C2MC::get_state_variables(void) {
    return state;
}

void C2MC::set_state_variables(State new_state) {
    state = new_state;
}

Constitutive C2MC::compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon) {

    /* USER DEFINED CODE STARTS HERE */
    double K = E/(3.0*(1-2.0*nu));
    double G = E/(2.0*(1.0+nu));
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K, G);
    /* USER DEFINED CODE ENDS HERE */

    return D_e;
}

double C2MC::compute_f(Cauchy sigma_prime, State state) {

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime);

    // State variables.
    //double psi = state[0]; if we consider Bolton's dilation related formulation


    /* USER DEFINED CODE STARTS HERE */
    double J_2, J_3, sign_theta, k_theta;
    double theta_r = 1.0/3.0*(std::asin(-3.0*std::sqrt(3.0)/2.0*J_3/(std::pow(J_2,3.0/2.0))));
    if (theta_r >= 0) {
        sign_theta = 1.0;
    } else { 
        sign_theta = -1.0;
    }
    double theta_tr = to_radians(theta_t);
    double phi_r = to_radians(phi);
    double psi_r = to_radians(psi);

    double coeff_B = (-std::cos(3.0*theta_tr)*(std::cos(theta_tr)-(1.0/std::sqrt(3.0))*std::sin(phi_r)*(sign_theta)*std::sin(theta_tr))
        -3.0*(sign_theta)*std::sin(3.0*theta_tr)*((sign_theta)*std::sin(theta_tr)+1.0/std::sqrt(3.0)*std::sin(phi_r)*std::cos(theta_tr)))/
        (18.0*(std::pow(std::cos(3.0*theta_tr),3.0)));
    double coeff_C = ((sign_theta)*std::sin(6.0*theta_tr)*(std::cos(theta_tr)-(1.0/std::sqrt(3.0))*std::sin(phi_r)*(sign_theta)*std::sin(theta_tr))
        -6.0*std::cos(6.0*theta_tr)*((sign_theta)*std::sin(theta_tr)+1.0/std::sqrt(3.0)*std::sin(phi_r)*std::cos(theta_tr)))/
        (18.0*(std::pow(std::cos(3.0*theta_tr),3.0)));
    double coeff_A = -(1.0/std::sqrt(3.0))*std::sin(phi_r)*(sign_theta)*std::sin(theta_tr)-coeff_B*(sign_theta)*std::sin(3.0*theta_tr)
        -coeff_C*std::pow(std::sin(3.0*theta_tr),2.0)+std::cos(theta_tr);
    
    if (std::abs(theta_r) > theta_tr) {
        k_theta = coeff_A+coeff_B*std::sin(3.0*theta_r)+coeff_B*(std::pow(std::sin(3.0*theta_r),2));
    } else { 
        k_theta = std::cos(theta_r)-(1.0/std::sqrt(3.0))*std::sin(phi_r)*std::sin(theta_r);
    }   
     
    double f = p_prime*std::sin(phi_r)+std::sqrt(std::pow((q*k_theta),2.0)+std::pow(a_h*std::sin(phi_r),2.0))-cohs*std::cos(phi_r);
    /* USER DEFINED CODE ENDS HERE */

    return f;
}

void C2MC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H) {
    // State variables.
    //double e = state[0];
    //double p_c = state[1];
    //double s_ep = state[2];

    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double q, p_prime, I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar;
    double k_theta, coeff_A, coeff_B, coeff_C, theta_r, theta_tr, phi_r, dk_dtheta_f, dk_dtheta_g;
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
    double var_alpha_f = (q*k_theta)/sqrt(std::pow((q*k_theta),2.0)+std::pow((a_h*sin(phi_r)),2.0));

    if (std::abs(theta_r) > theta_tr) {
        dk_dtheta_f = 3.0*coeff_B*std::cos(3.0*theta_r)+3.0*coeff_C*std::sin(6.0*theta_r);
    } else { 
        dk_dtheta_f = -std::sin(theta_r)-(1.0/std::sqrt(3.0))*std::sin(phi_r)*std::cos(theta_r);
    }

    double df_dp_prime = std::sin(phi_r);
    double df_dq = var_alpha_f*k_theta;
    double df_dtheta = dk_dtheta_f;
    /* USER DEFINED CODE ENDS HERE */
    double psi_r = to_radians(psi);
    Cauchy one = Cauchy::Constant(1.0); 
    if (q > 0.0 && df_dtheta != 0.0) {
        df_dsigma_prime = (df_dp_prime*dp_dsigma_prime) + ((df_dq - df_dtheta*tan(3.0*theta_s))*dq_dsigma_prime) 
            + (one*(sqrt(3.0)/(2.0*pow(q,2.0)*cos(3.0*theta_s)))*df_dtheta);
    } else { 
        df_dsigma_prime = (df_dp_prime*dp_dsigma_prime) + (df_dq*dq_dsigma_prime);
    }
    
    /* USER DEFINED CODE STARTS HERE */
    double var_alpha_g = (q*k_theta)/sqrt(std::pow((q*k_theta),2.0)+std::pow((a_h*sin(psi_r)),2.0));

    if (std::abs(theta_r) > theta_tr) {
        dk_dtheta_g = 3.0*coeff_B*std::cos(3.0*theta_r)+3.0*coeff_C*std::sin(6.0*theta_r);
    } else { 
        dk_dtheta_g = -std::sin(theta_r)-(1.0/std::sqrt(3.0))*std::sin(psi_r)*std::cos(theta_r);
    }
    double dg_dp_prime = std::sin(psi_r);
    double dg_dq = var_alpha_f*k_theta;
    double dg_dtheta = dk_dtheta_g;
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
    // double H_p_c = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*s_ep*df_dsigma_prime.trace();
    // double H_s_ep = (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*-k*(s_ep-1.0)*std::sqrt((1-A)*std::pow(df_dsigma_prime.trace(),2)
    //     + (A*2.0/3.0*(df_dsigma_prime*(df_dsigma_prime.transpose())).trace()));
    H = 0.0;
    /* USER DEFINED CODE ENDS HERE */

}

State C2MC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    State elastic_state(state.size());
    return elastic_state;
}

State C2MC::compute_plastic_state_variable_increment(double delta_lambda, Cauchy df_dsigma_prime, double H, Voigt Delta_epsilon_tilde_p) {
    State delta_state(state.size());
    return delta_state;
}