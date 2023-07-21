#include "C2MC.hpp"

C2MC::C2MC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state), Elastoplastic::Elastoplastic() {   
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
    double K = compute_K(E, nu);
    double G = compute_G(E, nu);
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K, G);
    /* USER DEFINED CODE ENDS HERE */

    return D_e;
}

double C2MC::compute_f(Cauchy sigma_prime, State state) {

    // Stress invariants.
    double q = compute_q(sigma_prime);
    double p_prime = compute_p_prime(sigma_prime); // C2MC definition is tension positive.

    // State variables.
    // No state variables for this model.

    /* USER DEFINED CODE STARTS HERE */
    using namespace std;
    double f, I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar;
    double A, B, C, k_theta;
    {   
        compute_stress_invariants(sigma_prime, I_1, I_2, I_3, J_1, J_2, J_3);
        compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
        compute_coefficients(phi_r, theta_s_bar, A, B, C, k_theta);
        f = -p_prime*sin(phi_r) + sqrt(pow((q/sqrt(3.0)),2.0)*pow((k_theta),2.0) + pow(a_h,2.0)*pow(sin(phi_r),2.0)) - cohs*cos(phi_r); 
    }
    /* USER DEFINED CODE ENDS HERE */
    return f;
}

void C2MC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H) {
    // State variables.
    // No state variables for this model.

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
    double df_dp_prime, df_dq, df_dtheta;
    {
        using namespace std;
        double A, B, C, k_theta, dk_dtheta_f;
        compute_coefficients(phi_r, theta_s_bar, A, B, C, k_theta);
        double alpha_f = ((q/sqrt(3.0))*k_theta)/sqrt(pow(((q/sqrt(3.0))*k_theta),2.0) + pow((a_h*sin(phi_r)),2.0));
        if (abs(theta_s_bar) > theta_tr) {
            dk_dtheta_f = 3.0*B*cos(3.0*theta_s_bar) + 3.0*C*sin(6.0*theta_s_bar);
        } else { 
            dk_dtheta_f = -sin(theta_s_bar)-(1.0/sqrt(3.0))*sin(phi_r)*cos(theta_s_bar);
        }
        df_dp_prime = -sin(phi_r);
        df_dq = alpha_f*k_theta;
        df_dtheta = dk_dtheta_f;
    }
    /* USER DEFINED CODE ENDS HERE */

    Cauchy one = Cauchy::Constant(1.0); 
    if (q > 0.0 && df_dtheta != 0.0) {
        df_dsigma_prime = (df_dp_prime*dp_dsigma_prime) + ((sqrt(3.0)*df_dq - df_dtheta*tan(3.0*theta_s_bar))*dq_dsigma_prime) 
            - (one*(sqrt(3.0)/(2.0*pow((q/sqrt(3.0)),2.0)*cos(3.0*theta_s_bar)))*df_dtheta);
    } else { 
        df_dsigma_prime = (df_dp_prime*dp_dsigma_prime) + (df_dq*dq_dsigma_prime);
    }
    
    /* USER DEFINED CODE STARTS HERE */
    double dg_dp_prime, dg_dq, dg_dtheta;
    {   
        // If associated flow, there is no need to compute these derivatives again.
        if (phi != psi) {
            using namespace std;
            double A, B, C, k_theta, dk_dtheta_g;
            compute_coefficients(psi_r, theta_s_bar, A, B, C, k_theta);
            double alpha_g = ((q/sqrt(3.0))*k_theta)/sqrt(pow(((q/sqrt(3.0))*k_theta),2.0)+pow((a_h*sin(psi_r)),2.0));
            if (abs(theta_s_bar) > theta_tr) {
                dk_dtheta_g = 3.0*B*cos(3.0*theta_s_bar)+3.0*C*sin(6.0*theta_s_bar);
            } else { 
                dk_dtheta_g = -sin(theta_s_bar)-(1.0/sqrt(3.0))*sin(psi_r)*cos(theta_s_bar);
            }
            dg_dp_prime = -sin(psi_r);
            dg_dq = alpha_g*k_theta;
            dg_dtheta = dk_dtheta_g;
        } else {
            dg_dp_prime = df_dp_prime;
            dg_dq = df_dq;
            dg_dtheta = df_dtheta;
        }
    }
    /* USER DEFINED CODE ENDS HERE */

    if (q > 0.0 && dg_dtheta != 0.0) {
        dg_dsigma_prime = (dg_dp_prime*dp_dsigma_prime) + ((sqrt(3.0)*dg_dq - dg_dtheta*tan(3.0*theta_s_bar))*dq_dsigma_prime) 
            - (one*(sqrt(3.0)/(2.0*pow((q/sqrt(3.0)),2.0)*cos(3.0*theta_s_bar)))*dg_dtheta);
    } else { 
        dg_dsigma_prime = (dg_dp_prime*dp_dsigma_prime) + (dg_dq*dq_dsigma_prime);
    }
    a = to_voigt(df_dsigma_prime);
    b = to_voigt(dg_dsigma_prime);

    /* USER DEFINED CODE STARTS HERE */
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

void C2MC::compute_coefficients(double angle_r, double theta_s_bar, double &A, double &B, double &C, double &k_theta) {
    using namespace std;
    double A_1, B_1, C_1, A_2, B_2, C_2, theta_tr, theta_r, sign_theta;

    // Only calculate A, B and C if required.
    if (abs(theta_s_bar) > theta_tr) {
        sign_theta = (theta_s_bar < 0.0) ? -1.0 : 1.0;

        // Coefficient C.
        C_1 = (-cos(3.0*theta_tr)*cos(theta_tr) - 3.0*sin(3.0*theta_tr)*sin(theta_tr))/(18.0*pow(cos(3*theta_tr),3.0));
        C_2 = (1.0/sqrt(3.0))*((cos(3*theta_tr)*sin(theta_tr) - 3.0*sin(3.0*theta_tr)*cos(theta_tr))/(18.0*pow(cos(3.0*theta_tr),3.0)));
        C = C_1 + C_2*sign_theta*sin(angle_r);

        // Coefficient B.
        B_1 = (cos(theta_tr)*sin(6.0*theta_tr) - 6.0*cos(6.0*theta_tr)*sin(theta_tr))/(18.0*pow(cos(3.0*theta_tr),3.0));
        B_2 = -(sin(theta_tr)*sin(6.0*theta_tr) + 6.0*cos(6.0*theta_tr)*cos(theta_tr))/(18.0*sqrt(3.0)*pow(cos(3.0*theta_tr),3.0));
        B = B_1*sign_theta + B_2*sin(angle_r);

        // Coefficient A.
        A_1 = cos(theta_tr) - B_1*sin(3.0*theta_tr) - C_1*pow(sin(3.0*theta_tr),2.0);
        A_2 = (-1.0/sqrt(3.0))*sin(theta_tr) - B_2*sin(3*theta_tr) - C_2*pow(sin(3.0*theta_tr),2.0);
        A = A_1 + A_2*sign_theta*sin(angle_r);

        // Coefficient k_theta. 
        k_theta = A + B*sin(3.0*theta_s_bar) + C*(pow(sin(3.0*theta_s_bar),2.0));
        PLOG_DEBUG << "A = " << A << "; B = " << B << "; C = " << C << "; k_theta = " << k_theta;
    } else { 
        k_theta = cos(theta_s_bar)-(1.0/sqrt(3.0))*sin(angle_r)*sin(theta_s_bar);
        PLOG_DEBUG << "k_theta = " << k_theta;
    } 
}