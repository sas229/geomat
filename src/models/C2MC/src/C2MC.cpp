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
    // Compute elastic constants.
    double K = compute_K(E, nu);
    double G = compute_G(E, nu);

    // Compute elastic matrix.
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K, G);
    return D_e;
}

double C2MC::compute_f(Cauchy sigma_prime, State state) {

    // Stress invariants.
    double p_prime = compute_p_prime(sigma_prime); // C2MC definition is tension positive.

    // State variables.
    // No state variables for this model.

    using namespace std;
    double f, I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar;
    double A, B, C, K_theta;
    {   
        compute_stress_invariants(sigma_prime, I_1, I_2, I_3, J_1, J_2, J_3);
        compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
        double sigma_bar = compute_sigma_bar(J_2);
        compute_coefficients(phi_r, theta_s_bar, A, B, C, K_theta);
        f = -p_prime*sin(phi_r) + sqrt(pow(sigma_bar,2.0)*pow((K_theta),2.0) + pow(a_h,2.0)*pow(sin(phi_r),2.0)) - cohs*cos(phi_r); 
    }
    return f;
}

void C2MC::compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, HardeningModuli  &H_s) {
    // State variables.
    // No state variables for this model.

    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double p_prime, I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar, sigma_bar;
    Cauchy s, dq_dsigma_prime, dJ_3_dsigma_prime, sigma, dsigma_bar_dsigma_prime;
    p_prime = compute_p_prime(sigma_prime);
    s = compute_s(sigma_prime, p_prime);
    compute_stress_invariants(sigma_prime, I_1, I_2, I_3, J_1, J_2, J_3);
    compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
    sigma_bar = compute_sigma_bar(J_2);
    dsigma_bar_dsigma_prime = compute_dsigma_bar_dsigma_prime(sigma_prime, s, sigma_bar);
    dJ_3_dsigma_prime = compute_dJ_3_dsigma_prime(sigma_prime, s, sigma_bar);

    double df_dp_prime, df_dsigma_bar, df_dtheta, dg_dp_prime, dg_dsigma_bar, dg_dtheta, C_1, C_2, C_3;
    {
        using namespace std;
        double A, B, C, K_theta, dK_dtheta;
        compute_coefficients(phi_r, theta_s_bar, A, B, C, K_theta);
        double alpha = (sigma_bar*K_theta)/sqrt(pow(sigma_bar,2.0)*pow(K_theta,2.0) + pow(a_h,2.0)*pow(sin(phi_r),2.0));
        if (abs(theta_s_bar) > theta_tr) {
            dK_dtheta = 3.0*B*cos(3.0*theta_s_bar) + 3.0*C*sin(6.0*theta_s_bar);
        } else { 
            dK_dtheta = -sin(theta_s_bar)-(1.0/sqrt(3.0))*sin(phi_r)*cos(theta_s_bar);
        }
        C_1 = -sin(phi_r);
        if (abs(theta_s_bar) > theta_tr) {
            C_2 = alpha*(A - 2*B*sin(3.0*theta_s_bar) - 5.0*C*pow(sin(3.0*theta_s_bar),2.0));
            C_3 = alpha*(-(3.0*sqrt(3.0)/(2.0*pow(sigma_bar,2)))*(B + 2.0*C*sin(3.0*theta_s_bar)));
        } else {
            C_2 = alpha*(K_theta - dK_dtheta*(tan(3.0*theta_s_bar)));
            C_3 = alpha*(-(sqrt(3.0)/(2.0*pow(sigma_bar,2)*cos(3.0*theta_s_bar)))*dK_dtheta);
        }
    }
    df_dsigma_prime = (C_1*dp_prime_dsigma_prime) + (C_2*dsigma_bar_dsigma_prime) + (pow(sigma_bar,2.0)*C_3)*((1.0/(pow(sigma_bar,2.0)))*dJ_3_dsigma_prime);

    // If non-associated flow, compute the gradient of the plastic potential function.
    if (phi_r != psi_r) {
        {
            using namespace std;
            double A, B, C, K_theta, dK_dtheta;
            compute_coefficients(psi_r, theta_s_bar, A, B, C, K_theta);
            double alpha = (sigma_bar*K_theta)/sqrt(pow(sigma_bar,2.0)*pow(K_theta,2.0) + pow(a_h,2.0)*pow(sin(psi_r),2.0));
            if (abs(theta_s_bar) > theta_tr) {
                dK_dtheta = 3.0*B*cos(3.0*theta_s_bar) + 3.0*C*sin(6.0*theta_s_bar);
            } else { 
                dK_dtheta = -sin(theta_s_bar)-(1.0/sqrt(3.0))*sin(psi_r)*cos(theta_s_bar);
            }
            C_1 = -sin(psi_r);
            if (abs(theta_s_bar) > theta_tr) {
                C_2 = alpha*(A - 2.0*B*sin(3.0*theta_s_bar) - 5.0*C*pow(sin(3.0*theta_s_bar),2.0));
                C_3 = alpha*(-(3.0*sqrt(3.0)/(2.0*pow(sigma_bar,2.0)))*(B + 2.0*C*sin(3.0*theta_s_bar)));
            } else {
                C_2 = alpha*(K_theta - dK_dtheta*(tan(3.0*theta_s_bar)));
                C_3 = alpha*(-(sqrt(3.0)/(2.0*pow(sigma_bar,2)*cos(3.0*theta_s_bar)))*dK_dtheta);
            }
        }
        dg_dsigma_prime = (C_1*dp_prime_dsigma_prime) + (C_2*dsigma_bar_dsigma_prime) + (pow(sigma_bar,2.0)*C_3)*((1.0/(pow(sigma_bar,2.0)))*dJ_3_dsigma_prime);
    } else {
        dg_dsigma_prime = df_dsigma_prime;
    }

    // Derivatives in Voigt form.
    a = to_voigt(df_dsigma_prime);
    b = to_voigt(dg_dsigma_prime);

    // Hardening modulus.
    H = 0.0;
}

State C2MC::compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) {
    State elastic_state(state.size());
    return elastic_state;
}

State C2MC::compute_plastic_state_variable_increment(double delta_lambda, Cauchy df_dsigma_prime, HardeningModuli  H_s, Voigt Delta_epsilon_tilde_p) {
    State delta_state(state.size());
    return delta_state;
}

void C2MC::compute_coefficients(double angle_r, double theta_s_bar, double &A, double &B, double &C, double &K_theta) {
    using namespace std;
    double A_1, B_1, C_1, A_2, B_2, C_2, sign_theta;

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

        // Coefficient K_theta. 
        K_theta = A + B*sin(3.0*theta_s_bar) + C*(pow(sin(3.0*theta_s_bar),2.0));
        PLOG_DEBUG << "Coefficients: A = " << A << "; B = " << B << "; C = " << C << "; K_theta = " << K_theta;
    } else { 
        K_theta = cos(theta_s_bar)-(1.0/sqrt(3.0))*sin(angle_r)*sin(theta_s_bar);
        PLOG_DEBUG << "Coefficients: K_theta = " << K_theta;
    } 
}