#include "C2MC.hpp"

C2MC::C2MC(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state), Elastoplastic::Elastoplastic() {   
    PLOG_INFO << "In C2MC init...";
    set_model_name("C2MC");
    set_model_type("Elastoplastic");

    // Initialise logger.
    initialise_log(log_severity);

    // Check inputs.
    check_inputs(get_model_name(), (int)parameters.size(), (int)state.size(), parameters_required, state_required);
}

Constitutive C2MC::compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon) {
    // Compute elastic constants.
    double K = compute_K_given_E_nu(E, nu);
    double G = compute_G_given_E_nu(E, nu);

    // Compute elastic matrix.
    Constitutive D_e = compute_isotropic_linear_elastic_matrix(K, G);
    return D_e;
}

double C2MC::compute_f(Cauchy sigma_prime, State state) {
    // Stress invariants.
    double p_prime = compute_p_prime(sigma_prime);

    using namespace std;
    double f, I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar;
    double A, B, C, K_theta;
    {   
        compute_stress_invariants(sigma_prime, I_1, I_2, I_3, J_1, J_2, J_3);
        compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
        double sigma_bar = compute_sigma_bar(J_2);
        compute_A_B_C(phi_r, theta_s_bar, A, B, C);
        compute_K_theta(phi_r, theta_s_bar, A, B, C, K_theta);
        f = sqrt(pow(sigma_bar,2.0)*pow((K_theta),2.0) + pow(a_h,2.0)*pow(sin(phi_r),2.0)) - c_prime*cos(phi_r) - p_prime*sin(phi_r); 
    }
    return f;
}

Derivatives C2MC::compute_derivatives(Cauchy sigma_prime, State state) {
    using namespace std;
    // Compute mean effective stress, deviatoric stress tensor and derivatives of the stress state for current stress state.
    double p_prime, I_1, I_2, I_3, J_1, J_2, J_3, theta_c, theta_s, theta_s_bar, sigma_bar;
    Cauchy df_dsigma_prime, dg_dsigma_prime, s, dq_dsigma_prime, dJ_3_dsigma_prime, sigma, dsigma_bar_dsigma_prime;
    p_prime = compute_p_prime(sigma_prime);
    s = compute_s(sigma_prime, p_prime);
    compute_stress_invariants(sigma_prime, I_1, I_2, I_3, J_1, J_2, J_3);
    compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
    sigma_bar = compute_sigma_bar(J_2);
    dsigma_bar_dsigma_prime = compute_dsigma_bar_dsigma_prime(sigma_prime, s, sigma_bar);
    dJ_3_dsigma_prime = compute_dJ_3_dsigma_prime(sigma_prime, s, sigma_bar);

    // Derivative of yield surface.
    double C_1, C_2, C_3;
    compute_derivative_coefficients(sigma_bar, phi_r, theta_s_bar, C_1, C_2, C_3);
    df_dsigma_prime = (C_1*dp_prime_dsigma_prime) + (C_2*dsigma_bar_dsigma_prime) + (pow(sigma_bar,2.0)*C_3)*((1.0/(pow(sigma_bar,2.0)))*dJ_3_dsigma_prime);

    // If non-associated flow, compute the gradient of the plastic potential function.
    if (phi_r != psi_r) {
        compute_derivative_coefficients(sigma_bar, psi_r, theta_s_bar, C_1, C_2, C_3);
        dg_dsigma_prime = (C_1*dp_prime_dsigma_prime) + (C_2*dsigma_bar_dsigma_prime) + (pow(sigma_bar,2.0)*C_3)*((1.0/(pow(sigma_bar,2.0)))*dJ_3_dsigma_prime);
    } else {
        dg_dsigma_prime = df_dsigma_prime;
    }

    // Return Derivatives object.
    Derivatives derivatives;
    derivatives.df_dsigma_prime = df_dsigma_prime;
    derivatives.dg_dsigma_prime = dg_dsigma_prime;
    return derivatives;
}

void C2MC::compute_K_theta(double angle_r, double theta_s_bar, double A, double B, double C, double &K_theta) {
    using namespace std;
    if (abs(theta_s_bar) > theta_t_r) {
        // C2 rounding.
        K_theta = A + B*sin(3.0*theta_s_bar) + C*(pow(sin(3.0*theta_s_bar),2.0));
        PLOG_DEBUG << "Coefficients: A = " << A << "; B = " << B << "; C = " << C << "; K_theta = " << K_theta;
    } else { 
        // No rounding required.
        K_theta = cos(theta_s_bar)-(1.0/sqrt(3.0))*sin(angle_r)*sin(theta_s_bar);
        PLOG_DEBUG << "Coefficients: K_theta = " << K_theta;
    } 
}

void C2MC::compute_A_B_C(double angle_r, double theta_s_bar, double &A, double &B, double &C) {
    using namespace std;
    double A_1, A_2, B_1, B_2, C_1, C_2, sign_theta; 

    // Sign of Lode's angle.
    sign_theta = (theta_s_bar < 0.0) ? -1.0 : 1.0;

    // Coefficient C.
    C_1 = (-cos(3.0*theta_t_r)*cos(theta_t_r) - 3.0*sin(3.0*theta_t_r)*sin(theta_t_r))/(18.0*pow(cos(3*theta_t_r),3.0));
    C_2 = (1.0/sqrt(3.0))*((cos(3*theta_t_r)*sin(theta_t_r) - 3.0*sin(3.0*theta_t_r)*cos(theta_t_r))/(18.0*pow(cos(3.0*theta_t_r),3.0)));
    C = C_1 + C_2*sign_theta*sin(angle_r);

    // Coefficient B.
    B_1 = (cos(theta_t_r)*sin(6.0*theta_t_r) - 6.0*cos(6.0*theta_t_r)*sin(theta_t_r))/(18.0*pow(cos(3.0*theta_t_r),3.0));
    B_2 = -(sin(theta_t_r)*sin(6.0*theta_t_r) + 6.0*cos(6.0*theta_t_r)*cos(theta_t_r))/(18.0*sqrt(3.0)*pow(cos(3.0*theta_t_r),3.0));
    B = B_1*sign_theta + B_2*sin(angle_r);

    // Coefficient A.
    A_1 = cos(theta_t_r) - B_1*sin(3.0*theta_t_r) - C_1*pow(sin(3.0*theta_t_r),2.0);
    A_2 = (-1.0/sqrt(3.0))*sin(theta_t_r) - B_2*sin(3*theta_t_r) - C_2*pow(sin(3.0*theta_t_r),2.0);
    A = A_1 + A_2*sign_theta*sin(angle_r);
}

void C2MC::compute_derivative_coefficients(double sigma_bar, double angle_r, double theta_s_bar, double &C_1, double &C_2, double &C_3) {
    using namespace std;
    
    // Coefficients.
    double alpha, K_theta, A, B, C, dK_dtheta;
    compute_A_B_C(angle_r, theta_s_bar, A, B, C);
    compute_K_theta(angle_r, theta_s_bar, A, B, C, K_theta);

    // Derivative coefficients.
    alpha = (sigma_bar*K_theta)/sqrt(pow(sigma_bar,2.0)*pow(K_theta,2.0) + pow(a_h,2.0)*pow(sin(psi_r),2.0));
    C_1 = -sin(psi_r);
    if (abs(theta_s_bar) > theta_t_r) {
        // C2 rounding.
        dK_dtheta = 3.0*B*cos(3.0*theta_s_bar) + 3.0*C*sin(6.0*theta_s_bar);
        C_2 = alpha*(A - 2.0*B*sin(3.0*theta_s_bar) - 5.0*C*pow(sin(3.0*theta_s_bar),2.0));
        C_3 = alpha*(-(3.0*sqrt(3.0)/(2.0*pow(sigma_bar,2.0)))*(B + 2.0*C*sin(3.0*theta_s_bar)));
    } else {
        // No rounding required.
        dK_dtheta = -sin(theta_s_bar)-(1.0/sqrt(3.0))*sin(psi_r)*cos(theta_s_bar);
        C_2 = alpha*(K_theta - dK_dtheta*(tan(3.0*theta_s_bar)));
        C_3 = alpha*(-(sqrt(3.0)/(2.0*pow(sigma_bar,2)*cos(3.0*theta_s_bar)))*dK_dtheta);
    }
}