#include "Elastic.hpp"

Elastic::Elastic(std::string log_severity) {
    // Initialise log.
    initialise_log(log_severity);
}

double Elastic::compute_K_given_E_nu(double E, double nu) {
    double K = E/(3.0*(1.0-2.0*nu));
    PLOG_VERBOSE << "Bulk modulus, K = " << K;
    return K;
}

double Elastic::compute_K_Butterfield(double p_prime, double Delta_epsilon_e_vol, double kappa_star, double tolerance) {
    Delta_epsilon_e_vol = abs(Delta_epsilon_e_vol) > tolerance ? Delta_epsilon_e_vol : 0.0;
    double K;
    if (Delta_epsilon_e_vol != 0.0) {
        // Secant bulk modulus.
        K = (p_prime/Delta_epsilon_e_vol)*(exp(Delta_epsilon_e_vol/kappa_star)-1.0);
    } else {
        // Tangent bulk modulus.
        K = p_prime/kappa_star;
    }
    PLOG_VERBOSE << "Bulk modulus, K = " << K;
    return K;
}

double Elastic::compute_G_given_E_nu(double E, double nu) {
    double G = E/(2.0*(1.0+nu));
    PLOG_VERBOSE << "Shear modulus, G = " << G; 
    return G;
}

double Elastic::compute_G_given_K_nu(double K, double nu) {
    double G = (3.0*(1.0-2.0*nu)*K)/(2.0*(1.0+nu));
    PLOG_VERBOSE << "Shear modulus, G = " << G; 
    return G;
}

Constitutive Elastic::compute_isotropic_linear_elastic_matrix(double K, double G) {
    // Check elastic paramaters are valid.
    check_elastic_parameters(K, G);

    // Fill elastic matrix with isotropic linear elastic coefficients.
    Constitutive D_e = Constitutive::Zero();
    D_e(0,0) = D_e(1,1) = D_e(2,2) += K + 4.0/3.0*G; 
    D_e(0,1) = D_e(0,2) = D_e(1,2) = D_e(1,0) = D_e(2,0) = D_e(2,1) += K - 2.0/3.0*G;
    D_e(3,3) = D_e(4,4) = D_e(5,5) += G; 
    return D_e;
}

Voigt Elastic::compute_elastic_stress_increment(Constitutive D_e, Voigt Delta_epsilon_tilde) {
    return D_e*Delta_epsilon_tilde;
}

Cauchy Elastic::compute_elastic_stress(Cauchy sigma_prime, Voigt Delta_epsilon_tilde) {
    Constitutive D_e = compute_D_e(sigma_prime, to_cauchy(Delta_epsilon_tilde));
    Voigt Delta_sigma_prime_tilde = compute_elastic_stress_increment(D_e, Delta_epsilon_tilde);
    return sigma_prime + to_cauchy(Delta_sigma_prime_tilde);
}

void Elastic::solve(void) {
    // Compute elastic matrix. 
    D_e = compute_D_e(sigma_prime, Delta_epsilon);

    // Update stress state.
    Voigt Delta_sigma_prime_tilde = D_e*Delta_epsilon_tilde;
    sigma_prime_tilde += Delta_sigma_prime_tilde;
    sigma_prime = to_cauchy(sigma_prime_tilde);
    p_prime = compute_p_prime(sigma_prime);
    q = compute_q(sigma_prime);

    // Take the Jacobian as the tangent stiffness.
    jacobian = D_e;
}

void Elastic::check_elastic_parameters(double K, double G) {
    PLOG_FATAL_IF(std::isnan(K)) << "Bulk modulus is " << K << ".";
    PLOG_FATAL_IF(std::isnan(G)) << "Shear modulus is " << G << ".";
    PLOG_FATAL_IF(K <= 0.0) << "Bulk modulus less than or equal to zero.";
    PLOG_FATAL_IF(G <= 0.0) << "Shear modulus less than or equal to zero.";
    if (K > 0.0 && G > 0.0) {
        return; 
    } else {
        assert(false);
        throw std::invalid_argument("Invalid elastic parameters.");
    }
}