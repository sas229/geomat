#include "Elastic.hpp"

Constitutive Elastic::compute_isotropic_linear_elastic_matrix(double K, double G) {
    // Check elastic paramaters are valid.
    PLOG_ERROR_IF(K <= 0.0) << "Bulk modulus less than or equal to zero.";
    PLOG_ERROR_IF(G <= 0.0) << "Shear modulus less than or equal to zero.";
    assert(K > 0.0 && G > 0.0);

    // Fill elastic matrix with isotropic linear elastic coefficients.
    Constitutive D_e = Constitutive::Zero();
    D_e(0,0) = D_e(1,1) = D_e(2,2) += K + 4.0/3.0*G; 
    D_e(0,1) = D_e(0,2) = D_e(1,2) = D_e(1,0) = D_e(2,0) = D_e(2,1) += K - 2.0/3.0*G;
    D_e(3,3) = D_e(4,4) = D_e(5,5) += G; 
    return D_e;
}

Cauchy Elastic::compute_isotropic_linear_elastic_trial_stress(Cauchy sigma_prime, double alpha, Voigt delta_epsilon_tilde) {
    Voigt delta_epsilon_tilde_trial = alpha*delta_epsilon_tilde;
    double delta_epsilon_e_vol = delta_epsilon_tilde_trial.cauchy().trace();
    double p_prime_trial = 1.0/3.0*sigma_prime.trace();
    double K_trial = compute_K(delta_epsilon_e_vol, p_prime_trial);
    double G_trial = compute_G(K_trial);
    Constitutive D_e_trial = compute_isotropic_linear_elastic_matrix(K_trial, G_trial);
    Voigt delta_sigma_prime_tilde_trial = compute_elastic_stress_increment(D_e_trial, delta_epsilon_tilde_trial);
    return sigma_prime + delta_sigma_prime_tilde_trial.cauchy();
}

Voigt Elastic::compute_elastic_stress_increment(Constitutive D_e, Voigt delta_epsilon_tilde) {
    return D_e*delta_epsilon_tilde;
}

void Elastic::solve(void) {
    double K = compute_K();
    double G = compute_G();
    D_e = compute_isotropic_linear_elastic_matrix(K, G);
    Voigt delta_sigma_tilde = D_e*delta_epsilon_tilde;
    sigma_prime += delta_sigma_tilde.cauchy();
    std::cout << D_e << "\n";
    std::cout << delta_sigma_tilde << "\n";
    std::cout << delta_epsilon_tilde << "\n";
    std::cout << sigma_prime << "\n";
}

