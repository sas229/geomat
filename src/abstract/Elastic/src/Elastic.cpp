#include "Elastic.hpp"

Constitutive Elastic::compute_isotropic_linear_elastic_matrix(double K, double G) {
    // Check elastic paramaters are valid.
    Checks::check_elastic_parameters(K, G);

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

