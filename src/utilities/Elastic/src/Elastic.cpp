#include "Elastic.hpp"

void Elastic::compute_isotropic_linear_elastic_matrix(void) {
    // Fill elastic matrix with isotropic linear elastic coefficients.
    D_e(0,0) = D_e(1,1) = D_e(2,2) += K + 4.0/3.0*G; 
    D_e(0,1) = D_e(0,2) = D_e(1,2) = D_e(1,0) = D_e(2,0) = D_e(2,1) += K - 2.0/3.0*G;
    D_e(3,3) = D_e(4,4) = D_e (5,5) += G;
}

Eigen::Matrix<double, 6, 6> Elastic::get_elastic_matrix(void) {
    return D_e;
}

void Elastic::compute_G_given_E_and_nu(void) {
    G = E/(2.0*(1.0+nu));
}

void Elastic::compute_E_given_G_and_nu(void) {
    E = G*(2.0*(1.0+nu));
}

void Elastic::compute_K_given_E_and_nu(void) {
    K = E/(3.0*(1.0-2.0*nu));
}
