#include "Elastic.hpp"

void Elastic::compute_isotropic_linear_elastic_matrix(void) {
    // Check elastic paramaters are initialised.
    PLOG_ERROR_IF(E == 0.0 || nu == 0.0) << "One or more elasticity parameters is zero.";

    // Fill elastic matrix with isotropic linear elastic coefficients.
    double C = E/((1.0+nu)*(1.0-2.0*nu));
    D_e(0,0) = D_e(1,1) = D_e(2,2) += C*(1.0-nu); 
    D_e(0,1) = D_e(0,2) = D_e(1,2) = D_e(1,0) = D_e(2,0) = D_e(2,1) += C*nu;
    D_e(3,3) = D_e(4,4) = D_e(5,5) += C*(1.0-2.0*nu)/2.0; 
    PLOG_INFO << "Isotropic linear elastic matrix computed.";
}

void Elastic::compute_anisotropic_linear_elastic_matrix(void) {
    // Check elastic paramaters are initialised.
    PLOG_ERROR_IF(E_v == 0.0 || E_h == 0.0 || G_v == 0.0 || nu_h == 0.0 || nu_v == 0.0) << "One or more elasticity parameters is zero.";

    // Fill elastic matrix with anisotropic linear elastic coefficients.
    D_e(0,0) = D_e(1,1) += 1/E_h; 
    D_e(2,2) += 1/E_v; 
    D_e(0,1) = D_e(1,0) = -nu_h/E_h;
    D_e(2,0) = D_e(0,2) = D_e(2,1) = D_e(1,2) = -nu_v/E_v;
    D_e(3,3) = D_e(4,4) += 1/G_v;
    D_e(5,5) = 2.0*(1+nu_h)/E_h;
    PLOG_INFO << "Anisotropic linear elastic matrix computed.";
}

void Elastic::compute_simplified_anisotropic_linear_elastic_matrix(void) {
    // Check elastic paramaters are initialised.
    PLOG_ERROR_IF(alpha == 0.0 || E_v == 0.0 || nu_h == 0.0) << "One or more elasticity parameters is zero.";

    // Fill elastic matrix with anisotropic linear elastic coefficients.
    double C = 1/E_v;
    D_e(0,0) = D_e(1,1) += C*1/std::pow(alpha,2); 
    D_e(2,2) = C*1;
    D_e(0,1) = D_e(1,0) = C*-nu_h/std::pow(alpha,2);
    D_e(2,0) = D_e(0,2) = D_e(2,1) = D_e(1,2) = C*-nu_h/alpha;
    D_e(3,3) = D_e(4,4) += C*2.0*(1.0+nu_h)/alpha;
    D_e(5,5) = C*2.0*(1.0+nu_h)/std::pow(alpha,2);
    PLOG_INFO << "Simplified anisotropic linear elastic matrix computed.";
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
