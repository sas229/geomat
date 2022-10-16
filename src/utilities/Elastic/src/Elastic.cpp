#include "Elastic.hpp"

void Elastic::compute_isotropic_linear_elastic_matrix(double E, double nu, Eigen::Matrix<double, 6, 6> &D_e) {
    // Check elastic paramaters are initialised.
    PLOG_ERROR_IF(std::isnan(E)) << "Young's modulus, E, not initialised.";
    PLOG_ERROR_IF(std::isnan(nu)) << "Poisson's ratio, nu, not initialised.";
    assert(!std::isnan(E) && !std::isnan(nu));

    // Fill elastic matrix with isotropic linear elastic coefficients.
    double C = E/((1.0+nu)*(1.0-2.0*nu));
    D_e(0,0) = D_e(1,1) = D_e(2,2) += C*(1.0-nu); 
    D_e(0,1) = D_e(0,2) = D_e(1,2) = D_e(1,0) = D_e(2,0) = D_e(2,1) += C*nu;
    D_e(3,3) = D_e(4,4) = D_e(5,5) += C*(1.0-2.0*nu)/2.0; 
    PLOG_INFO << "Isotropic linear elastic matrix computed.";
}

void Elastic::compute_anisotropic_linear_elastic_matrix(double E_h, double E_v, double G_h, double nu_v, double nu_h, Eigen::Matrix<double, 6, 6> &D_e) {
    // Check elastic paramaters are initialised.
    PLOG_ERROR_IF(std::isnan(E_h)) << "Young's modulus in vertical direction, E_h, not initialised.";
    PLOG_ERROR_IF(std::isnan(E_v)) << "Young's modulus in horizontal direction, E_v, not initialised.";
    PLOG_ERROR_IF(std::isnan(G_h)) << "Shear modulus in horizontal direction, G_h, not initialised.";
    PLOG_ERROR_IF(std::isnan(nu_v)) << "Poisson's ratio in vertical direction, nu_v, not initialised.";
    PLOG_ERROR_IF(std::isnan(nu_h)) << "Poisson's ratio in horizontal direction, nu_h, not initialised.";
    assert(!std::isnan(E_v) && !std::isnan(E_h) && !std::isnan(G_h) && !std::isnan(nu_v) && !std::isnan(nu_v));

    // Fill elastic matrix with anisotropic linear elastic coefficients.
    D_e(0,0) = D_e(1,1) += 1/E_h; 
    D_e(2,2) += 1/E_v; 
    D_e(0,1) = D_e(1,0) = -nu_h/E_h;
    D_e(2,0) = D_e(0,2) = D_e(2,1) = D_e(1,2) = -nu_v/E_v;
    D_e(3,3) = D_e(4,4) += 1/G_v;
    D_e(5,5) = 2.0*(1+nu_h)/E_h;
    PLOG_INFO << "Anisotropic linear elastic matrix computed.";
}

void Elastic::compute_simplified_anisotropic_linear_elastic_matrix(double alpha, double E_v, double nu_h, Eigen::Matrix<double, 6, 6> &D_e) {
    // Check elastic paramaters are initialised.
    PLOG_ERROR_IF(std::isnan(alpha)) << "Square root of ratio of Young's moduli, alpha, not initialised.";
    PLOG_ERROR_IF(std::isnan(E_v)) << "Young's modulus in horizontal direction, E_v, not initialised.";
    PLOG_ERROR_IF(std::isnan(nu_h)) << "Poisson's ratio in horizontal direction, nu_h, not initialised.";
    assert(!std::isnan(alpha) && !std::isnan(E_v) && !std::isnan(nu_h));

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

void Elastic::update_isotropic_linear_elastic_matrix(void) {
    Elastic::compute_isotropic_linear_elastic_matrix(this->E, this->nu, this->D_e);
}

void Elastic::update_anisotropic_linear_elastic_matrix(void) {
    Elastic::compute_anisotropic_linear_elastic_matrix(this->E_h, this->E_v, this->G_h, this->nu_h, this->nu_v, this->D_e);
}

void Elastic::update_simplified_anisotropic_linear_elastic_matrix(void) {
    Elastic::compute_simplified_anisotropic_linear_elastic_matrix(this->alpha, this->E_v, this->nu_h, this->D_e);
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
