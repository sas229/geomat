#include "Model.hpp"

void Model::set_name(std::string s) {
    name = s;
}

void Model::set_n_parameters(int i) {
    n_parameters = i;
}

void Model::set_n_state_variables(int i) {
    n_state_variables = i;
}

void Model::set_sigma(Eigen::VectorXd s) {
    sigma_v = s;
    // sigma(0) = 100.12; // Modified just to show the maps work - to be deleted later.
    
    // Assemble stress tensor for invariant calculation.
    sigma(0,0) = sigma_v(0);
    sigma(1,1) = sigma_v(1);
    sigma(2,2) = sigma_v(2);
    // The ordering for the terms below is specific to Abaqus. 
    // Consider something more conventional by adding an interface to 
    // convert for Abaqus usage.
    sigma(0,1) = sigma(1,0) = sigma_v(3); 
    sigma(0,2) = sigma(2,0) = sigma_v(4);
    sigma(1,2) = sigma(2,1) = sigma_v(5);
}

void Model::set_jacobian(Eigen::MatrixXd j) {
    jacobian = j;
    jacobian(0,2) = 10.99; // Modified just to show the maps work - to be deleted later.
}

std::string Model::get_name(void) {
    return name;
}

int Model::get_n_parameters(void) {
    return n_parameters;
}

int Model::get_n_state_variables(void) {
    return n_state_variables;
}

Eigen::VectorXd Model::get_sigma(void) {
    return sigma_v;
}

Eigen::MatrixXd Model::get_jacobian(void) {
    return jacobian;
}

std::vector<double> Model::get_state_variables(void) {
    return state;
}

void Model::compute_isotropic_linear_elastic_matrix(void) {
    // Fill elastic matrix with isotropic linear elastic coefficients.
    D_e(0,0) = D_e(1,1) = D_e(2,2) += K + 4.0/3.0*G; 
    D_e(0,1) = D_e(0,2) = D_e(1,2) = D_e(1,0) = D_e(2,0) = D_e(2,1) += K - 2.0/3.0*G;
    D_e(3,3) = D_e(4,4) = D_e (5,5) += G;
}

Eigen::Matrix<double, 6, 6> Model::get_elastic_matrix(void) {
    return D_e;
}

void Model::compute_G(void) {
    G = E/(2.0*(1.0+nu));
}

void Model::compute_E(void) {
    E = G*(2.0*(1.0+nu));
}

void Model::compute_K(void) {
    K = E/(3.0*(1.0-2.0*nu));
}

void Model::compute_invariants(void) {
    // Stress invariants.
    I_1 = sigma.trace();
    std::cout << "I_1 = " << I_1 << "\n";
    I_2 = 1.0/2.0*(std::pow(sigma.trace(),2) - ((sigma.array().pow(2)).matrix().trace()));
    std::cout << "I_2 = " << I_2 << "\n";
    I_2 = (sigma(1,1)*sigma(2,2)) + (sigma(2,2)*sigma(0,0)) + (sigma(0,0)*sigma(1,1));
    std::cout << "I_2 = " << I_2 << "\n";
    I_3 = sigma.determinant();
    std::cout << "I_3 = " << I_3 << "\n";
    
    // Mean stress.
    p = I_1/3.0;
    std::cout << "p = " << p << "\n";

    // Deviatoric stress tensor, s_t.
    s = (sigma.array()-(delta*p)).matrix();

    // Deviatoric stress invariants.
    J_1 = s.trace(); // Ought to be zero.
    std::cout << "J_1 = " << J_1 << "\n";
    J_2 = (std::pow(I_1,2)/3.0) - I_2;
    std::cout << "J_2 = " << J_2 << "\n";
    J_2 = 1.0/2.0*(((s.array()).pow(2)).matrix().trace());
    std::cout << "J_2 = " << J_2 << "\n";
    J_2 = 1.0/2.0*(std::pow(sigma(0,0)-p,2) + std::pow(sigma(1,1)-p,2) + std::pow(sigma(2,2)-p,2));
    std::cout << "J_2 = " << J_2 << "\n";
    J_3 = s.determinant();
    std::cout << "J_3 = " << J_3 << "\n";

    // Deviatoric stress.
    q = std::pow(3*J_2, 0.5);
    std::cout << "q = " << q << "\n";
    
    // Principal stresses.
    Eigen::VectorXcd principal_stresses = sigma.eigenvalues();
    sigma_1 = principal_stresses(0).real();
    sigma_2 = principal_stresses(1).real();
    sigma_3 = principal_stresses(2).real();
    std::cout << sigma_1 << "\n";
    std::cout << sigma_2 << "\n";
    std::cout << sigma_3 << "\n";

    mises = std::sqrt(3.0*J_2);
    std::cout << "mises = " << mises << "\n";

    max_shear = std::max({std::abs(sigma_1-sigma_2), std::abs(sigma_2-sigma_3), std::abs(sigma_1-sigma_3)})/2.0;
    std::cout << "max_shear = " << max_shear << "\n";
}



