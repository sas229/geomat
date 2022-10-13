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

void Model::set_sigma(Eigen::VectorXd v) {
    // Stress in Voigt notation form - change sign to use compression positive soil mechanics convention.
    sigma_prime_v = -v;
    
    // Assemble stress tensor for invariant calculation.
    sigma_prime(0,0) = sigma_prime_v(0);
    sigma_prime(1,1) = sigma_prime_v(1);
    sigma_prime(2,2) = sigma_prime_v(2);
    // The ordering for the terms below is specific to Abaqus. 
    // Consider something more conventional by adding an interface to 
    // convert for Abaqus usage.
    sigma_prime(0,1) = sigma_prime(1,0) = sigma_prime_v(3); 
    sigma_prime(0,2) = sigma_prime(2,0) = sigma_prime_v(4);
    sigma_prime(1,2) = sigma_prime(2,1) = sigma_prime_v(5);

    // Total stresses given pore pressure, u.
    sigma = (sigma_prime.array() + (delta.array()*u)).matrix();
}

void Model::set_jacobian(Eigen::MatrixXd m) {
    jacobian = m;
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
    // Change sign back to tension positive sign convention.
    return -sigma_prime_v;
}

Eigen::MatrixXd Model::get_jacobian(void) {
    return jacobian;
}

std::vector<double> Model::get_state_variables(void) {
    return state;
}

void Model::compute_stress_invariants(void) {
    // Stress invariants.
    I_1 = sigma.trace();
    I_2 = 1.0/2.0*(std::pow(sigma.trace(),2) - ((sigma.array().pow(2)).matrix().trace()));
    I_3 = sigma.determinant();
    
    // Mean stress.
    p = 1.0/3.0*sigma.trace();
    p_prime = 1.0/3.0*sigma_prime.trace();

    // Deviatoric stress tensor, s.
    s = (sigma.array()-(delta.array()*p)).matrix();

    // Deviatoric stress invariants.
    J_1 = s.trace(); // Is always zero...
    J_2 = (std::pow(I_1,2)/3.0) - I_2;
    J_3 = s.determinant();

    // Deviatoric stress using principal stresses.
    compute_principal_stresses();
    q = std::sqrt(1.0/2.0*((pow((sigma_1-sigma_2),2) + pow((sigma_2-sigma_3),2) + pow((sigma_1-sigma_3),2))));
    
    // von Mises criterion.
    mises = std::sqrt(3.0*J_2);
    
    // Maximum shear stress.
    max_shear = std::max({std::abs(sigma_1-sigma_2), std::abs(sigma_2-sigma_3), std::abs(sigma_1-sigma_3)})/2.0;
}

void Model::compute_principal_stresses(void) {
    // Solve eigenvalues and eigevectors..
    Eigen::EigenSolver<Eigen::MatrixXd> es(sigma_prime);
        
    // Principal stress magnitudes.
    Eigen::Vector3d principal_stresses = es.eigenvalues().real();
    S = principal_stresses.asDiagonal();

    // Sort principal stresses to allocate major, intermediate and minor appropriately.
    Eigen::Vector3d ordered_principal {principal_stresses.data()};
    std::sort(ordered_principal.data(), ordered_principal.data()+ordered_principal.size(), std::greater<int>());
    sigma_1 = ordered_principal(0);
    sigma_2 = ordered_principal(1);
    sigma_3 = ordered_principal(2);

    // Principal stress directions.
    T = es.eigenvectors().real();

    // Lode angle.
    theta_c = 1.0/3.0*std::acos(J_3/2.0*std::pow((3.0/J_2), 3.0/2.0));
    theta_s = pi/6.0 - theta_c;
    theta_s_bar = -theta_s;
}

void Model::compute_cartesian_stresses(void) {
    sigma = T*S*T.transpose();
}