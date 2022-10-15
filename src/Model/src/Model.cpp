#include "Model.hpp"

// Setters.

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

// Getters.

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

// Computers.

Eigen::Matrix3d Model::compute_cartesian_stresses(Eigen::Matrix3d T, Eigen::Matrix3d S) {
    return T*S*T.transpose();
}

void Model::compute_lode(double J_2, double J_3, double &theta_c, double &theta_s, double &theta_s_bar) {
    theta_c = 1.0/3.0*std::acos(J_3/2.0*std::pow((3.0/J_2), 3.0/2.0));
    theta_s = pi/6.0 - theta_c;
    theta_s_bar = -theta_s;
}

double Model::compute_max_shear(double sigma_1, double sigma_2, double sigma_3) {
    return std::max({std::abs(sigma_1-sigma_2), std::abs(sigma_2-sigma_3), std::abs(sigma_1-sigma_3)})/2.0;
}

double Model::compute_mises_stress(double J_2) {
    return std::sqrt(3.0*J_2);
}

double Model::compute_p(Eigen::Matrix3d sigma) {
    return 1.0/3.0*sigma.trace();
}

double Model::compute_p_prime(Eigen::Matrix3d sigma_prime) {
    return Model::compute_p(sigma_prime);
}

void Model::compute_principal_stresses(Eigen::Matrix3d sigma_prime, double &sigma_1, double &sigma_2, double &sigma_3, Eigen::Matrix3d &S, Eigen::Matrix3d &T) {
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
}

double Model::compute_q(Eigen::Matrix3d _sigma) {
    double q = std::sqrt(1.0/2.0*(
        (std::pow((_sigma(0,0)-_sigma(1,1)),2) + std::pow((_sigma(1,1)-_sigma(2,2)),2) + std::pow((_sigma(2,2)-_sigma(0,0)),2)) 
        + 6.0*(std::pow(_sigma(0,1),2) + std::pow(_sigma(0,2),2) + std::pow(_sigma(1,2),2))));
    return q;
}

Eigen::Matrix3d Model::compute_s(Eigen::Matrix3d _sigma, double _p) {
    return (_sigma.array()-(delta.array()*_p)).matrix();
}

void Model::compute_stress_invariants(Eigen::Matrix3d sigma, double &I_1, double &I_2, double &I_3, double &J_1, double &J_2, double &J_3) {
    // Stress invariants.
    I_1 = sigma.trace();
    I_2 = 1.0/2.0*(std::pow(sigma.trace(),2) - ((sigma.array().pow(2)).matrix().trace()));
    I_3 = sigma.determinant();

    // Mean stress.
    double p = compute_p_prime(sigma);
    
    // Deviatoric stress tensor, s.
    Eigen::Matrix3d s = compute_s(sigma, p);

    // Deviatoric stress invariants.
    J_1 = s.trace(); // Is always zero...
    J_2 = (std::pow(I_1,2)/3.0) - I_2;
    J_3 = s.determinant();
}

// Updaters.

void Model::update_cartesian_stresses(void) {
    this->sigma = Model::compute_cartesian_stresses(this->T, this->S);
}

void Model::update_lode(void) {
    Model::compute_lode(this->J_2, this->J_3, this->theta_c, this->theta_s, this->theta_s_bar);
}

void Model::update_max_shear(void) {
    this->max_shear = Model::compute_max_shear(this->sigma_1, this->sigma_2, this->sigma_3);
}

void Model::update_mises_stress(void) {
    this->mises_stress = Model::compute_mises_stress(this->J_2);
}

void Model::update_p(void) {
    this->p = Model::compute_p(this->sigma);
}

void Model::update_p_prime(void) {
    this->p_prime = Model::compute_p(this->sigma_prime);
}

void Model::update_principal_stresses(void) {
    Model::compute_principal_stresses(this->sigma_prime, this->sigma_1, this->sigma_2, this->sigma_3, this->S, this->T);
}   

void Model::update_q(void) {
    this->q = Model::compute_q(this->sigma);
}

void Model::update_s(void) {
    Model::compute_s(this->sigma, this->p);
}

void Model::update_stress_invariants(void) {
    Model::compute_stress_invariants(this->sigma, this->I_1, this->I_2, this->I_3, this->J_1, this->J_2, this->J_3);
}