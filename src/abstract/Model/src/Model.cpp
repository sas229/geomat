#include "Model.hpp"

void Model::set_sigma_prime_tilde(Voigt sigma_prime_tilde) {
    // Stress in Voigt notation form - change sign to use compression positive soil mechanics convention.
    this->sigma_prime_tilde = -sigma_prime_tilde;
    this->sigma_prime = to_cauchy(this->sigma_prime_tilde);

    //Debug log output.
    PLOG_INFO << "Initial stress state set.";
    PLOG_DEBUG << "sigma_prime_tilde = " << this->sigma_prime_tilde;
    PLOG_DEBUG << "sigma_prime = \n" << this->sigma_prime;

    // Compute stress variables.
    compute_stress_variables();
}

void Model::set_Delta_epsilon_tilde(Voigt Delta_epsilon_tilde) {
    // Strain increment in Voigt notation form - change sign to use compression positive soil mechanics convention.
    this->Delta_epsilon_tilde = -Delta_epsilon_tilde;
    this->Delta_epsilon = to_cauchy(this->Delta_epsilon_tilde);
    solved = false;

    //Debug log output.
    PLOG_INFO << "Strain increment set.";
    PLOG_DEBUG << "Delta_epsilon_tilde = " << this->Delta_epsilon_tilde;
    PLOG_DEBUG << "Delta_epsilon = \n" << this->Delta_epsilon;
}

void Model::compute_stress_variables(void) {
    // Total stresses given pore pressure, u.
    sigma = compute_sigma(sigma_prime, u);

    // Mean stress.
    p = compute_p(sigma);
    p_prime = compute_p_prime(sigma_prime);
    
    // Deviatoric stress tensor, s.
    s = compute_s(sigma, p);

    // Stress invariants and other stress measures.
    compute_stress_invariants(sigma, p, s, I_1, I_2, I_3, J_1, J_2, J_3);
    
    q = compute_q(sigma_prime);
    mises_stress = compute_mises_stress(J_2);
    max_shear = compute_max_shear(sigma_1, sigma_2, sigma_3);

    // Lode andlge and principal stresses.
    compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
    compute_principal_stresses(sigma_prime, sigma_1, sigma_2, sigma_3, R, S);
}

Cauchy Model::compute_dq_dsigma_prime(Cauchy sigma_prime, Cauchy s, double q) {
    Cauchy dq_dsigma_prime = Cauchy::Zero();
    if (q != 0.0) {
        dq_dsigma_prime(0,0) = s(0,0);
        dq_dsigma_prime(1,1) = s(1,1);
        dq_dsigma_prime(2,2) = s(2,2);
        dq_dsigma_prime(0,1) = 2*sigma_prime(0,1);
        dq_dsigma_prime(1,0) = 2*sigma_prime(1,0);
        dq_dsigma_prime(0,2) = 2*sigma_prime(0,2);
        dq_dsigma_prime(2,0) = 2*sigma_prime(2,0);
        dq_dsigma_prime(1,2) = 2*sigma_prime(1,2);
        dq_dsigma_prime(2,1) = 2*sigma_prime(2,1);
        dq_dsigma_prime *= 3.0/(2.0*q);
    }

    // Debug output.
    PLOG_DEBUG << "dq_dsigma_prime = \n" << dq_dsigma_prime;
    
    return dq_dsigma_prime;
}

Voigt Model::to_voigt(Cauchy cauchy) {
    Voigt voigt = Voigt::Zero();
    voigt(0) = cauchy(0,0);
    voigt(1) = cauchy(1,1);
    voigt(2) = cauchy(2,2);
    voigt(3) = cauchy(0,1);
    voigt(4) = cauchy(0,2);
    voigt(5) = cauchy(1,2);
    return voigt;
}

Cauchy Model::to_cauchy(Voigt voigt) {
    Cauchy cauchy = Cauchy::Zero();
    cauchy(0,0) = voigt(0);
    cauchy(1,1) = voigt(1);
    cauchy(2,2) = voigt(2);
    cauchy(0,1) = cauchy(1,0) = voigt(3);
    cauchy(0,2) = cauchy(2,0) = voigt(4);
    cauchy(1,2) = cauchy(2,1) = voigt(5);
    return cauchy;
}

Cauchy Model::compute_cartesian_stresses(Cauchy R, Cauchy S) {
    Cauchy sigma_prime = R*S*R.transpose();

    // Debug output.
    PLOG_DEBUG << "sigma_prime = \n" << sigma_prime;

    return sigma_prime;
}

double Model::compute_Delta_epsilon_vol(Cauchy Delta_epsilon) {
    double Delta_epsilon_vol = Delta_epsilon.trace();

    // Debug output.
    PLOG_DEBUG << "Delta_epsilon_vol = " << Delta_epsilon_vol;

    return Delta_epsilon_vol;
}

void Model::compute_lode(double J_2, double J_3, double &theta_c, double &theta_s, double &theta_s_bar) {
    // Note: need to carefully check these definitions.
    if (J_2 == 0.0 && J_3 == 0.0) {
        theta_c = 0.0;
    } else {
        theta_c = 1.0/3.0*std::acos(J_3/2.0*std::pow((3.0/J_2), 3.0/2.0));
    }
    theta_s = pi/6.0 - theta_c;
    theta_s_bar = -theta_s;

    // Debug output.
    PLOG_DEBUG << "Lode's angles: theta_c = " << theta_c << "; theta_s = " << theta_s << "; theta_s_bar = " << theta_s_bar;
}

double Model::compute_max_shear(double sigma_1, double sigma_2, double sigma_3) {
    double max_shear = std::max({std::abs(sigma_1-sigma_2), std::abs(sigma_2-sigma_3), std::abs(sigma_1-sigma_3)})/2.0;

    // Debug output.
    PLOG_DEBUG << "max_shear = " << max_shear;
    
    return max_shear;
}

double Model::compute_mises_stress(double J_2) {
    double mises_stress = std::sqrt(3.0*J_2);

    // Debug output.
    PLOG_DEBUG << "mises_stress = " << mises_stress;

    return mises_stress;
}

double Model::compute_p(Cauchy sigma) {
    double p = 1.0/3.0*sigma.trace();

    // Debug output.
    LOG_DEBUG << "p = " << p;

    return p;
}

double Model::compute_p_prime(Cauchy sigma_prime) {
    double p_prime = 1.0/3.0*sigma_prime.trace();

    // Debug output.
    PLOG_DEBUG << "p_prime = " << p_prime;

    return p_prime;
}

void Model::compute_principal_stresses(Cauchy sigma_prime, double &sigma_1, double &sigma_2, double &sigma_3, Cauchy &R, Cauchy &S) {
    // Solve eigenvalues and eigevectors..
    Eigen::EigenSolver<Eigen::MatrixXd> es(sigma_prime);
        
    // Principal stress magnitudes.
    Eigen::Vector3d principal_stresses = es.eigenvalues().real();
    S(0,0) = principal_stresses(0);
    S(1,1) = principal_stresses(1);
    S(2,2) = principal_stresses(2);

    // Sort principal stresses to allocate major, intermediate and minor appropriately.
    Eigen::Vector3d ordered_principal {principal_stresses.data()};
    std::sort(ordered_principal.data(), ordered_principal.data()+ordered_principal.size(), std::greater<int>());
    sigma_1 = ordered_principal(0);
    sigma_2 = ordered_principal(1);
    sigma_3 = ordered_principal(2);

    // Principal stress directions.
    R = es.eigenvectors().real();

    // Debug output.
    PLOG_DEBUG << "Principal stresses: sigma_1 = " << sigma_1 << "; sigma_2 = " << sigma_2 << "; sigma_3 = " << sigma_3;
    PLOG_DEBUG << "R = \n" << R; 
    PLOG_DEBUG << "S = \n" << S; 
}

double Model::compute_q(Cauchy sigma) {
    double q = std::sqrt(1.0/2.0*(
        (std::pow((sigma(0,0)-sigma(1,1)),2) + std::pow((sigma(1,1)-sigma(2,2)),2) + std::pow((sigma(2,2)-sigma(0,0)),2)) 
        + 6.0*(std::pow(sigma(0,1),2) + std::pow(sigma(0,2),2) + std::pow(sigma(1,2),2))));

    // Debug output.
    PLOG_DEBUG << "q = " << q;

    return q;
}

Cauchy Model::compute_s(Cauchy sigma, double p) {
    Cauchy s = sigma - p*eye;

    // Debug log output.
    PLOG_DEBUG << "s = \n" << s;

    return s;
}

Cauchy Model::compute_sigma(Cauchy sigma_prime, double u) {
    Cauchy sigma = sigma_prime + u*eye;

    // Debug log output.
    PLOG_DEBUG << "sigma = \n" << sigma;

    return sigma;
}

void Model::compute_stress_invariants(Cauchy sigma, double p, Cauchy s, double &I_1, double &I_2, double &I_3, double &J_1, double &J_2, double &J_3) {
    // Stress invariants.
    I_1 = sigma.trace();
    I_2 = 1.0/2.0*(std::pow(sigma.trace(),2) - (sigma.cwiseProduct(sigma)).trace());
    I_3 = sigma.determinant();

    // Deviatoric stress invariants.
    J_1 = s.trace(); // Is always zero...
    J_2 = (std::pow(I_1,2)/3.0) - I_2;
    J_3 = s.determinant();

    //Debug log output.
    PLOG_DEBUG << "I_1 = " << I_1 << "; I_2 = " << I_2 << "; I_3 = " << I_3 << "; J_1 = " << J_1 << "; J_2 = " << J_2 << "; J_3 = " << J_3;
}

Cauchy Model::compute_dtheta_dsigma_prime(Cauchy sigma_prime) {
    // Requires implementation for Lode angle dependent models.
    Cauchy dtheta_dsigma_prime;
    dtheta_dsigma_prime(0,0) = 0.0;
    dtheta_dsigma_prime(1,0) = 0.0;
    dtheta_dsigma_prime(2,0) = 0.0;
    dtheta_dsigma_prime(0,1) = 0.0;
    dtheta_dsigma_prime(1,1) = 0.0;
    dtheta_dsigma_prime(2,1) = 0.0;
    dtheta_dsigma_prime(0,2) = 0.0;
    dtheta_dsigma_prime(1,2) = 0.0;
    dtheta_dsigma_prime(2,2) = 0.0;

    // Debug output.
    PLOG_DEBUG << "dtheta_dsigma_prime = \n" << dtheta_dsigma_prime;

    return dtheta_dsigma_prime;
}

Cauchy Model::compute_dJ_3_dsigma_prime(Cauchy sigma_prime, Cauchy s, double q) {
    Cauchy dJ_3_dsigma_prime;
    dJ_3_dsigma_prime(0,0) = (s(1,1)*s(2,2) - std::pow(sigma_prime(1,2),2)) + (1.0*std::pow(q,2)/3.0);
    dJ_3_dsigma_prime(1,1) = (s(0,0)*s(2,2) - std::pow(sigma_prime(0,2),2)) + (1.0*std::pow(q,2)/3.0);
    dJ_3_dsigma_prime(2,2) = (s(0,0)*s(1,1) - std::pow(sigma_prime(0,1),2)) + (1.0*std::pow(q,2)/3.0);
    dJ_3_dsigma_prime(0,1) = 2.0*(sigma_prime(1,2)*sigma_prime(0,2) - s(2,2)*sigma_prime(0,1));
    dJ_3_dsigma_prime(0,2) = 2.0*(sigma_prime(0,2)*sigma_prime(0,1) - s(0,0)*sigma_prime(1,2));
    dJ_3_dsigma_prime(1,2) = 2.0*(sigma_prime(0,1)*sigma_prime(1,2) - s(1,1)*sigma_prime(0,2));
    dJ_3_dsigma_prime(1,0) = dJ_3_dsigma_prime(0,1);
    dJ_3_dsigma_prime(2,0) = dJ_3_dsigma_prime(0,2);
    dJ_3_dsigma_prime(2,1) = dJ_3_dsigma_prime(1,2);

    // Debug output.
    PLOG_DEBUG << "dJ_3_dsigma_prime = \n" << dJ_3_dsigma_prime;

    return dJ_3_dsigma_prime;
}

void Model::set_name(std::string s) {
    name = s;
}

void Model::set_model_type(std::string t) {
    model_type = t;
}

void Model::set_IP_number(int n) {
    IP_number = n;
}

int Model::get_IP_number(void) {
    return IP_number;
}

double Model::get_p_prime(void) {
    return p_prime;
}

double Model::get_q(void) {
    return q;
}

double Model::get_I_1(void) {
    return I_1;
}

double Model::get_I_2(void) {
    return I_2;
}

double Model::get_I_3(void) {
    return I_3;
}

double Model::get_J_1(void) {
    return J_1;
}

double Model::get_J_2(void) {
    return J_2;
}

double Model::get_J_3(void) {
    return J_3;
}

bool Model::get_solved(void) {
    return solved;
}

double Model::get_mises_stress(void) {
    return mises_stress;
}

double Model::get_max_shear(void) {
    return max_shear;
}

std::string Model::get_name(void) {
    return name;
}

std::string Model::get_model_type(void) {
    return model_type;
}

Voigt Model::get_sigma_prime_tilde(void) {
    // Note: change sign back to tension positive sign convention.
    return -sigma_prime_tilde;
}

Cauchy Model::get_sigma_prime(void) {
    // Note: change sign back to tension positive sign convention.
    return -sigma_prime;
}

Jacobian Model::get_jacobian(void) {
    return jacobian;
}