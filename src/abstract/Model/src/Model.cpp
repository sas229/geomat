#include "Model.hpp"

Model::Model(std::string log_severity) {
    // Initialise log.
    initialise_log(log_severity);
}

bool Model::check_inputs(std::string name, int parameters_size, int state_size, int parameters_required, int state_required) {
    PLOG_FATAL_IF(parameters_size != parameters_required) << parameters_size << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state_size != state_required) << state_size << " state variables supplied when " << state_required << " expected.";
    if (parameters_size == parameters_required && state_size == state_required) {
        PLOG_INFO << name << " model instantiated with " << parameters_size << " parameters and " << state_size << " state variables.";  
        return true;
    } else {
        assert(false);
        throw std::invalid_argument("Incorrect number of parameters or state variables supplied.");
        return false;
    }
}

// Setters.

void Model::set_model_name(std::string s) {
    name = s;
}

void Model::set_model_type(std::string t) {
    model_type = t;
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

Cauchy Model::compute_dq_dsigma_prime(Cauchy sigma_prime, Cauchy s, double q) {
    if (q == 0.0) {
        return Cauchy::Zero();
    } else {
        Cauchy dq_dsigma_prime;
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
        return dq_dsigma_prime;
    }
}

Cauchy Model::compute_dsigma_bar_dsigma_prime(Cauchy sigma_prime, Cauchy s, double sigma_bar) {
    if (sigma_bar == 0.0) {
        return Cauchy::Zero();
    } else {
        Cauchy dsigma_bar_dsigma_prime;
        dsigma_bar_dsigma_prime(0,0) = s(0,0);
        dsigma_bar_dsigma_prime(1,1) = s(1,1);
        dsigma_bar_dsigma_prime(2,2) = s(2,2);
        dsigma_bar_dsigma_prime(0,1) = 2*sigma_prime(0,1);
        dsigma_bar_dsigma_prime(1,0) = 2*sigma_prime(1,0);
        dsigma_bar_dsigma_prime(0,2) = 2*sigma_prime(0,2);
        dsigma_bar_dsigma_prime(2,0) = 2*sigma_prime(2,0);
        dsigma_bar_dsigma_prime(1,2) = 2*sigma_prime(1,2);
        dsigma_bar_dsigma_prime(2,1) = 2*sigma_prime(2,1);
        dsigma_bar_dsigma_prime *= 1.0/(2.0*sigma_bar);
        return dsigma_bar_dsigma_prime;
    }
}

void Model::set_sigma_prime_tilde(Voigt sigma_prime_tilde) {
    // Stress in Voigt notation form - using compression positive soil mechanics convention.
    this->sigma_prime_tilde = sigma_prime_tilde;
    this->sigma_prime = to_cauchy(this->sigma_prime_tilde);

    // Total stresses given pore pressure, u.
    sigma = compute_sigma(sigma_prime, u);

    // Stress invariants and other stress measures.
    compute_stress_invariants(sigma, I_1, I_2, I_3, J_1, J_2, J_3);
    compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);
    compute_principal_stresses(sigma_prime, sigma_1, sigma_2, sigma_3, R, S);
    p_prime = compute_p_prime(sigma_prime);
    p = compute_p(sigma);
    q = compute_q(sigma_prime);
    mises_stress = compute_mises_stress(J_2);
    max_shear = compute_max_shear(sigma_1, sigma_2, sigma_3);
}

void Model::set_Delta_epsilon_tilde(Voigt Delta_epsilon_tilde) {
    // Strain increment in Voigt notation form - change sign to use compression positive soil mechanics convention.
    this->Delta_epsilon_tilde = Delta_epsilon_tilde;
    this->Delta_epsilon = to_cauchy(Delta_epsilon_tilde);
    solved = false;   
}

// Getters.

std::string Model::get_model_name(void) {
    return name;
}

std::string Model::get_model_type(void) {
    return model_type;
}

Voigt Model::get_sigma_prime_tilde(void) {
    return sigma_prime_tilde;
}

Jacobian Model::get_jacobian(void) {
    return jacobian;
}

// Compute methods.

Cauchy Model::compute_cartesian_stresses(Cauchy R, Cauchy S) {
    return R*S*R.transpose();
}

double Model::compute_Delta_epsilon_vol(Cauchy Delta_epsilon) {
    return Delta_epsilon.trace();
}

void Model::compute_lode(double J_2, double J_3, double &theta_c, double &theta_s, double &theta_s_bar) {
    double val;
    if (J_2 != 0) {
        val = J_3/2.0*std::pow((3.0/J_2), 3.0/2.0); // Must be in the range of +/- 1.0.
        val = std::min(val, 1.0);
        val = std::max(val, -1.0);
    } else {
        val = 0.0;
    }
    theta_c = (1.0/3.0)*std::acos(val);
    theta_s = (1.0/3.0)*std::asin(val);
    theta_s_bar = -theta_s;
}

double Model::compute_max_shear(double sigma_1, double sigma_2, double sigma_3) {
    return std::max({std::abs(sigma_1-sigma_2), std::abs(sigma_2-sigma_3), std::abs(sigma_1-sigma_3)})/2.0;
}

double Model::compute_mises_stress(double J_2) {
    return std::sqrt(3.0*J_2);
}

double Model::compute_p(Cauchy sigma) {
    return 1.0/3.0*sigma.trace();
}

double Model::compute_p_prime(Cauchy sigma_prime) {
    return Model::compute_p(sigma_prime);
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
}

double Model::compute_q(Cauchy sigma) {
    double q = std::sqrt(1.0/2.0*(
        (std::pow((sigma(0,0)-sigma(1,1)),2) + std::pow((sigma(1,1)-sigma(2,2)),2) + std::pow((sigma(2,2)-sigma(0,0)),2)) 
        + 6.0*(std::pow(sigma(0,1),2) + std::pow(sigma(0,2),2) + std::pow(sigma(1,2),2))));
    return q;
}

Cauchy Model::compute_s(Cauchy sigma, double p) {
    return sigma - p*eye;
}

Cauchy Model::compute_sigma(Cauchy sigma_prime, double u) {
    return sigma_prime + u*eye;
}

void Model::compute_stress_invariants(Cauchy sigma, double &I_1, double &I_2, double &I_3, double &J_1, double &J_2, double &J_3) {
    // Stress invariants.
    I_1 = sigma.trace();
    I_2 = 1.0/2.0*(std::pow(sigma.trace(),2) - (sigma.cwiseProduct(sigma)).trace());
    I_3 = sigma.determinant();

    // Mean stress.
    double p = compute_p(sigma);
    
    // Deviatoric stress tensor, s.
    Cauchy s = compute_s(sigma, p);

    // Deviatoric stress invariants.
    J_1 = s.trace(); // Is always zero...
    J_2 = (std::pow(I_1,2)/3.0) - I_2;
    J_3 = s.determinant();
}

double Model::compute_sigma_bar(double J_2) {
    return std::sqrt(J_2);
}

Cauchy Model::compute_dJ_3_dsigma_prime(Cauchy sigma_prime, Cauchy s, double sigma_bar) {
    dJ_3_dsigma_prime(0,0) = (s(1,1)*s(2,2) - std::pow(sigma_prime(1,2),2.0)) + (1.0*std::pow(sigma_bar,2)/3.0);
    dJ_3_dsigma_prime(1,1) = (s(0,0)*s(2,2) - std::pow(sigma_prime(0,2),2.0)) + (1.0*std::pow(sigma_bar,2)/3.0);
    dJ_3_dsigma_prime(2,2) = (s(0,0)*s(1,1) - std::pow(sigma_prime(0,1),2.0)) + (1.0*std::pow(sigma_bar,2)/3.0);
    dJ_3_dsigma_prime(0,1) = 2.0*(sigma_prime(1,2)*sigma_prime(0,2) - s(2,2)*sigma_prime(0,1));
    dJ_3_dsigma_prime(0,2) = 2.0*(sigma_prime(0,2)*sigma_prime(0,1) - s(0,0)*sigma_prime(1,2));
    dJ_3_dsigma_prime(1,2) = 2.0*(sigma_prime(0,1)*sigma_prime(1,2) - s(1,1)*sigma_prime(0,2));
    dJ_3_dsigma_prime(1,0) = dJ_3_dsigma_prime(0,1);
    dJ_3_dsigma_prime(2,0) = dJ_3_dsigma_prime(0,2);
    dJ_3_dsigma_prime(2,1) = dJ_3_dsigma_prime(1,2);
    return dJ_3_dsigma_prime;
}

void Model::initialise_log(std::string severity) {
    // If no logger found, initialise the logger.
    auto log = plog::get();
    if (log == NULL) {
        static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
        if (severity == "verbose") {
            plog::init(plog::verbose, &consoleAppender);
        } else if (severity == "debug") {
            plog::init(plog::debug, &consoleAppender);
        } else if (severity == "info") {
            plog::init(plog::info, &consoleAppender);
        } else if (severity == "warning") {
            plog::init(plog::warning, &consoleAppender);
        } else if (severity == "error") {
            plog::init(plog::error, &consoleAppender);
        } else if (severity == "fatal") {
            plog::init(plog::fatal, &consoleAppender);
        } else if (severity == "none") {
            plog::init(plog::none, &consoleAppender);
        } else {
            plog::init(plog::warning, &consoleAppender);
            PLOG_WARNING << "Unrecognised log severity set. Possible options include \"verbose\", \"debug\", \"info\", \"warning\", \"error\", \"fatal\", or \"none\". Defaulting to \"error\".";
            plog::get()->setMaxSeverity(plog::error);
        }
    } else {
        // Otherwise adjust the severity.
        log->setMaxSeverity(plog::severityFromString(severity.c_str()));
    }
}