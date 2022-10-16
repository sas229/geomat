#include "MCC.hpp"

MCC::MCC() {
    set_name("MCC");
    set_n_parameters(6);
    set_n_state_variables(2);
    PLOG_INFO << name << " model instantiated with " << n_parameters << " parameters and " << n_state_variables << " state variables.";
    nu = 0.3;
    E = 50000.0;
    update_isotropic_linear_elastic_matrix();
}

void MCC::set_state_variables(std::vector<double> state) {
    this->state = state;
}

void MCC::set_parameters(std::vector<double> parameters) {
    nu = parameters[0];
    M = parameters[1];
    N = parameters[2];
    lambda_star = parameters[3];
    kappa_star = parameters[4];
}

