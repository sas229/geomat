#include "MCC.hpp"

MCC::MCC() {
    set_name("MCC");
    set_n_parameters(6);
    set_n_state_variables(2);
    PLOG_DEBUG << name << " model instantiated with " << n_parameters << " parameters and " << n_state_variables << " state variables.";
}

void MCC::set_state_variables(std::vector<double> s) {
    state = s;
}

void MCC::set_parameters(std::vector<double> p) {
    nu = p[0];
    M = p[1];
    N = p[2];
    lambda_star = p[3];
    kappa_star = p[4];
}

