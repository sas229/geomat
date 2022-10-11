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

void Model::set_stress(Eigen::VectorXd s) {
    stress = s;
    stress(0) = 100.12; // Modified just to show the maps work - to be deleted later.
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

Eigen::VectorXd Model::get_stress(void) {
    return stress;
}

Eigen::MatrixXd Model::get_jacobian(void) {
    return jacobian;
}
