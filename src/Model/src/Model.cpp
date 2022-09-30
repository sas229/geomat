#include "Model.hpp"

void Model::set_name(std::string name) {
    _name = name;
}

void Model::set_nparams(int nparams) {
    _nparams = nparams;
}

void Model::set_nstatev(int nstatev) {
    _nstatev = nstatev;
}

void Model::set_stress(Eigen::VectorXd stress) {
    _ntens = stress.size();
    _stress = stress;
    _stress(0) = 100.12; // Modified just to show the maps work - to be deleted later.
}

std::string Model::get_name(void) {
    return _name;
}

int Model::get_nparams(void) {
    return _nparams;
}

int Model::get_nstatev(void) {
    return _nstatev;
}

Eigen::VectorXd Model::get_stress(void) {
    return _stress;
}
