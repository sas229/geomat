#include "Model.hpp"

// Setters.
void Model::set_name(std::string name) {
    _name = name;
}

void Model::set_nparams(int nparams) {
    _nparams = nparams;
    PLOG_DEBUG << _name << " model instantiated with " << _nparams << " parameters.";
}

void Model::set_stress(Vector6d stress) {
    _stress = stress;
    _stress(0) = 100.12;
}

// Getters.
std::string Model::get_name(void) {
    return _name;
}

int Model::get_nparams(void) {
    return _nparams;
}

Vector6d Model::get_stress(void) {
    return _stress;
}
