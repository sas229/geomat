#include "Checks.hpp"

void Checks::check_inputs(std::string name, int parameters_size, int state_size, int parameters_required, int state_required) {
    PLOG_FATAL_IF(parameters_size != parameters_required) << parameters_size << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state_size != state_required) << state_size << " parameters supplied when " << state_required << " expected.";
    assert(parameters_size == parameters_required && state_size == state_required);
    PLOG_INFO << name << " model instantiated with " << parameters_size << " parameters and " << state_size << " state variables.";  
}

void Checks::check_elastic_parameters(double K, double G) {
    PLOG_FATAL_IF(std::isnan(K)) << "Bulk modulus is " << K << ".";
    PLOG_FATAL_IF(std::isnan(G)) << "Shear modulus is " << G << ".";
    PLOG_FATAL_IF(K <= 0.0) << "Bulk modulus less than or equal to zero.";
    PLOG_FATAL_IF(G <= 0.0) << "Shear modulus less than or equal to zero.";
    assert(K > 0.0 && G > 0.0);
}