#include "LinearElastic.hpp"
#include "LinearElastic_Definition.hpp"

LinearElastic::LinearElastic(Parameters parameters, State state) : parameters(parameters), state(state) {
    set_model_name("LinearElastic");
    int parameters_required = 2;
    int state_required = 0;
    PLOG_FATAL_IF(parameters.size() != parameters_required) << parameters.size() << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state.size() != state_required) << state.size() << " parameters supplied when " << state_required << " expected.";
    assert(parameters.size() == parameters_required && state.size() == state_required);
    PLOG_INFO << name << " model instantiated with " << parameters.size() << " parameters and " << state.size() << " state variables."; 
}

double LinearElastic::compute_K(double Delta_epsilon_e_vol, double p_prime) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    if (Delta_epsilon_e_vol != 0.0) {
        return LINEAR_ELASTIC_SECANT_BULK_MODULUS;
    } else {
        return LINEAR_ELASTIC_TANGENT_BULK_MODULUS;
    }
}

double LinearElastic::compute_G(double K) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    return LINEAR_ELASTIC_SHEAR_MODULUS;
}