#include "LinearElastic.hpp"
#include "LinearElastic_Definition.hpp"

LinearElastic::LinearElastic(std::vector<double> parameters, std::vector<double> state) : parameters(parameters), state(state) {
    set_name("LinearElastic");
    int parameters_required = 2;
    int state_required = 0;
    PLOG_FATAL_IF(parameters.size() != parameters_required) << parameters.size() << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state.size() != state_required) << state.size() << " parameters supplied when " << state_required << " expected.";
    assert(parameters.size() == parameters_required && state.size() == state_required);
    PLOG_INFO << name << " model instantiated with " << parameters.size() << " parameters and " << state.size() << " state variables."; 
    std::cout << K << "\n";
}

double LinearElastic::compute_K(double delta_epsilon_e_vol, double p_prime) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    if (delta_epsilon_e_vol != 0.0) {
        return BULK_MODULUS_SECANT;
    } else {
        return BULK_MODULUS_TANGENT;
    }
}

double LinearElastic::compute_G(double K) {
    using namespace std; /* Use std namespace for eye-pleasing model definitions. */
    return SHEAR_MODULUS;
}