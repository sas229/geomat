#include "LinearElastic.hpp"

LinearElastic::LinearElastic(Parameters parameters, State state) : parameters(parameters), state(state) {
    set_name("LinearElastic");
    int parameters_required = 2;
    int state_required = 0;
    PLOG_FATAL_IF(parameters.size() != parameters_required) << parameters.size() << " parameters supplied when " << parameters_required << " expected.";
    PLOG_FATAL_IF(state.size() != state_required) << state.size() << " parameters supplied when " << state_required << " expected.";
    assert(parameters.size() == parameters_required && state.size() == state_required);
    PLOG_INFO << name << " model instantiated with " << parameters.size() << " parameters and " << state.size() << " state variables."; 
}

double LinearElastic::compute_K(double Delta_epsilon_e_vol, double p_prime) {
    if (Delta_epsilon_e_vol != 0.0) {
        return K;
    } else {
        return K;
    }
}

double LinearElastic::compute_G(double K) {
    return G;
}