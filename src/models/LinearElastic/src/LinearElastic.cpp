#include "LinearElastic.hpp"

LinearElastic::LinearElastic(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state) {
    set_model_name("LinearElastic");
    set_model_type("Elastic");

    // Initialise logger.
    Logging::initialise_log(log_severity);

    // Check inputs.
    Checks::check_inputs(name, parameters.size(), state.size(), parameters_required, state_required);
}

State LinearElastic::get_state_variables(void) {
    return state;
}

void LinearElastic::set_state_variables(State new_state) {
    state = new_state;
}

Constitutive LinearElastic::compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon) {

    /* USER DEFINED CODE STARTS HERE */
    D_e = compute_isotropic_linear_elastic_matrix(K, G);
    /* USER DEFINED CODE ENDS HERE */
    
    return D_e;
}