#include "LinearElastic.hpp"

LinearElastic::LinearElastic(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state) {
    set_model_name("LinearElastic");
    set_model_type("Elastic");

    // Initialise logger.
    initialise_log(log_severity);

    // Check inputs.
    check_inputs(get_model_name(), (int)parameters.size(), (int)state.size(), parameters_required, state_required);
}

Constitutive LinearElastic::compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon) {
    D_e = compute_isotropic_linear_elastic_matrix(K, G);   
    return D_e;
}