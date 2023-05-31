#include "LinearElastic.hpp"
#include "LinearElastic_Definition.hpp"

LinearElastic::LinearElastic(Parameters parameters, State state, std::string log_severity) : parameters(parameters), state(state) {
    set_model_name("LinearElastic");
    set_model_type("Elastic");

    // Initialise logger.
    Logging::initialise_log(log_severity);

    // Check inputs.
    Checks::check_inputs(name, parameters.size(), state.size(), parameters_required, state_required);
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