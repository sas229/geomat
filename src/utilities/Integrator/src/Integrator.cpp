#include "Integrator.hpp"

Integrator::Integrator(Settings *settings, ModelFunctions *mf) : me(settings, mf) {
    this->settings = settings;
    this->mf = mf;
}

void Integrator::solve(
    Cauchy &sigma_prime_ep, 
    State &state_ep, 
    Voigt Delta_epsilon_tilde_p
    ) {
    if (settings->solver == "Explicit") {
        if (settings->method == "ModifiedEuler") {
            me.solve(sigma_prime_ep, state_ep, Delta_epsilon_tilde_p);
        }
    } else {
        PLOG_FATAL << "Invalid solver defined. This option is not yet implemented.";
        assert(false);
    } 
}