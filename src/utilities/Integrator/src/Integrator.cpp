#include "Integrator.hpp"

Integrator::Integrator(Settings *settings, ModelFunctions *mf) : me(settings, mf), rkdp(settings, mf) {
    this->settings = settings;
    this->mf = mf;
}

void Integrator::solve(
    Cauchy &sigma_prime, 
    State &state, 
    Voigt Delta_epsilon_tilde
    ) {
    if (settings->solver == "Explicit") {
        if (settings->method == "ModifiedEuler") {
            me.solve(sigma_prime, state, Delta_epsilon_tilde);
        } else if (settings->method == "RKDP") {
            rkdp.solve(sigma_prime, state, Delta_epsilon_tilde);
        }
    } else {
        PLOG_FATAL << "Invalid solver defined. This option is not yet implemented.";
        assert(false);
    } 
}