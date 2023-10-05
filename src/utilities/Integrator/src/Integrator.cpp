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
    PLOG_DEBUG << "Solver: " << settings->solver;
    PLOG_DEBUG << "Method: " << settings->method;
    if (settings->solver == std::string("Explicit")) {
        if (settings->method == std::string("ModifiedEuler")) {
            me.solve(sigma_prime, state, Delta_epsilon_tilde);
        } else if (settings->method == std::string("RKDP")) {
            rkdp.solve(sigma_prime, state, Delta_epsilon_tilde);
        } else {
            PLOG_FATAL << "Invalid method defined. This option is not yet implemented.";
            assert(false);
        }
    } else {
        PLOG_FATAL << "Invalid solver defined. This option is not yet implemented.";
        assert(false);
    } 
}