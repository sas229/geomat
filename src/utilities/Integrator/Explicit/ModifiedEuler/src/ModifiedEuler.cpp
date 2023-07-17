#include "ModifiedEuler.hpp"

ModifiedEuler::ModifiedEuler(Settings *settings, ModelFunctions *mf) {
    this->settings = settings;
    this->mf = mf;

    // Set order of method.
    this->order = 2.0;
}

void ModifiedEuler::compute_initial_estimate(void) {
    // Calculate estimate of increment.
    Voigt Delta_sigma_prime_1, Delta_sigma_prime_2;
    State Delta_state_1, Delta_state_2;
    mf->compute_plastic_increment(sigma_prime_ep, state_ep, Delta_epsilon_tilde_dT, Delta_sigma_prime_1, Delta_state_1);
    Cauchy sigma_prime_1 = sigma_prime_ep + to_cauchy(Delta_sigma_prime_1);
    State state_1 = state_ep + Delta_state_1;
    mf->compute_plastic_increment(sigma_prime_1, state_ep, Delta_epsilon_tilde_dT, Delta_sigma_prime_2, Delta_state_2);
    PLOG_DEBUG << "State variable increment 1, Delta_state1 = " << Delta_state_1;
    PLOG_DEBUG << "State variable increment 2, Delta_state2 = " << Delta_state_2;

    // Calculate modified Euler stresses and state variables.
    sigma_prime_ini = sigma_prime_ep + to_cauchy(1.0/2.0*(Delta_sigma_prime_1 + Delta_sigma_prime_2));
    state_ini = state_ep + 1.0/2.0*(Delta_state_1+Delta_state_2);
    PLOG_DEBUG << "Initial stress estimate, sigma_prime_ini_tilde = " << to_voigt(sigma_prime_ini);
    PLOG_DEBUG << "State variable estimate, state_ini = " << state_ini;

    // Compute error estimate.
    int size_state = state_ini.size();
    State error(size_state);
    error[0] = (to_cauchy(Delta_sigma_prime_2 - Delta_sigma_prime_1)).norm()/(2.0*sigma_prime_ini.norm());
    for (int i=1; i<size_state; i++) {
        error[i] = std::abs((Delta_state_2[i] - Delta_state_1[i]))/(2.0*state_ini[i]);
    }
    R_n = error.maxCoeff();

    // Check error estimate against machine tolerance.
    R_n = std::max(R_n, settings->EPS);

    // Check acceptance of estimate.
    if (R_n < settings->STOL) {
        accepted = true;
    } else {
        accepted = false;
    }
}