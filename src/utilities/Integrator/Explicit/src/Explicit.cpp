#include "Explicit.hpp"

void Explicit::solve(
    Cauchy &sigma_prime, 
    State &state, 
    Voigt Delta_epsilon_tilde 
) { 
    // Store effective stress and state variables locally.
    sigma_prime_ep  = sigma_prime;
    state_ep = state;
    
    // Sloan et al. (2001) substepping with automatic error control.
    int substeps = 0;
    int corrections = 0;
    dT = 1.0;
    T = 0.0;
    while (T < 1.0) {
        // Compute initial estimate and error estimate.
        PLOG_DEBUG << "Plastic increment, dT = " << dT;
        Delta_epsilon_tilde_dT = Delta_epsilon_tilde*dT;
        PLOG_DEBUG << "Computing initial stress increment estimate.";
        compute_initial_estimate();
        if (accepted) {
            // Increment accepted; update stress and state variables.
            PLOG_DEBUG << "Increment accepted.";
            sigma_prime_u = sigma_prime_c = sigma_prime_ini;
            state_u = state_c = state_ini;

            // Correct stresses and state variables back to yield surface, if required. 
            int ITS_YSC = 0;
            mf->check_yield_surface_drift(sigma_prime_u, state_u, sigma_prime_c, state_c, ITS_YSC);
            corrections += ITS_YSC;

            // Update stress state and state variables.
            sigma_prime_ep = sigma_prime_c;
            state_ep = state_c;

            // Increment substep and pseudotime T.
            substeps += 1;
            T += dT;
        }
        dT = compute_new_substep_size();
    }
    // Strain increment solved.
    PLOG_INFO << "Elastoplastic stress increment integrated to within a tolerance FTOL = " << settings->FTOL << " via " << substeps << " substeps and " << corrections << " drift corrections.";
    sigma_prime = sigma_prime_ep;
    state = state_ep;
}

double Explicit::compute_new_substep_size(void) {
    if (dT == settings->DT_MIN) {
        PLOG_FATAL << "Minimum step size DT_MIN = " << settings->DT_MIN << " failed to generated an accepted increment.";
        assert(false);
    }
    double q_step;
    if (accepted) {
        // Step size accepted: allow substep size to grow.
        q_step = std::min(0.9*std::pow(settings->STOL/R_n, 1.0/order), 1.1);
    } else {
        // Step size rejected: calculate reduced substep size factor and limit substep size growth factor.
        q_step = std::max(0.9*std::pow(settings->STOL/R_n, 1.0/order), 0.1);
        q_step = std::min(q_step, 1.0);
    }
    dT *= q_step;

    // Delimit substep size to pseudomtime remaining and minimum substep size.
    dT = std::max(dT, settings->DT_MIN);
    dT = std::min(dT, 1.0-T);
    return dT;
}