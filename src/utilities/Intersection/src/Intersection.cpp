#include "Intersection.hpp"

double Intersection::compute_alpha(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde, double FTOL){
    PLOG_DEBUG << "Computing alpha.";
    double alpha;

    // Current stress state (i.e. alpha = 0.0).
    double f_0 = compute_f_alpha(sigma_prime, state, 0.0, Delta_epsilon_tilde);

    // Fully elastic increment trial stress state.
    double f_1 = compute_f_alpha(sigma_prime, state, 1.0, Delta_epsilon_tilde);

    PLOG_DEBUG << "f_0 = " << f_0 << " and f_1 = " << f_1;

    // Check increment type by finding alpha.
    if (f_1 <= FTOL) {
        PLOG_DEBUG << "Fully elastic increment; alpha = 1.0.";
        alpha = 1.0;
    } else if (std::abs(f_0) <= FTOL && f_1> FTOL) {
        PLOG_DEBUG << "Check potential for elastoplastic unloading-reloading.";
        if (!check_unload_reload(sigma_prime, state)) {
            // Compute bounds on alpha.
            double alpha_0 = 0.0;
            double alpha_1 = 1.0;
            compute_alpha_bounds(alpha_0, alpha_1);

            // Trial stress state for alpha_0.
            double f_0 = compute_f_alpha(sigma_prime, state, alpha_0, Delta_epsilon_tilde);

            // Trial stress state for alpha_1.
            double f_1 = compute_f_alpha(sigma_prime, state, alpha_1, Delta_epsilon_tilde);

            // Perform Pegasus intersection within the bounds of alpha_0 and alpha_1.
            PLOG_DEBUG << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = " << alpha_0 << " and alpha_1 = " << alpha_1 <<".";
            alpha = pegasus_regula_falsi(alpha_0, alpha_1, f_0, f_1);
        } else {
            PLOG_DEBUG << "Fully plastic increment; alpha = 0.0.";
            alpha = 0.0;
        } 
    } else if (f_0 < -FTOL && f_1 > FTOL) {
        PLOG_DEBUG << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = 0.0 and alpha_1 = 1.0.";
        alpha = pegasus_regula_falsi(0.0, 1.0, f_0, f_1);
    } else {
        PLOG_FATAL << "Illegal stress state.";
        assert(false);
    }
    return alpha;
}