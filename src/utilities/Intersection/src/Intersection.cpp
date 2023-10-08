#include "Intersection.hpp"

Intersection::Intersection(Settings *settings, ModelFunctions *mf) {
    this->settings = settings;
    this->mf = mf;
}

double Intersection::solve(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde) {  
    // Store current model state.
    this->sigma_prime = sigma_prime;
    this->state = state;
    this->Delta_epsilon_tilde = Delta_epsilon_tilde;

    // Reset bounds on alpha.
    alpha_0 = 0.0;
    alpha_1 = 1.0;

    // Yield function value for current stress state (i.e. alpha = 0.0).
    f_0 = mf->compute_f(sigma_prime, state);
    PLOG_DEBUG << "f_0 = " << f_0;

    // Compute trial stress state with full strain increment.
    Cauchy sigma_prime_1, Delta_sigma_prime_1;
    State state_1(state.size());
    State Delta_state_1(state.size());
    mf->compute_trial_increment(sigma_prime, state, 1.0 * Delta_epsilon_tilde, Delta_sigma_prime_1, Delta_state_1);
    sigma_prime_1 = sigma_prime + Delta_sigma_prime_1;
    state_1 = state + Delta_state_1;
    f_1 = mf->compute_f(sigma_prime_1, state_1);
    PLOG_DEBUG << "f_1 = " << f_1;

    // Check increment type by finding alpha.
    if (f_1 <= settings->FTOL) {
        PLOG_INFO << "Fully elastic increment; alpha = 1.0";
        return 1.0;
    } else if (std::abs(f_0) <= settings->FTOL && f_1 > settings->FTOL) {
        PLOG_INFO << "Check potential for elastoplastic unloading-reloading.";
        if (!check_unload_reload()) {
            // Refine bounds on alpha.
            refine_alpha_bounds();

            // Perform Pegasus intersection within the bounds of alpha_0 and alpha_1.
            PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = " << alpha_0 << " and alpha_1 = " << alpha_1 << ".";
            return pegasus_regula_falsi();
        } else {
            PLOG_INFO << "Fully plastic increment; alpha = 0.0.";
            return 0.0;
        }
    } else if (f_0 < -settings->FTOL && f_1 > settings->FTOL) {
        PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = 0.0 and alpha_1 = 1.0.";
        return pegasus_regula_falsi();
    } else {
        PLOG_FATAL << "Illegal stress state.";
        assert(false);
    }
    return 0.0;
}
    
bool Intersection::check_unload_reload(void) {
    // Compute required derivatives for the given stress state.
    Cauchy df_dsigma_prime_c, dg_dsigma_prime_c;
    Voigt a_c, b_c;
    HardeningModuli H_s_c(state.size());
    StateFactors B_s_c(state.size());
    Derivatives derivatives = mf->compute_derivatives(sigma_prime, state);
    df_dsigma_prime_c = derivatives.df_dsigma_prime;
    dg_dsigma_prime_c = derivatives.dg_dsigma_prime;
    a_c = to_voigt(df_dsigma_prime_c);
    b_c = to_voigt(dg_dsigma_prime_c);

    // Compute the elastic matrix using tangent moduli.
    Constitutive D_e_tan = mf->compute_D_e(sigma_prime, Cauchy::Zero());

    // Compute elastic stress increment.
    Voigt Delta_sigma_prime_e = D_e_tan * Delta_epsilon_tilde;

    // Check unloading-reloading criterion.
    double cos_theta = (double)(a_c.transpose() * Delta_sigma_prime_e) / (double)(a_c.squaredNorm() * Delta_sigma_prime_e.squaredNorm());

    // Check against tolerance.
    if (cos_theta >= -settings->LTOL) {
        PLOG_INFO << "No unload-reload in current increment.";
        return false;
    } else {
        PLOG_INFO << "Elastoplastic unload-reload increment.";
        return true;
    }
}

void Intersection::refine_alpha_bounds(void) {
    double alpha_n, d_alpha;
    int i = 0;
    int j = 0;
    while (alpha_n <= alpha_1) {
        d_alpha = (alpha_1 - alpha_0) / settings->NSUB;
        alpha_n = alpha_0 + d_alpha;
        while (j < settings->NSUB) {
            // Compute the elastic matrix.
            Cauchy Delta_epsilon = to_cauchy(Delta_epsilon_tilde);

            // Compute elastic stress increment.
            Cauchy sigma_prime_n, Delta_sigma_prime_n;
            State state_n(state.size());
            State Delta_state_n(state.size());
            mf->compute_trial_increment(sigma_prime, state, alpha_n * Delta_epsilon_tilde, Delta_sigma_prime_n, Delta_state_n);
            sigma_prime_n = sigma_prime + Delta_sigma_prime_n;
            state_n = state + state_n;

            // Check yield function.
            double f_n = mf->compute_f(sigma_prime_n, state_n);

            // Check criterion.
            if (f_n > settings->FTOL) {
                // Intersection found.
                alpha_1 = alpha_n;
                
                // Update yield function value for current lower bound stress state (i.e. alpha = alpha_0).
                Cauchy sigma_prime_0, Delta_sigma_prime_0;
                State state_0, Delta_state_0;
                mf->compute_trial_increment(sigma_prime, state, alpha_0 * Delta_epsilon_tilde, Delta_sigma_prime_0, Delta_state_0);
                sigma_prime_0 = sigma_prime + Delta_sigma_prime_0;
                state_0 = state + Delta_state_0;
                f_0 = mf->compute_f(sigma_prime_0, state_0);

                // Update yield function value for current upper bound stress state (i.e. alpha = alpha_1).
                Cauchy sigma_prime_1, Delta_sigma_prime_1;
                State state_1, Delta_state_1;
                mf->compute_trial_increment(sigma_prime, state, alpha_1 * Delta_epsilon_tilde, Delta_sigma_prime_1, Delta_state_1);
                sigma_prime_1 = sigma_prime + Delta_sigma_prime_1;
                state_1 = state + Delta_state_1;
                f_1 = mf->compute_f(sigma_prime_1, state_1);
                
                // Break from loops.
                return;
            } else {
                // Continue iterating.
                alpha_0 = alpha_n;
                alpha_n += d_alpha;
                j += 1;
            }
        }
    }
    PLOG_FATAL << "Transition not found therefore failed to compute valid bounds on alpha.";
    assert(false);
}

double Intersection::pegasus_regula_falsi(void) {
    // Pegasus algorithm.
    double alpha_n;
    double f_n = settings->FTOL;
    Cauchy sigma_prime_n, Delta_sigma_prime_n;
    State state_n(state.size());
    State Delta_state_n(state.size());

    // Iterate to find optimal alpha if a plastic increment.
    double ITS_YSI = 0;
    while (ITS_YSI < settings->MAXITS_YSI && std::abs(f_n) >= settings->FTOL) {
        alpha_n = alpha_1 - f_1 * (alpha_1 - alpha_0) / (f_1 - f_0);

        mf->compute_trial_increment(sigma_prime, state, alpha_n * Delta_epsilon_tilde, Delta_sigma_prime_n, Delta_state_n);
        sigma_prime_n = sigma_prime + Delta_sigma_prime_n;
        state_n = state + state_n;
        f_n = mf->compute_f(sigma_prime_n, state_n);

        // Update trial using Pegasus method rules.
        if (f_n * f_1 < 0) {
            alpha_1 = alpha_0;
            f_1 = f_0;
        } else {
            f_1 = f_1 * f_0 / (f_0 + f_n);
        }
        alpha_0 = alpha_n;
        f_0 = f_n;
        ITS_YSI += 1;
    }
    alpha = alpha_n;
    alpha = std::max(alpha, 0.0);
    alpha = std::min(alpha, 1.0);
    double f = f_n;
    if (std::abs(f) >= settings->FTOL) {
        PLOG_FATAL << "Performed " << settings->MAXITS_YSI << " Pegasus iteration(s): alpha = " << alpha << "; |f| = " << std::abs(f) << " > tolerance = " << settings->FTOL << ".";
        assert(false);
    } else {
        PLOG_INFO << "Performed " << ITS_YSI << " Pegasus iteration(s): "
                  << "alpha = " << alpha << "; "
                  << "|f| = " << std::abs(f) << " < tolerance = " << settings->FTOL << ".";
        assert(alpha >= 0.0 && alpha <= 1.0);
    }
    return alpha;
}