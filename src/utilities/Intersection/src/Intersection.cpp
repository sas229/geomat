#include "Intersection.hpp"

double Intersection::compute_alpha(
    Cauchy sigma_prime,
    State state,
    Voigt Delta_epsilon_tilde,
    double FTOL,
    double LTOL,
    int MAXITS_YSI,
    int NSUB,
    YieldFunction compute_f,
    TrialFunction compute_trial_stress,
    ConstitutiveMatrixFunction compute_D_e,
    DerivativeFunction compute_derivatives) {
    // Current stress state (i.e. alpha = 0.0).
    double f_0 = compute_f(sigma_prime, state);

    // Compute trial stress state with full strain increment.
    Cauchy sigma_prime_1 = compute_trial_stress(sigma_prime, 1.0 * Delta_epsilon_tilde);
    double f_1 = compute_f(sigma_prime_1, state);

    // Check increment type by finding alpha.
    if (f_1 <= FTOL) {
        PLOG_INFO << "Fully elastic increment; alpha = 1.0.";
        return 1.0;
    } else if (std::abs(f_0) <= FTOL && f_1 > FTOL) {
        PLOG_INFO << "Check potential for elastoplastic unloading-reloading.";
        if (!Intersection::check_unload_reload(sigma_prime, state, Delta_epsilon_tilde, LTOL, compute_f, compute_D_e, compute_derivatives)) {
            // // Compute bounds on alpha.
            double alpha_0 = 0.0;
            double alpha_1 = 1.0;
            Intersection::compute_alpha_bounds(sigma_prime, state, Delta_epsilon_tilde, FTOL, MAXITS_YSI, NSUB, compute_f, compute_D_e, compute_trial_stress, alpha_0, alpha_1);

            // Trial stress state for alpha_0.
            Cauchy sigma_prime_0 = compute_trial_stress(sigma_prime, alpha_0 * Delta_epsilon_tilde);
            double f_0 = compute_f(sigma_prime_0, state);

            // Trial stress state for alpha_1.
            Cauchy sigma_prime_1 = compute_trial_stress(sigma_prime, alpha_1 * Delta_epsilon_tilde);
            double f_1 = compute_f(sigma_prime_1, state);

            // Perform Pegasus intersection within the bounds of alpha_0 and alpha_1.
            PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = " << alpha_0 << " and alpha_1 = " << alpha_1 << ".";
            return Intersection::pegasus_regula_falsi(sigma_prime, state, Delta_epsilon_tilde, alpha_0, alpha_1, f_0, f_1, FTOL, MAXITS_YSI, compute_f, compute_trial_stress);
        } else {
            PLOG_INFO << "Fully plastic increment; alpha = 0.0.";
            return 0.0;
        }
    } else if (f_0 < -FTOL && f_1 > FTOL) {
        PLOG_INFO << "Plastic increment: finding alpha via the Pegasus algorithm using alpha_0 = 0.0 and alpha_1 = 1.0.";
        return Intersection::pegasus_regula_falsi(sigma_prime, state, Delta_epsilon_tilde, 0.0, 1.0, f_0, f_1, FTOL, MAXITS_YSI, compute_f, compute_trial_stress);
    } else {
        PLOG_FATAL << "Illegal stress state.";
        assert(false);
    }
    return 0.0;
}

bool Intersection::check_unload_reload(
    Cauchy sigma_prime,
    State state,
    Voigt Delta_epsilon_tilde,
    double LTOL,
    YieldFunction compute_f,
    ConstitutiveMatrixFunction compute_D_e,
    DerivativeFunction compute_derivatives) {
    // Compute required derivatives for the given stress state.
    Cauchy df_dsigma_prime_check, dg_dsigma_prime_check;
    Voigt a_check, b_check;
    double H_check;
    compute_derivatives(sigma_prime, state, df_dsigma_prime_check, a_check, dg_dsigma_prime_check, b_check, H_check);

    // Compute the elastic matrix using tangent moduli.
    Constitutive D_e_tan = compute_D_e(sigma_prime, Cauchy::Zero());

    // Compute elastic stress increment.
    Voigt Delta_sigma_e = D_e_tan * Delta_epsilon_tilde;

    // Check unloading-reloading criterion.
    double cos_theta = (double)(a_check.transpose() * Delta_sigma_e) / (double)(a_check.squaredNorm() * Delta_sigma_e.squaredNorm());

    // Check against tolerance.
    if (cos_theta >= -LTOL) {
        return false;
    } else {
        return true;
    }
}

void Intersection::compute_alpha_bounds(
    Cauchy sigma_prime,
    State state,
    Voigt Delta_epsilon_tilde,
    double FTOL,
    int MAXITS_YSI,
    int NSUB,
    YieldFunction compute_f,
    ConstitutiveMatrixFunction compute_D_e,
    TrialFunction compute_trial_stress,
    double &alpha_0,
    double &alpha_1) {
    double alpha_n, d_alpha;
    int i = 0;
    int j = 0;
    while (i < MAXITS_YSI) {
        d_alpha = (alpha_1 - alpha_0) / NSUB;
        alpha_n = alpha_0 + d_alpha;
        while (j < NSUB) {
            // Compute the elastic matrix.
            Cauchy Delta_epsilon = to_cauchy(Delta_epsilon_tilde);
            Constitutive D_e_trial = compute_D_e(sigma_prime, alpha_n * Delta_epsilon);

            // Compute elastic stress increment.
            Cauchy sigma_prime_trial = compute_trial_stress(sigma_prime, alpha_n * Delta_epsilon_tilde);

            // Check yield function.
            double f_trial = compute_f(sigma_prime_trial, state);

            // Check criterion.
            if (f_trial > FTOL) {
                // Break from loops.
                alpha_1 = alpha_n;
                return;
            } else {
                // Continue iterating.
                alpha_0 = alpha_n;
                alpha_n = alpha_n + d_alpha;
                j += 1;
            }
        }
        i += 1;
    }
    PLOG_FATAL << "Maximum number of iterations MAXITS_YSI = " << MAXITS_YSI << " therefore failed to compute valid bounds on alpha.";
    assert(false);
}

double Intersection::pegasus_regula_falsi(
    Cauchy sigma_prime,
    State state,
    Voigt Delta_epsilon_tilde,
    double alpha_0,
    double alpha_1,
    double f_0,
    double f_1,
    double FTOL,
    int MAXITS_YSI,
    YieldFunction compute_f,
    TrialFunction compute_trial_stress) {
    // Pegasus algorithm.
    double alpha, alpha_n;
    double f_n = FTOL;
    Cauchy sigma_prime_n = sigma_prime;
    State state_n = state;

    // Iterate to find optimal alpha if a plastic increment.
    double ITS_YSI = 0;
    while (ITS_YSI < MAXITS_YSI && std::abs(f_n) >= FTOL) {
        alpha_n = alpha_1 - f_1 * (alpha_1 - alpha_0) / (f_1 - f_0);

        sigma_prime_n = compute_trial_stress(sigma_prime, alpha_n * Delta_epsilon_tilde);
        f_n = compute_f(sigma_prime_n, state_n);

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
    if (std::abs(f) >= FTOL) {
        PLOG_FATAL << "Performed " << MAXITS_YSI << " Pegasus iteration(s): alpha = " << alpha << "; |f| = " << std::abs(f) << " > tolerance = " << FTOL << ".";
        assert(false);
    } else {
        PLOG_INFO << "Performed " << ITS_YSI << " Pegasus iteration(s): "
                  << "alpha = " << alpha << "; "
                  << "|f| = " << std::abs(f) << " < tolerance = " << FTOL << ".";
        assert(alpha >= 0.0 && alpha <= 1.0);
    }
    return alpha;
}