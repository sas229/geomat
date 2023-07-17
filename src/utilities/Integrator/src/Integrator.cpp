#include "Integrator.hpp"

Integrator::Integrator(Settings settings, ModelFunctions mf) {
    this->settings = settings;
    this->mf = mf;
}

void Integrator::solve(
    Cauchy &sigma_prime_ep, 
    State &state_ep, 
    Voigt Delta_epsilon_tilde_p
) {
    if (settings.method == "Explicit") {
        // Sloan et al. (2001) substepping with automatic error control.
        int substeps = 0;
        int corrections = 0;
        double dT = 1.0;
        double T = 0.0;
        while (T < 1.0) {
            PLOG_DEBUG << "Plastic increment, dT = " << dT;
            Voigt Delta_epsilon_tilde_p_dT = Delta_epsilon_tilde_p*dT;
            double Delta_epsilon_vol_p_dT = to_cauchy(Delta_epsilon_tilde_p_dT).trace();

            // Compute stress and state variable increment estimates.
            Voigt Delta_sigma_prime_1, Delta_sigma_prime_2;
            State delta_state_1, delta_state_2;
            Explicit::compute_plastic_increment(sigma_prime_ep, state_ep, Delta_epsilon_tilde_p_dT, Delta_sigma_prime_1, delta_state_1, mf);
            Cauchy sigma_prime_1 = sigma_prime_ep + to_cauchy(Delta_sigma_prime_1);
            State state_1 = state_ep + delta_state_1;
            Explicit::compute_plastic_increment(sigma_prime_1, state_ep, Delta_epsilon_tilde_p_dT, Delta_sigma_prime_2, delta_state_2, mf);
            PLOG_DEBUG << "State variable increment 1, delta_state_1 = " << delta_state_1;
            PLOG_DEBUG << "State variable increment 2, delta_state_2 = " << delta_state_2;

            // Calculate modified Euler stresses and state variables.
            Cauchy sigma_prime_ini = sigma_prime_ep + to_cauchy(1.0/2.0*(Delta_sigma_prime_1 + Delta_sigma_prime_2));
            State state_ini = state_ep + 1.0/2.0*(delta_state_1+delta_state_2);
            PLOG_DEBUG << "Initial stress estimate, sigma_prime_ini_tilde = " << to_voigt(sigma_prime_ini);
            PLOG_DEBUG << "State variable estimate, state_ini = " << state_ini;

            // Compute error estimate.
            bool accepted = false;
            double R_n = Explicit::compute_error_estimate(sigma_prime_ini, state_ini, Delta_sigma_prime_1, Delta_sigma_prime_2, delta_state_1, delta_state_2, settings);
            if (R_n < settings.STOL) {
                // Accept increment and correct stresses and state variables back to yield surface.
                accepted = true;
                Cauchy sigma_prime_u, sigma_prime_c;
                sigma_prime_u = sigma_prime_c = sigma_prime_ini;
                State state_u, state_c;
                state_u = state_c = state_ini;

                // Correct stresses and state variables back to yield surface. 
                double f_u, f_c;
                f_u = f_c = mf.compute_f(sigma_prime_u, state_u);
                int ITS_YSC = 0;
                while (ITS_YSC < settings.MAXITS_YSC && std::abs(f_c) > settings.FTOL) {
                    // If yield surface drift correction unsuccessful, log fault.
                    if (ITS_YSC >= settings.MAXITS_YSC && std::abs(f_c) > settings.FTOL) {
                        PLOG_FATAL << "Maximum number of yield surface correction iterations performed and f = " << f_c << " > FTOL = " << settings.FTOL << ".";
                        assert(false);
                    } 

                    // Calculate uncorrected elastic constitutive matrix using tangent moduli and elastic stress increment.
                    Constitutive D_e_u = mf.compute_D_e(sigma_prime_u, Cauchy::Zero());

                    // Calculate uncorrected derivatives.
                    Cauchy df_dsigma_prime_u, dg_dsigma_prime_u;
                    Voigt a_u, b_u;
                    double H_u;
                    mf.compute_derivatives(sigma_prime_u, state_u, df_dsigma_prime_u, a_u, dg_dsigma_prime_u, b_u, H_u);

                    // Compute consistent correction to yield surface.
                    Explicit::compute_yield_surface_correction(sigma_prime_u, state_u, f_u, H_u, a_u, b_u, D_e_u, sigma_prime_c, state_c, mf);

                    // Check yield surface function.
                    f_c = mf.compute_f(sigma_prime_c, state_c);
                    f_u = mf.compute_f(sigma_prime_u, state_u);
                    if (std::abs(f_c) > std::abs(f_u)) {
                        // Apply normal correction instead.
                        Explicit::compute_normal_yield_surface_correction(sigma_prime_u, state_u, f_u, a_u, sigma_prime_c, state_c);
                    }

                    // Correct stress and state variables.
                    sigma_prime_u = sigma_prime_c;
                    state_u = state_c;
                    ITS_YSC += 1;
                }

                // Update stress state and state variables.
                corrections += ITS_YSC;
                sigma_prime_ep = sigma_prime_c;
                state_ep = state_c;

                // Increment substep and pseudotime T.
                substeps += 1;
                T += dT;
            }
            dT = Explicit::compute_new_substep_size(accepted, dT, T, R_n, settings);
        }
        // Strain increment solved.
        PLOG_INFO << "Elastoplastic stress increment integrated to within a tolerance FTOL = " << settings.FTOL << " via " << substeps << " substeps and " << corrections << " drift corrections.";
    }
}