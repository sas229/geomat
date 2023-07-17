#include "ModifiedEuler.hpp"

ModifiedEuler::ModifiedEuler(Settings *settings, ModelFunctions *mf) {
    this->settings = settings;
    this->mf = mf;
}

void ModifiedEuler::solve(
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
        compute_initial_estimate();
        if (accepted) {
            // Correct stresses and state variables back to yield surface.
            sigma_prime_u = sigma_prime_c = sigma_prime_ini;
            state_u = state_c = state_ini;

            // Correct stresses and state variables back to yield surface. 
            f_u = f_c = mf->compute_f(sigma_prime_u, state_u);
            int ITS_YSC = 0;
            while (ITS_YSC < settings->MAXITS_YSC && std::abs(f_c) > settings->FTOL) {
                // If yield surface drift correction unsuccessful, log fault.
                if (ITS_YSC >= settings->MAXITS_YSC && std::abs(f_c) > settings->FTOL) {
                    PLOG_FATAL << "Maximum number of yield surface correction iterations performed and f = " << f_c << " > FTOL = " << settings->FTOL << ".";
                    assert(false);
                } 

                // Compute consistent correction to yield surface.
                compute_yield_surface_correction();

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
        dT = compute_new_substep_size();
    }
    // Strain increment solved.
    PLOG_INFO << "Elastoplastic stress increment integrated to within a tolerance FTOL = " << settings->FTOL << " via " << substeps << " substeps and " << corrections << " drift corrections.";
    sigma_prime = sigma_prime_ep;
    state = state_ep;
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

void ModifiedEuler::compute_yield_surface_correction(void) {
    // Calculate uncorrected elastic constitutive matrix using tangent moduli and elastic stress increment.
    Constitutive D_e_u = mf->compute_D_e(sigma_prime_u, Cauchy::Zero());

    // Calculate uncorrected derivatives.
    Cauchy df_dsigma_prime_u, dg_dsigma_prime_u;
    Voigt a_u, b_u;
    double H_u;
    mf->compute_derivatives(sigma_prime_u, state_u, df_dsigma_prime_u, a_u, dg_dsigma_prime_u, b_u, H_u);

    // Compute correction factor.
    double delta_lambda_c = f_u/(H_u + a_u.transpose()*D_e_u*b_u);

    // Update stress and state variables using corrections.
    Voigt Delta_sigma_prime_c = -delta_lambda_c*D_e_u*b_u;
    State Delta_state_c = mf->compute_state_increment(delta_lambda_c, df_dsigma_prime_u, H_u, Voigt::Zero());
    sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
    state_c = state_u + Delta_state_c;

    // Check yield surface function value.
    double f_c = mf->compute_f(sigma_prime_c, state_c);
    f_u = mf->compute_f(sigma_prime_u, state_u);
    if (std::abs(f_c) > std::abs(f_u)) {
        // Apply normal correction instead.
        double delta_lambda_c = f_u/(a_u.transpose()*a_u);
        
        //Update stress and state variables using correction.
        Voigt Delta_sigma_prime_c = -delta_lambda_c*a_u;
        sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
        state_c = state_u; /* i.e. no correction to state variables. */
    }
}

double ModifiedEuler::compute_new_substep_size(void) {
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