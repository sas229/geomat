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
            // Correct stresses and state variables back to yield surface.
            PLOG_DEBUG << "Increment accepted. Checking for yield surface drift.";
            sigma_prime_u = sigma_prime_c = sigma_prime_ini;
            state_u = state_c = state_ini;

            // Correct stresses and state variables back to yield surface. 
            f_u = f_c = mf->compute_f(sigma_prime_u, state_u);
            int ITS_YSC = 0;
            while (std::abs(f_c) > settings->FTOL) {
                // If yield surface drift correction unsuccessful, log fault.
                if (ITS_YSC >= settings->MAXITS_YSC && std::abs(f_c) > settings->FTOL) {
                    PLOG_FATAL << "Maximum number of yield surface correction iterations performed and |f| = " << std::abs(f_c) << " > FTOL = " << settings->FTOL << ".";
                    assert(false);
                } 

                // Compute consistent correction to yield surface.
                compute_yield_surface_correction();
                PLOG_DEBUG << "f_c = " << f_c;

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

void Explicit::compute_yield_surface_correction(void) {
    // Calculate uncorrected elastic constitutive matrix using tangent moduli and elastic stress increment.
    Constitutive D_e_u = mf->compute_D_e(sigma_prime_u, Cauchy::Zero());

    // Calculate uncorrected derivatives.
    Cauchy df_dsigma_prime_u, dg_dsigma_prime_u;
    Voigt a_u, b_u;
    HardeningModuli H_s_u(state_ep.size());
    StateFactors B_s_u(state_ep.size());
    mf->compute_derivatives(sigma_prime_u, state_u, df_dsigma_prime_u, dg_dsigma_prime_u, H_s_u, B_s_u);
    a_u = to_voigt(df_dsigma_prime_u);
    b_u = to_voigt(dg_dsigma_prime_u);

    // Compute correction factor.
    double H_u = H_s_u.sum();
    double delta_lambda_c = f_u/(H_u + a_u.transpose()*D_e_u*b_u);

    // Update stress and state variables using corrections.
    Voigt Delta_sigma_prime_c = -delta_lambda_c*D_e_u*b_u;
    State Delta_state_c = -delta_lambda_c*B_s_u;
    sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
    state_c = state_u + Delta_state_c;

    // Check yield surface function value.
    f_c = mf->compute_f(sigma_prime_c, state_c);
    f_u = mf->compute_f(sigma_prime_u, state_u);
    if (std::abs(f_c) > std::abs(f_u)) {
        // Apply normal correction instead.
        double delta_lambda_c = f_u/(a_u.transpose()*a_u);
        
        //Update stress and state variables using correction.
        Voigt Delta_sigma_prime_c = -delta_lambda_c*a_u;
        sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
        state_c = state_u; /* i.e. no correction to state variables. */

        // Check yield surface function value.
        f_c = mf->compute_f(sigma_prime_c, state_c);
    }    
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