#include "Elastoplastic.hpp"

Elastoplastic::Elastoplastic(std::string log_severity) : intersection(&settings, &mf), integrator(&settings, &mf) {
    // Initialise log.
    initialise_log(log_severity);

    // Define binds to model functions.
    using namespace std::placeholders;
    mf.compute_f = std::bind(&Elastoplastic::compute_f, this, _1, _2);
    mf.compute_trial_increment = std::bind(&Elastoplastic::compute_elastic_increment, this, _1, _2, _3, _4, _5);
    mf.compute_D_e = std::bind(&Elastoplastic::compute_D_e, this, _1, _2);
    mf.compute_derivatives = std::bind(&Elastoplastic::compute_derivatives, this, _1, _2);
    mf.compute_plastic_increment = std::bind(&Elastoplastic::compute_plastic_increment, this, _1, _2, _3, _4, _5);
    mf.check_yield_surface_drift = std::bind(&Elastoplastic::check_yield_surface_drift, this, _1, _2, _3, _4, _5);
}

void Elastoplastic::solve(void) {
    if (!solved) {
        // Get the current state variables.
        State state = get_state_variables();
        
        // Initial conditions.
        PLOG_DEBUG << "Strain increment, Delta_epsilon_tilde = " << Delta_epsilon_tilde;
        PLOG_DEBUG << "Initial stress state, sigma_prime_tilde = " << to_voigt(sigma_prime);
        PLOG_DEBUG << "Initial state variables, state = " << state;
        
        // Compute alpha using bound functions.
        double alpha = intersection.solve(sigma_prime, state, Delta_epsilon_tilde);
        Cauchy sigma_prime_e, sigma_prime_ep;
        State state_e, state_ep;
        Voigt Delta_epsilon_tilde_e = alpha*Delta_epsilon_tilde;
        if (alpha == 0.0) {
            // Fully plastic increment. Update stress and state variables.
            sigma_prime_e = sigma_prime;
            state_e = state;
        }
        if (alpha > 0.0) {
            // Update stress based on elastic portion of strain increment.
            Voigt Delta_epsilon_tilde_e = alpha*Delta_epsilon_tilde;
            Cauchy Delta_epsilon_e = to_cauchy(Delta_epsilon_tilde_e);
            D_e = compute_D_e(sigma_prime, Delta_epsilon_e);
            sigma_prime_e = sigma_prime + to_cauchy(D_e*alpha*Delta_epsilon_tilde);
            state_e = state;    // No elastic change in state variables.
        } 
        PLOG_DEBUG << "Elastic strain increment, Delta_epsilon_e = " << Delta_epsilon_tilde_e;
        PLOG_DEBUG << "Stress state after elastic increment, sigma_prime_e_tilde = " << to_voigt(sigma_prime_e);
        PLOG_DEBUG << "State variables after elastic increment, state_e = " << state_e;
        if (alpha < 1.0) {
            // Perform elastoplastic stress integration on plastic portion of strain increment.
            PLOG_INFO << "Plastic increment; alpha = " << alpha;
            Voigt Delta_epsilon_tilde_p = (1.0-alpha)*Delta_epsilon_tilde;
            sigma_prime_ep = sigma_prime_e; 
            state_ep = state_e;
            PLOG_DEBUG << "Plastic strain increment, Delta_epsilon_tilde_p = " << Delta_epsilon_tilde_p;
            PLOG_DEBUG << "Initial elastoplastic stress state, sigma_prime_ep_tilde = " << to_voigt(sigma_prime_ep);
            PLOG_DEBUG << "Initial elastoplastic state variables, state_ep = " << state_ep;
            integrator.solve(sigma_prime_ep, state_ep, Delta_epsilon_tilde_p);
            sigma_prime = sigma_prime_ep;
            set_state_variables(state_ep);
        } else {
            // Fully elastic increment. Update stress and state variables.
            sigma_prime = sigma_prime_e;
            set_state_variables(state_e);
            PLOG_INFO << "Fully elastic stress increment integrated directly.";
        }
        PLOG_DEBUG << "Final stress state, sigma_prime_tilde = " << to_voigt(sigma_prime);
        PLOG_DEBUG << "Final state variables, state = " << state;

        // Compute final stress invariants.
        p_prime = compute_p_prime(sigma_prime);
        q = compute_q(sigma_prime);
        compute_principal_stresses(sigma_prime, sigma_1, sigma_2, sigma_3, R, S);
        sigma_prime_tilde = to_voigt(sigma_prime);
        solved = true;
    }
}

void Elastoplastic::check_yield_surface_drift(Cauchy sigma_prime_u, State state_u, Cauchy &sigma_prime_c, State &state_c, int &ITS_YSC) {
    PLOG_DEBUG << "Checking for yield surface drift.";
    double f_u, f_c;
    f_u = f_c = compute_f(sigma_prime_u, state_u);
    ITS_YSC = 0;
    while (std::abs(f_c) > settings.FTOL) {
        // If yield surface drift correction unsuccessful, log fault.
        if (ITS_YSC >= settings.MAXITS_YSC && std::abs(f_c) > settings.FTOL) {
            PLOG_FATAL << "Maximum number of " << settings.MAXITS_YSC << " yield surface correction iterations performed and |f| = " << std::abs(f_c) << " > FTOL = " << settings.FTOL << ".";
            assert(false);
            throw std::range_error("Maximum number of yield surface correction iterations performed.");
        } else {
            ITS_YSC++;
        }

        // Calculate uncorrected elastic constitutive matrix using tangent moduli and elastic stress increment.
        Constitutive D_e_u = compute_D_e(sigma_prime_u, Cauchy::Zero());

        // Calculate uncorrected derivatives.
        Derivatives derivatives = compute_derivatives(sigma_prime_u, state_u);
        Cauchy df_dsigma_prime_u = derivatives.df_dsigma_prime;
        Cauchy dg_dsigma_prime_u = derivatives.dg_dsigma_prime;
        Voigt a_u = to_voigt(df_dsigma_prime_u);
        Voigt b_u = to_voigt(dg_dsigma_prime_u);
        HardeningModuli H_s_u = derivatives.H_s;
        StateFactors B_s_u = derivatives.B_s;

        // Compute correction plastic multiplier.
        double H_u = H_s_u.sum();
        double delta_lambda_c = f_u/(H_u + a_u.transpose()*D_e_u*b_u);

        // Apply consistent correction to stress state and state variables.
        Voigt Delta_sigma_prime_c = -delta_lambda_c*D_e_u*b_u;
        State Delta_state_c = -delta_lambda_c*B_s_u;
        sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
        state_c = state_u + Delta_state_c;

        // Check corrected yield surface function value.
        f_c = compute_f(sigma_prime_c, state_c);
        if (std::abs(f_c) > std::abs(f_u)) {
            // Apply normal correction instead.
            double delta_lambda_c = f_u/(a_u.transpose()*a_u);
            
            // Update stress and state variables using correction.
            Voigt Delta_sigma_prime_c = -delta_lambda_c*a_u;
            sigma_prime_c = sigma_prime_u + to_cauchy(Delta_sigma_prime_c);
            state_c = state_u; /* i.e. no correction to state variables. */

            // Check yield surface function value again.
            f_c = compute_f(sigma_prime_c, state_c);
        }

        // Update uncorrected values for next iteration.
        f_u = f_c;
        sigma_prime_u = sigma_prime_c;
        state_u = state_c;
    }
    PLOG_DEBUG << "Yield surface drift correction converged in " << ITS_YSC << " iterations.";
}

State Elastoplastic::compute_elastic_state_variable_increment(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde_e) {
    State Delta_state = State::Zero(state.size());
    PLOG_VERBOSE << "In compute_elastic_state_variable_increment().";
    return Delta_state;
}

void Elastoplastic::compute_elastic_increment(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde_e, Cauchy &Delta_sigma_prime_e, State &Delta_state_e) {
    // Update state variables, if required (method does nothing if not overriden by model implementation).
    Delta_state_e = compute_elastic_state_variable_increment(sigma_prime, state, Delta_epsilon_tilde_e);

    // Elastic constitutive matrix.
    Constitutive D_e = compute_D_e(sigma_prime, to_cauchy(Delta_epsilon_tilde_e));

    // Elastic stress increment.
    Voigt Delta_sigma_prime_tilde_e = compute_elastic_stress_increment(D_e, Delta_epsilon_tilde_e);
    Delta_sigma_prime_e =  to_cauchy(Delta_sigma_prime_tilde_e);
}

void Elastoplastic::compute_plastic_increment(
    Cauchy sigma_prime, 
    State state, 
    Voigt Delta_epsilon_tilde_p_dT, 
    Voigt &Delta_sigma_prime, 
    State &delta_state
) {   
    // Calculate elastic constitutive matrix using tangent moduli and elastic stress increment.
    Constitutive D_e = compute_D_e(sigma_prime, Cauchy::Zero());

    // Derivatives of yield surface and plastic potential function.
    Derivatives derivatives = compute_derivatives(sigma_prime, state);
    Cauchy df_dsigma_prime = derivatives.df_dsigma_prime;
    Cauchy dg_dsigma_prime = derivatives.dg_dsigma_prime;
    Voigt a = to_voigt(df_dsigma_prime);
    Voigt b = to_voigt(dg_dsigma_prime);
    HardeningModuli H_s = derivatives.H_s;
    StateFactors B_s = derivatives.B_s;

    // Compute elastoplastic constitutive matrix.
    double H = H_s.sum();
    Constitutive D_ep = D_e-(D_e*b*a.transpose()*D_e)/(H + a.transpose()*D_e*b);

    // Compute elastoplastic multiplier. 
    Voigt Delta_sigma_prime_e = D_e*Delta_epsilon_tilde_p_dT;
    double delta_lambda = (double)(a.transpose()*Delta_sigma_prime_e)/(double)(H + a.transpose()*D_e*b);
    PLOG_DEBUG << "Plastic multiplier, delta_lambda = " << delta_lambda;

    // Update stress and state variable.
    Delta_sigma_prime = D_ep*Delta_epsilon_tilde_p_dT;
    delta_state = -delta_lambda*B_s;
}