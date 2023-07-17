#include "Elastoplastic.hpp"

Elastoplastic::Elastoplastic() : intersection(&settings, &mf), integrator(&settings, &mf) {
    // Define binds to model functions.
    using namespace std::placeholders;
    mf.compute_f = std::bind(&Elastoplastic::compute_f, this, _1, _2);
    mf.compute_trial_stress = std::bind(&Elastoplastic::compute_elastic_stress, this, _1, _2);
    mf.compute_D_e = std::bind(&Elastoplastic::compute_D_e, this, _1, _2);
    mf.compute_derivatives = std::bind(&Elastoplastic::compute_derivatives, this, _1, _2, _3, _4, _5, _6, _7);
    mf.compute_state_increment = std::bind(&Elastoplastic::compute_plastic_state_variable_increment, this, _1, _2, _3, _4);
    mf.compute_plastic_increment = std::bind(&Elastoplastic::compute_plastic_increment, this, _1, _2, _3, _4, _5);
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
        PLOG_INFO << "Plastic increment; alpha = " << alpha;
        Cauchy sigma_prime_e, sigma_prime_ep;
        State state_e, state_ep;
        Voigt Delta_epsilon_tilde_e = alpha*Delta_epsilon_tilde;
        if (alpha == 0.0) {
            // Fully plastic increment. Update stress and state variables.
            sigma_prime_e = sigma_prime;
            state_e = state;
        }
        if (alpha > 0.0) {
            // Update stress and state variables based on elastic portion of strain increment.
            sigma_prime_e = sigma_prime + to_cauchy(D_e*alpha*Delta_epsilon_tilde);
            state_e = compute_elastic_state_variable(Delta_epsilon_tilde_e);
        } 
        PLOG_DEBUG << "Elastic strain increment, Delta_epsilon_e = " << Delta_epsilon_tilde_e;
        PLOG_DEBUG << "Stress state after elastic increment, sigma_prime_e_tilde = " << to_voigt(sigma_prime_e);
        PLOG_DEBUG << "State variables after elastic increment, state_e = " << state_e;
        if (alpha < 1.0) {
            // Perform elastoplastic stress integration on plastic portion of strain increment.
            Voigt Delta_epsilon_tilde_p = (1.0-alpha)*Delta_epsilon_tilde;
            sigma_prime_ep = sigma_prime_e; 
            state_ep = state_e;
            PLOG_DEBUG << "Plastic strain increment, Delta_epsilon_tilde_p = " << Delta_epsilon_tilde_p;
            PLOG_DEBUG << "Initial elastoplastic stress state, sigma_prime_ep_tilde = " << to_voigt(sigma_prime_ep);
            PLOG_DEBUG << "Initial elastoplastic state variables, state_ep = " <<state_ep;
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

void Elastoplastic::compute_plastic_increment(
    Cauchy sigma_prime, 
    State state, 
    Voigt Delta_epsilon_tilde_p_dT, 
    Voigt &Delta_sigma_prime, 
    State &delta_state
    ) {   
    // Calculate elastic constitutive matrix using tangent moduli and elastic stress increment.
    Constitutive D_e = compute_D_e(sigma_prime, Cauchy::Zero());

    // Compute elastoplastic constitutive matrix and elastoplastic multiplier. 
    Cauchy df_dsigma_prime, dg_dsigma_prime;
    Voigt a, b;
    double H;
    compute_derivatives(sigma_prime, state, df_dsigma_prime, a, dg_dsigma_prime, b, H);
    Constitutive D_ep = D_e-(D_e*b*a.transpose()*D_e)/(H+a.transpose()*D_e*b);
    Voigt Delta_sigma_prime_e = D_e*Delta_epsilon_tilde_p_dT;
    double delta_lambda = (double)(a.transpose()*Delta_sigma_prime_e)/(double)(H + a.transpose()*D_e*b);
    PLOG_DEBUG << "Plastic multiplier, delta_lambda = " << delta_lambda;

    // Update stress and state variable increment by reference.
    Delta_sigma_prime = D_ep*Delta_epsilon_tilde_p_dT;
    delta_state = compute_plastic_state_variable_increment(delta_lambda, df_dsigma_prime, H, Delta_epsilon_tilde_p_dT);
}