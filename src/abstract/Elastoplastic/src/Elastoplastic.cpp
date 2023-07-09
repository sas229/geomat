#include "Elastoplastic.hpp"

void Elastoplastic::solve(void) {
    if (!solved) {
        // Get the current state variables.
        State state = get_state_variables();

        // Define binds to class methods.
        auto compute_f_func = std::bind(&Elastoplastic::compute_f, this, std::placeholders::_1, std::placeholders::_2);
        auto compute_trial_stress_func = std::bind(&Elastoplastic::compute_elastic_stress, this, std::placeholders::_1, std::placeholders::_2);
        auto compute_D_e_func = std::bind(&Elastoplastic::compute_D_e, this, std::placeholders::_1, std::placeholders::_2);
        auto compute_derivatives_func = std::bind(&Elastoplastic::compute_derivatives, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7);
        auto compute_state_func = std::bind(&Elastoplastic::compute_plastic_state_variable_increment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

        // Compute alpha using bound functions.
        double alpha = Intersection::compute_alpha(sigma_prime, state, Delta_epsilon_tilde, FTOL, LTOL, MAXITS_YSI, NSUB, compute_f_func, compute_trial_stress_func, compute_D_e_func, compute_derivatives_func);
        
        Voigt Delta_epsilon_tilde_e = alpha*Delta_epsilon_tilde;
        PLOG_DEBUG << "Strain increment, Delta_epsilon_tilde = " << Delta_epsilon_tilde;
        PLOG_DEBUG << "Initial stress state, sigma_prime_tilde = " << to_voigt(sigma_prime);
        PLOG_DEBUG << "Initial state variables, state = " << state;
        PLOG_INFO << "Plastic increment; alpha = " << alpha;
        Cauchy sigma_prime_e, sigma_prime_ep;
        State state_e, state_ep;
        if (alpha == 0.0) {
            // Fully plastic increment. Update stress and state variables.
            sigma_prime_e = sigma_prime;
            state_e = get_state_variables();
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
            if (method == "Sloan") {
                // Perform elastoplastic stress integration using Sloan's method.
                Sloan::solve(sigma_prime_ep, state_ep, Delta_epsilon_tilde_p, FTOL, STOL, DT_MIN, EPS, MAXITS_YSC, compute_f_func, compute_D_e_func, compute_derivatives_func, compute_state_func);
            } else {
                PLOG_FATAL << "Invalid elastoplastic integration method.";
                assert(false);
            }
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