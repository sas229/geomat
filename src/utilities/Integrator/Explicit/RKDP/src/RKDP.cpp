#include "RKDP.hpp"

RKDP::RKDP(Settings *settings, ModelFunctions *mf) {
    this->settings = settings;
    this->mf = mf;

    // Set order of method.
    this->order = 5.0;
}

void RKDP::compute_initial_estimate(void) {
    // Initial state.
    Cauchy sigma_prime_tilde_1 = sigma_prime_ep;
    State state_tilde_1 = state_ep;
    Voigt Delta_sigma_prime_1;
    State Delta_state_1;
    mf->compute_plastic_increment(sigma_prime_tilde_1, state_tilde_1, Delta_epsilon_tilde_dT, Delta_sigma_prime_1, Delta_state_1);
    PLOG_DEBUG << "State variable increment 1, Delta_state_1 = " << Delta_state_1;
    
    // 1st order estimate.
    Cauchy sigma_prime_tilde_2 = sigma_prime_tilde_1 + to_cauchy((1.0/5.0)*Delta_sigma_prime_1);
    State state_tilde_2 = state_tilde_1 + (1.0/5.0)*Delta_state_1;
    Voigt Delta_sigma_prime_2;
    State Delta_state_2;
    mf->compute_plastic_increment(sigma_prime_tilde_2, state_tilde_2, Delta_epsilon_tilde_dT, Delta_sigma_prime_2, Delta_state_2);
    PLOG_DEBUG << "State variable increment 2, Delta_state_2 = " << Delta_state_2;
    
    // 2nd order estimate.
    Cauchy sigma_prime_tilde_3 = sigma_prime_tilde_1 + to_cauchy((3.0/40.0)*Delta_sigma_prime_1 + (9.0/40.0)*Delta_sigma_prime_2);
    State state_tilde_3 = state_tilde_1 + (3.0/40.0)*Delta_state_1 + (9.0/40.0)*Delta_state_2;
    Voigt Delta_sigma_prime_3;
    State Delta_state_3;
    mf->compute_plastic_increment(sigma_prime_tilde_3, state_tilde_3, Delta_epsilon_tilde_dT, Delta_sigma_prime_3, Delta_state_3);
    PLOG_DEBUG << "State variable increment 3, Delta_state_3 = " << Delta_state_3;
    
    // 3rd order estimate.
    Cauchy sigma_prime_tilde_4 = sigma_prime_tilde_1 + to_cauchy((3.0/10.0)*Delta_sigma_prime_1 - (9.0/10.0)*Delta_sigma_prime_2 + (6.0/5.0)*Delta_sigma_prime_3);
    State state_tilde_4 = state_tilde_1 + (3.0/10.0)*Delta_state_1 - (9.0/10.0)*Delta_state_2 + (6.0/5.0)*Delta_state_3;
    Voigt Delta_sigma_prime_4;
    State Delta_state_4;
    mf->compute_plastic_increment(sigma_prime_tilde_4, state_tilde_4, Delta_epsilon_tilde_dT, Delta_sigma_prime_4, Delta_state_4);
    PLOG_DEBUG << "State variable increment 4, Delta_state_4 = " << Delta_state_4;
    
    // 4th order estimate.
    Cauchy sigma_prime_tilde_5 = sigma_prime_tilde_1 + to_cauchy((226.0/729.0)*Delta_sigma_prime_1 - (25.0/27.0)*Delta_sigma_prime_2 + (880.0/729.0)*Delta_sigma_prime_3 + (55.0/729.0)*Delta_sigma_prime_4);
    State state_tilde_5 = state_tilde_1 + (226.0/729.0)*Delta_state_1 - (25.0/27.0)*Delta_state_2 + (880.0/729.0)*Delta_state_3 + (55.0/729.0)*Delta_state_4;
    Voigt Delta_sigma_prime_5;
    State Delta_state_5;
    mf->compute_plastic_increment(sigma_prime_tilde_5, state_tilde_5, Delta_epsilon_tilde_dT, Delta_sigma_prime_5, Delta_state_5);
    PLOG_DEBUG << "State variable increment 5, Delta_state_5 = " << Delta_state_5;
    
    // 5th order estimate.
    Cauchy sigma_prime_tilde_6 = sigma_prime_tilde_1 + to_cauchy(- (181.0/270.0)*Delta_sigma_prime_1 + (5.0/2.0)*Delta_sigma_prime_2 - (226.0/297.0)*Delta_sigma_prime_3 - (91.0/27.0)*Delta_sigma_prime_4 + (189.0/55.0)*Delta_sigma_prime_5);
    State state_tilde_6 = state_tilde_1 - (181.0/270.0)*Delta_state_1 + (5.0/2.0)*Delta_state_2 - (226.0/297.0)*Delta_state_3 - (91.0/27.0)*Delta_state_4 + (189.0/55.0)*Delta_state_5;
    Voigt Delta_sigma_prime_6;
    State Delta_state_6;
    mf->compute_plastic_increment(sigma_prime_tilde_6, state_tilde_6, Delta_epsilon_tilde_dT, Delta_sigma_prime_6, Delta_state_6);
    PLOG_DEBUG << "State variable increment 6, Delta_state_6 = " << Delta_state_6;

    // Calculate RKDP estimate of effective stresses and state variables.
    sigma_prime_ini = sigma_prime_ep + to_cauchy(
        (19.0/216.0)*Delta_sigma_prime_1 + 
        (1000.0/2079.0)*Delta_sigma_prime_3 - 
        (125.0/216.0)*Delta_sigma_prime_4 + 
        (81.0/88.0)*Delta_sigma_prime_5 + 
        (5.0/56.0)*Delta_sigma_prime_6
    );
    state_ini = state_ep +
        (19.0/216.0)*Delta_state_1 + 
        (1000.0/2079.0)*Delta_state_3 - 
        (125.0/216.0)*Delta_state_4 + 
        (81.0/88.0)*Delta_state_5 + 
        (5.0/56.0)*Delta_state_6;
    PLOG_DEBUG << "Initial stress estimate, sigma_prime_ini_tilde = " << to_voigt(sigma_prime_ini);
    PLOG_DEBUG << "State variable estimate, state_ini = " << state_ini;

    // Compute error estimate.
    int size_error = 1 + state_ini.size();
    Eigen::VectorXd error(size_error);
    error[0] = (to_cauchy(
        (11.0/360.0)*Delta_sigma_prime_1 - 
        (10.0/63.0)*Delta_sigma_prime_3 +
        (55.0/72.0)*Delta_sigma_prime_4 -
        (27.0/40.0)*Delta_sigma_prime_5 +
        (11.0/280.0)*Delta_sigma_prime_6
    )).norm()/(sigma_prime_ini.norm());
    for (int i=1; i<size_error; ++i) {
        error[i] = std::abs(
            (11.0/360.0)*Delta_state_1[i-1] - 
            (10.0/63.0)*Delta_state_3[i-1] +
            (55.0/72.0)*Delta_state_4[i-1] -
            (27.0/40.0)*Delta_state_5[i-1] +
            (11.0/280.0)*Delta_state_6[i-1]
        )/state_ini[i-1];
    }
    PLOG_DEBUG << "Initial estimate error vector: " << error;
    R_n = error.maxCoeff();

    // Check error estimate against machine tolerance.
    R_n = std::max(R_n, settings->EPS);
    PLOG_DEBUG << "R_n = " << R_n;

    // Check acceptance of estimate.
    if (R_n < settings->STOL) {
        accepted = true;
    } else {
        accepted = false;
    }
}