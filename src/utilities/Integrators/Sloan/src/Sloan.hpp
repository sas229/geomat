#ifndef SLOAN_H
#define SLOAN_H

#include <plog/Log.h>
#include <functional>
#include <string>
#include "Types.hpp"


namespace Sloan {

    void solve(
        Cauchy &sigma_prime_ep, 
        State &state_ep, 
        Voigt Delta_epsilon_tilde_p, 
        double FTOL,
        double STOL,
        double DT_MIN,
        double EPS,
        int MAXITS_YSC,
        YieldFunction compute_f,
        ConstitutiveMatrixFunction compute_D_e,
        DerivativeFunction compute_derivatives,
        StateIncrementFunction compute_plastic_state_variable_increment
    );

    double compute_new_substep_size(
        bool accepted, 
        double dT, 
        double T, 
        double R_n, 
        double DT_MIN, 
        double STOL
    );

    void compute_yield_surface_correction(
        Cauchy sigma_prime_u, 
        State state_u, 
        double f_u, 
        double H_u, 
        Voigt a_u, 
        Voigt b_u, 
        Constitutive D_e_u, 
        Cauchy &sigma_prime_c, 
        State &state_c, 
        StateIncrementFunction compute_plastic_state_variable_increment
    );

    void compute_normal_yield_surface_correction(
        Cauchy sigma_prime_u, 
        State state_u, 
        double f_u, 
        Voigt a_u, 
        Cauchy &sigma_prime_c, 
        State &state_c
    );

    void compute_plastic_increment(
        Cauchy sigma_prime, 
        State state, 
        Voigt Delta_epsilon_tilde_p_dT, 
        Voigt &Delta_sigma_prime, 
        State &delta_state, 
        ConstitutiveMatrixFunction compute_D_e, 
        DerivativeFunction compute_derivatives, 
        StateIncrementFunction compute_plastic_state_variable_increment
    );

    double compute_error_estimate(
        Cauchy sigma_prime_ini, 
        State state_ini, 
        Voigt Delta_sigma_prime_1, 
        Voigt Delta_sigma_prime_2, 
        State delta_state_1, 
        State delta_state_2, 
        double EPS
    );

}  // namespace Sloan

#endif
