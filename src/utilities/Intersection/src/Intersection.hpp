#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <plog/Log.h>
#include <string>
#include <functional>
#include "Types.hpp"

namespace Intersection {
    
    /**
     * @brief Method to compute the elastic fraction of the current strain increment following Sloan et al. (2001).
     * 
     * @param[in] sigma_prime Current stress state.
     * @param[in] state Current state variables.
     * @return double
     */
    double compute_alpha(
        Cauchy sigma_prime, 
        State state, 
        Voigt Delta_epsilon_tilde,
        YieldFunction compute_f,
        TrialFunction compute_trial_stress,
        ConstitutiveMatrixFunction compute_D_e,
        DerivativeFunction compute_derivatives
        ); 

    // /**
    //  * @brief Method to determine if an increment is an unload-reload plastic increment.
    //  * 
    //  * @param[in] sigma_prime Current stress state.
    //  * @param[in] state Current state variables.
    //  * @return true
    //  * @return false 
    //  */
    // bool check_unload_reload(Cauchy sigma_prime, State state);        

    // /**
    //  * @brief Method to compute bounds for alpha for elastoplastic unloading-reloading increment.
    //  * 
    //  * @param[in,out] alpha_0 Lower bound for alpha.
    //  * @param[in,out] alpha_1 Upper bound for alpha.
    //  */
    // void compute_alpha_bounds(double &alpha_0, double &alpha_1);

    // /**
    //  * @brief Pegasus regula falsi algorithm to find root of general non-linear system of equations.
    //  * 
    //  * @param[in] sigma_prime Current stress state.
    //  * @param[in] state Current state variables.
    //  * @param[in] alpha_0 Lower bound on alpha.
    //  * @param[in] alpha_1 Upper bound on alpha.
    //  * @param[in] f_0 Initial value of objective function with lower bound alpha.
    //  * @param[in] f_1 Initial value of objective function with upper bound alpha.
    //  * @return double
    //  */
    // double pegasus_regula_falsi(Cauchy sigma_prime, State state, double alpha_0, double alpha_1, double f_0, double f_1);

}

#endif