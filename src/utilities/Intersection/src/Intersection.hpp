#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <plog/Log.h>
#include <functional>
#include <string>
#include "Types.hpp"


namespace Intersection {

/**
 * @brief Yield function binding.
 */
typedef std::function<double(Cauchy sigma_prime_f, State state_f)> YieldFunction;

/**
 * @brief Trial stress function binding.
 */
typedef std::function<Cauchy(Cauchy sigma_prime_f, Voigt delta_epsilon_tilde)> TrialFunction;

/**
 * @brief Constitutive matrix function binding.
 */
typedef std::function<Constitutive(Cauchy sigma_prime_f, Cauchy Delta_epsilon)> ConstitutiveMatrixFunction;

/**
 * @brief Derivative computation function binding.
 */
typedef std::function<void(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H)> DerivativeFunction;


/**
 * @brief 
 *
 * @param[in] sigma_prime Current stress state.
 * @param[in] state Current state variables.
 * @param[in] Delta_epsilon_tilde Current strain increment.
 * @param[in] FTOL Yield surface tolerance.
 * @param[in] LTOL Unload-reload tolerance.
 * @param[in] MAXITS_YSI Maximum number of iterations for the yield surface intersection algorithm.
 * @param[in] NSUB Maximum number of subincrements for the yield surface intersection algorithm.
 * @param[in] compute_f Yield function binding.
 * @param[in] compute_trial_stress Trial stress function binding.
 * @param[in] compute_D_e Constitutive matrix function binding.
 * @param[in] compute_derivatives Derivative computation function binding.
 * @return double
 */
double compute_alpha(
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
    DerivativeFunction compute_derivatives);

/**
 * @brief Method to determine if an increment is an unload-reload plastic increment.
 *
 * @param[in] sigma_prime Current stress state.
 * @param[in] state Current state variables.
 * @param[in] Delta_epsilon_tilde Current strain increment.
 * @param[in] LTOL Unload-reload tolerance.
 * @param[in] compute_f Yield function binding.
 * @param[in] compute_D_e Constitutive matrix function binding.
 * @param[in] compute_derivatives Derivative computation function binding.
 * @return true
 * @return false
 */
bool check_unload_reload(
    Cauchy sigma_prime,
    State state,
    Voigt Delta_epsilon_tilde,
    double LTOL,
    YieldFunction compute_f,
    ConstitutiveMatrixFunction compute_D_e,
    DerivativeFunction compute_derivatives);

/**
 * @brief Method to compute bounds for alpha for elastoplastic unloading-reloading increment.
 *
 * @param[in] sigma_prime Current stress state.
 * @param[in] state Current state variables.
 * @param[in] Delta_epsilon_tilde Current strain increment.
 * @param[in] FTOL Yield surface tolerance.
 * @param[in] MAXITS_YSI Maximum number of iterations for the yield surface intersection algorithm.
 * @param[in] NSUB Maximum number of subincrements for the yield surface intersection algorithm.
 * @param[in] compute_f Yield function binding.
 * @param[in] compute_D_e Constitutive matrix function binding.
 * @param[in] compute_trial_stress Trial stress function binding.
 * @param[in,out] alpha_0 Lower bound for alpha.
 * @param[in,out] alpha_1 Upper bound for alpha.
 */
void compute_alpha_bounds(
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
    double &alpha_1);

/**
 * @brief Pegasus regula falsi algorithm to find root of general non-linear system of equations.
 *
 * @param[in] sigma_prime Current stress state.
 * @param[in] state Current state variables.
 * @param[in] Delta_epsilon_tilde Current strain increment.
 * @param[in] alpha_0 Lower bound on alpha.
 * @param[in] alpha_1 Upper bound on alpha.
 * @param[in] f_0 Initial value of objective function with lower bound alpha.
 * @param[in] f_1 Initial value of objective function with upper bound alpha.
 * @param[in] FTOL Yield surface tolerance.
 * @param[in] MAXITS_YSI Maximum number of iterations for the yield surface intersection algorithm.
 * @param[in] compute_f Yield function binding.
 * @param[in] compute_trial_stress Trial stress function binding.
 * @return double
 */
double pegasus_regula_falsi(
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
    TrialFunction compute_trial_stress);

}  // namespace Intersection

#endif
