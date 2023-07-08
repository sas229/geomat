#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <plog/Log.h>
#include <functional>
#include <string>
#include "Types.hpp"


namespace Intersection {

/**
 * @brief Yield function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam state_f State variables.
 */
typedef std::function<double(Cauchy sigma_prime_f, State state_f)> YieldFunction;

/**
 * @brief Trial stress function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam delta_epsilon_tilde Strain increment in Voigt form.
 */
typedef std::function<Cauchy(Cauchy sigma_prime_f, Voigt delta_epsilon_tilde)> TrialFunction;

/**
 * @brief Constitutive matrix function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam Delta_epsilon Strain increment in Cuachy form.
 */
typedef std::function<Constitutive(Cauchy sigma_prime_f, Cauchy Delta_epsilon)> ConstitutiveMatrixFunction;

/**
 * @brief Derivative computation function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam state_f State variables.
 * @tparam df_dsigma_prime_f Derivatives of the yield function with respect to effective stress tensor in Cauchy form.
 * @tparam a_f Derivatives of the yield function with respect to effective stress tensor in Voigt form.
 * @tparam dg_dsigma_prime_f Derivatives of the plastic potential function with respect to state variables in Cauchy form.
 * @tparam b_f Derivatives of the plastic potential function with respect to state variables in Voigt form.
 * @tparam H_f Hardening modulus.
 */
typedef std::function<void(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H)> DerivativeFunction;


/**
 * @brief Function to compute the fraction \f$ \alpha \f$ of the step that is elastic.
 * 
 * The approach uses the "Pegasus" algorithm of Dowell and Jarret (1972) via @ref pegasus_regula_falsi. 
 * Due to the specific nature of the problem the algorithm is preconditioned for various cases to ensure
 * unconditional convergence. Taking the effective stress state and state variables at the start of the 
 * current increment, the current yield function value is computed as:
 * 
 * \f[ f_{0} = f\left(\sigma_{0}, H_{0}\right)\f]
 *  
 * The trial stress state is then computed using an elastic predictor and the full 
 * strain increment (i.e. \f$ \alpha = 1.0 \f$):
 * 
 * \f[ \sigma_{1} = \sigma_{0} + \mathbf{D}_{e} \alpha \Delta \epsilon = \sigma_{0} + \mathbf{D}_{e} \Delta \epsilon \f]
 * 
 * The value of the yield function at this trial stress is then computed:
 * 
 * \f[ f_{1} = f\left(\sigma_{1}, H_{0}\right)\f]
 *
 * If \f$ f_{1} < \text{FTOL} \f$, where FTOL is a user-defined tolerance typically taken as \f$ 1 \cdot 10^{-8} \f$,
 * then the increment is fully elastic and \f$ \alpha = 1.0 \f$ is returned. 
 * 
 * Otherwise, if \f$ \lvert f_{0} \rvert \leq \text{FTOL} \f$ and \f$ f_{1} \gt \text{FTOL} \f$ then the increment is partially plastic.
 * In this instance, a check is made for unloading-reloading via @ref check_unload_reload. If this check is true, then the
 * unloading-reloading occurs and the "Pegasus" algorithm is used to compute \f$ \alpha \f$ with bounds on \f$ \alpha \f$ 
 * of \f$ 0.0 \f$ and \f$ 1.0 \f$ via @ref pegasus_regula_falsi. Otherwise, if unloading-reloading does not occur closer bounds on 
 * \f$ \alpha \f$ are computed via @ref compute_alpha_bounds prior to calling via @ref pegasus_regula_falsi. If the algorithm fails 
 * to converge within a user-defined MAXITS_YSI iterations, typically taken as 10, then the algorithm raises a false assertion and
 * accompanying log message. The technique is model agnostic and can be used with any model because of the use of various function
 * bindings rather than hard-coded functions.
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
 * @brief Function to determine if an increment is an unload-reload increment.
 * 
 * An elastic to plastic transition can occur if a stress state initially lying on the yield surface is 
 * subjected to unloading prior to reloading (Sloan et al. 2001). This situation occurs when the angle
 * \f$ \theta \f$ between the yield surface gradient and the tangential elastic stress increment
 * \f$ \Delta \sigma_{e} \f$ is greater than \f$ 90^{\circ} \f$, i.e. when:
 * 
 * \f[ \cos \theta = \frac{a_c^{T} \Delta \sigma_{e}}{\left\|a_c\right\|_2\left\|\Delta \sigma_{e}\right\|_2} < -LTOL \f]
 *
 * where \f$ a_{c} \f$ represents the derivatives of the yield function \f$ f \f$ with respect to the 
 * current effective stress state in Voigt form, \f$ \Delta \sigma_{e} \f$ is the increment in the effective stress 
 * state assuming the increment is fully elastic and LTOL is a user-defined tolerance typically taken as
 * \f$ 1 \cdot 10^{-6} \f$. Returns true if the above criteria is met and false otherwise.
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
 * @brief Function to compute bounds for alpha for elastoplastic unloading-reloading increment.
 * 
 * The approach used to break the current strain increment into NSUB subincrements (typically 10)
 * in order to obtain a closer estimate of the range of \f$ \alpha \f$ values that bracket the elastic
 * to plastic transition. The subincrements are computed as:
 * 
 * \f[ \alpha_{n} = 1.0/\text{NSUB} \f]
 * 
 * giving the effective stress state after application of the subincrement as:
 * 
 * \f[ \sigma_{n} = \sigma_{0} + D_{e} \alpha_{n} \Delta \epsilon \f]
 * 
 * and correspondingly the yield surface function value as:
 * 
 * \f[ f_{n} = f\left( \sigma_{n}, H_{n} \right) \f]
 * 
 * If \f$ f_{n} \gt \text{FTOL} \f$ then the current subincrement contains the transition and bounds on 
 * \f$ \alpha \f$ of \f$ \alpha_{n} \pm 1.0/\text{NSUB}\f$ are returned. Otherwise the process continues 
 * with \f$ \alpha_{n} = \alpha_{n} + 1.0/\text{NSUB} \f$ recursively until \f$ f_{n} > \text{FTOL} \f$.
 *
 * @param[in] sigma_prime Current stress state.
 * @param[in] state Current state variables.
 * @param[in] Delta_epsilon_tilde Current strain increment.
 * @param[in] FTOL Yield surface tolerance.
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
    int NSUB,
    YieldFunction compute_f,
    ConstitutiveMatrixFunction compute_D_e,
    TrialFunction compute_trial_stress,
    double &alpha_0,
    double &alpha_1);

/**
 * @brief Pegasus regula falsi algorithm to find root of general non-linear system of equations.
 * 
 * The "Pegasus" regula falsi algorithm of Dowell and Jarrett (1972) is used to find the root of a general 
 * non-linear system of equations - in this case the elastic to plastic transition for an elastoplastic 
 * constitutive model. The algorithm begins with two initial guesses for the root, \f$ \alpha_{0} \f$ and
 * \f$ \alpha_{1} \f$, which are used to compute the corresponding yield function values \f$ f_{0} \f$ and
 * \f$ f_{1} \f$. The algorithm then proceeds to compute a new estimate for the root, \f$ \alpha_{n} \f$ via:
 * 
 * \f[ \alpha_{n} = \alpha_{1} - f_{1} \left[ \frac{\left( \alpha_{1} - \alpha_{0} \right)}{\left( f_{1} - f_{0} \right)} \right] \f]
 * 
 * Using the new estimate for the root, the corresponding yield function value \f$ f_{n} \f$ is computed. If the
 * product of \f$ f_{n} \text{ and } f_{1} \f$ is less than zero then:
 * 
 * \f[ \alpha_{1} = \alpha_{0} \f]
 * \f[ f_{1} = f_{0} \f]
 * 
 * Otherwise:
 * 
 * \f[ f_{1} = \frac{f_{1} f_{0}}{\left( f_{0} + f_{n} \right)} \f]
 * 
 * In either case:
 * 
 * \f[ \alpha_{0} = \alpha_{n} \text{ and } f_{0} = f_{n} \f]
 * 
 * Updated estimates for \f$ \alpha \f$ are computed until the yield function value is less than the user-defined
 * tolerance FTOL or until the number of iterations equals the user-defined maximum number of iterations MAXITS_YSI.
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
