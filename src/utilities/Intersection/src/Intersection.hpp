#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <plog/Log.h>
#include "Tensor.hpp"
#include "Types.hpp"

class Intersection {

    public:

        /**
         * @brief Intersection class constructor.
         * 
         * @param settings Stress integration settings.
         * @param mf Model specific function bindings. 
         */
        Intersection(Settings *settings, ModelFunctions *mf);

        /** 
         * @brief Intersection class destructor. 
         */
        virtual ~Intersection() {}

        /**
         * @brief Method to compute the fraction \f$ \alpha \f$ of the step that is elastic.
         * 
         * The approach uses the "Pegasus" algorithm of Dowell and Jarret (1972) via @ref pegasus_regula_falsi. 
         * Due to the specific nature of the problem the algorithm is preconditioned for various cases to ensure
         * unconditional convergence. Taking the effective stress state and state variables at the start of the 
         * current increment, the current yield function value is computed as:
         * 
         * \f[ f_{0} = f\left(\sigma^{\prime}_{0}, H_{0}\right)\f]
         *  
         * where \f$ \sigma^{\prime}_{0} \f$ is the initial effective stress state and \f$ H_{0} \f$ 
         * is the initial vector of state variables. The trial stress state is then computed using an 
         * elastic predictor and the full strain increment (i.e. \f$ \alpha = 1.0 \f$):
         * 
         * \f[ \sigma^{\prime}_{1} = \sigma^{\prime}_{0} + \mathbf{D}_{e} \alpha \Delta \epsilon = \sigma^{\prime}_{0} + \mathbf{D}_{e} \Delta \epsilon \f]
         * 
         * The value of the yield function at this trial stress is then computed:
         * 
         * \f[ f_{1} = f\left(\sigma^{\prime}_{1}, H_{0}\right)\f]
         *
         * If \f$ f_{1} < \text{FTOL} \f$, which is a user-defined tolerance typically taken as \f$ 1 \cdot 10^{-8} \f$,
         * then the increment is fully elastic and \f$ \alpha = 1.0 \f$ is returned. 
         * 
         * Otherwise, if \f$ \lvert f_{0} \rvert \leq \text{FTOL} \f$ and \f$ f_{1} \gt \text{FTOL} \f$ then the increment is partially plastic.
         * In this instance, a check is made for unloading-reloading via @ref check_unload_reload. If this check is true, then
         * unloading-reloading occurs and the "Pegasus" algorithm is used to compute \f$ \alpha \f$ with bounds on \f$ \alpha \f$ 
         * of \f$ 0.0 \f$ and \f$ 1.0 \f$ via @ref pegasus_regula_falsi. Otherwise, if unloading-reloading does not occur closer bounds on 
         * \f$ \alpha \f$ are computed via @ref refine_alpha_bounds prior to calling via @ref pegasus_regula_falsi. If the algorithm fails 
         * to converge within a user-defined MAXITS_YSI iterations, typically taken as 10, then the algorithm raises a false assertion and
         * accompanying log message. The technique is model agnostic and can be used with any model because of the use of various function
         * bindings rather than hard-coded functions.
         * 
         * @param[in] sigma_prime Current stress state.
         * @param[in] state Current state variables.
         * @param[in] Delta_epsilon_tilde Current strain increment.
         * @return double
         */
        double solve(
            Cauchy sigma_prime,
            State state,
            Voigt Delta_epsilon_tilde
        );

    private:

        /**
         * @brief Method to determine if an increment is an unload-reload increment.
         * 
         * An elastic to plastic transition can occur if a stress state initially lying on the yield surface is 
         * subjected to unloading prior to reloading (Sloan et al. 2001). This situation occurs when the angle
         * \f$ \theta \f$ between the yield surface gradient and the tangential elastic stress increment
         * \f$ \Delta \sigma^{\prime}_{e} \f$ is greater than \f$ 90^{\circ} \f$, i.e. when:
         * 
         * \f[ \cos \theta = \frac{a_c^{T} \Delta \sigma^{\prime}_{e}}{\left\|a_c\right\|_2\left\|\Delta \sigma^{\prime}_{e}\right\|_2} < -LTOL \f]
         *
         * where \f$ a_{c} \f$ represents the derivatives of the yield function \f$ f \f$ with respect to the 
         * current effective stress state in Voigt form, \f$ \Delta \sigma^{\prime}_{e} \f$ is the increment in the effective stress 
         * state assuming the increment is fully elastic and LTOL is a user-defined tolerance typically taken as
         * \f$ 1 \cdot 10^{-6} \f$. Returns true if the above criteria is met and false otherwise.
         * 
         * @return true
         * @return false
         */
        bool check_unload_reload(void);

        /**
         * @brief Method to compute bounds for alpha for elastoplastic unloading-reloading increment.
         * 
         * The approach used to break the current strain increment into NSUB subincrements (typically 10)
         * in order to obtain a closer estimate of the range of \f$ \alpha \f$ values that bracket the elastic
         * to plastic transition. The subincrements are computed as:
         * 
         * \f[ \alpha_{n} = 1.0/\text{NSUB} \f]
         * 
         * giving the effective stress state after application of the subincrement as:
         * 
         * \f[ \sigma^{\prime}_{n} = \sigma^{\prime}_{0} + D_{e} \alpha_{n} \Delta \epsilon \f]
         * 
         * and correspondingly the yield surface function value as:
         * 
         * \f[ f_{n} = f\left( \sigma^{\prime}_{n}, H_{n} \right) \f]
         * 
         * where \f$ H_{n} \f$ represents the vector of state variables after application of the subincrement.
         * 
         * If \f$ f_{n} \gt \text{FTOL} \f$ then the current subincrement contains the transition and bounds on 
         * \f$ \alpha \f$ of \f$ \alpha_{n} \pm 1.0/\text{NSUB}\f$ are returned. Otherwise the process continues 
         * with \f$ \alpha_{n} = \alpha_{n} + 1.0/\text{NSUB} \f$ recursively until \f$ f_{n} > \text{FTOL} \f$.
         */
        void refine_alpha_bounds(void);

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
         * @return double
         */
        double pegasus_regula_falsi(void);

        /**
         * @brief Lower bound on alpha.
         */
        double alpha_0 = 0.0;

        /**
         * @brief Upper bound on alpha.
         */
        double alpha_1 = 1.0;

        /**
         * @brief Fraction of strain increment applied at point of intersection.
         */
        double alpha;

        /**
         * @brief Yield function value for initial state.
         */
        double f_0;

        /**
         * @brief Yield function value after application of strain increment.
         */
        double f_1;

        /**
         * @brief Initial effective stress state to solve for surface intersection.
         * 
         */
        Cauchy sigma_prime;

        /**
         * @brief Initial state variables to solve for surface intersection.
         */
        State state;

        /**
         * @brief Strain increment to solve for surface intersection.
         * 
         */
        Voigt Delta_epsilon_tilde;

        /**
         * @brief Stress integration settings.
         */
        Settings *settings;

        /**
         * @brief Model specfic function bindings.
         */
        ModelFunctions *mf;

};

#endif
