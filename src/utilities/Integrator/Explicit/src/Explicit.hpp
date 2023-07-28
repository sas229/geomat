#ifndef EXPLICIT_H
#define EXPLICIT_H

#include <plog/Log.h>
#include "Tensor.hpp"
#include "Types.hpp"

class Explicit {

    public:

        /** @brief Virtual destructor for Explicit integrator base class. */
        virtual ~Explicit() {}

        /**
         * @brief Method to solve the current incremement.
         * 
         * @param[in,out] sigma_prime Effective stress state in Cauchy form.
         * @param[in,out] state State variables.
         * @param[in] Delta_epsilon_tilde Strain increment in Voigt form.
         */
        void solve(
            Cauchy &sigma_prime, 
            State &state, 
            Voigt Delta_epsilon_tilde
        );

    protected:

        /**
         * @brief Pure virtual method to compute an initial estimate for the current increment.
         *
         * @note Must be overriden by the solver classes that inherit from it (e.g. ModifiedEuler).
         */
        virtual void compute_initial_estimate(void) = 0;

        /**
         * @brief Method to compute the new substep size for the next increment.
         * 
         * @return double 
         */
        double compute_new_substep_size(void);

        /**
         * @brief Stress integration settings.
         */
        Settings *settings;

        /**
         * @brief Model specfic function bindings.
         */
        ModelFunctions *mf;

        /**
         * @brief Order of integration method.
         */
        double order;

        /**
         * @brief Effective stress state in Cauchy form.
         */
        Cauchy sigma_prime_ep;

        /**
         * @brief State variables.
         */
        State state_ep;

        /**
         * @brief Strain increment for pseudo-time dT.
         */
        Voigt Delta_epsilon_tilde_dT;

        /**
         * @brief Initial effective stress state estimate.
         */
        Cauchy sigma_prime_ini;

        /**
         * @brief Initial state variable estimate.
         */
        State state_ini;

        /**
         * @brief Uncorrected effective stress state.
         */
        Cauchy sigma_prime_u;

        /**
         * @brief Uncorrected state variables.
         */
        State state_u;

        /**
         * @brief Yield function value for uncorrected state.
         */
        double f_u;

        /**
         * @brief Corrected effective stress state.
         */
        Cauchy sigma_prime_c;

        /**
         * @brief Corrected state variables.
         */
        State state_c;

        /**
         * @brief Yield function value for corrected state.
         */
        double f_c;

        /**
         * @brief Boolean to indicate if current increment is accepted.
         */
        bool accepted;

        /**
         * @brief Error estimate.
         */
        double R_n;

        /**
         * @brief Current pseudo-time increment.
         */
        double dT;

        /**
         * @brief Current pseudo-time in the range of 0.0-1.0.
         */
        double T;
};

#endif
