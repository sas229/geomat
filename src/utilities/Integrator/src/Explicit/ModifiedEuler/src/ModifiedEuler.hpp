#ifndef MODIFIEDEULER_H
#define MODIFIEDEULER_H

#include <plog/Log.h>
#include <vector>
#include <string>
#include "Types.hpp"


class ModifiedEuler {

    public:

        ModifiedEuler(Settings *settings, ModelFunctions *mf);

        ~ModifiedEuler() {}

        void solve(
            Cauchy &sigma_prime, 
            State &state, 
            Voigt Delta_epsilon_tilde
        );

    private:

        void compute_initial_estimate(void);

        void compute_yield_surface_correction(void);

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
        const double order = 2.0;

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
