#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <plog/Log.h>
#include "Types.hpp"

class Integrator {

    public:

        /**
         * @brief Integrator constructor.
         * 
         * @param settings Stress integration settings.
         * @param mf Model specific function bindings. 
         */
        Integrator(Settings settings, ModelFunctions mf);

        /** 
         * @brief Integrator class destructor. 
         */
        virtual ~Integrator() {}

        void solve(
            Cauchy &sigma_prime_ep, 
            State &state_ep, 
            Voigt Delta_epsilon_tilde_p, 
        );

    private:

        /**
         * @brief Method to compute the new substep size.
         * 
         * @return double 
         */
        double compute_new_substep_size(void);

        /**
         * @brief Method to compute the yield surface drift correction.
         */
        void compute_yield_surface_correction(void);
        
        /**
         * @brief Method to compute the normal yield surface drift correction.
         */
        void compute_normal_yield_surface_correction(void);

        /**
         * @brief Method to compute a plastic increment.
         */
        void compute_plastic_increment(void);

        /**
         * @brief Method to compute the error estimate.
         *
         * @return double 
         */
        double compute_error_estimate(void);

        /**
         * @brief Stress integration settings.
         */
        Settings settings;

        /**
         * @brief Model specfic function bindings.
         */
        ModelFunctions mf;
        
}

#endif
