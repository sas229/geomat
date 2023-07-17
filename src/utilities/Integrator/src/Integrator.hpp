#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <plog/Log.h>
#include "Types.hpp"
#include "ModifiedEuler.hpp"
#include "RKDP.hpp"

class Integrator {

    public:

        /**
         * @brief Integrator constructor.
         * 
         * @param settings Stress integration settings.
         * @param mf Model specific function bindings. 
         */
        Integrator(Settings *settings, ModelFunctions *mf);

        /** 
         * @brief Integrator class destructor. 
         */
        virtual ~Integrator() {}

        void solve(
            Cauchy &sigma_prime, 
            State &state, 
            Voigt Delta_epsilon_tilde
        );

    private:

        /**
         * @brief Stress integration settings.
         */
        Settings *settings;

        /**
         * @brief Model specfic function bindings.
         */
        ModelFunctions *mf;

        /**
         * @brief Modified Euler stress integration solver instance.
         */
        ModifiedEuler me;

        /**
         * @brief Runge-Kutta-Dormand-Prince integration solver instance.
         * 
         */
        RKDP rkdp;
        
};

#endif
