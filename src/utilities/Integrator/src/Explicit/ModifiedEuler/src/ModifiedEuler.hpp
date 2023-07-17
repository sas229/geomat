#ifndef MODIFIEDEULER_H
#define MODIFIEDEULER_H

#include <plog/Log.h>
#include <functional>
#include <string>
#include "Types.hpp"


class ModifiedEuler {

    public:

        ModifiedEuler(Settings *settings, ModelFunctions *mf);

        ~ModifiedEuler() {}

        void solve(
            Cauchy &sigma_prime_ep, 
            State &state_ep, 
            Voigt Delta_epsilon_tilde_p 
        );

    private:
        double compute_new_substep_size(
            bool accepted, 
            double dT, 
            double T, 
            double R_n
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
            State &state_c
        );

        void compute_normal_yield_surface_correction(
            Cauchy sigma_prime_u, 
            State state_u, 
            double f_u, 
            Voigt a_u, 
            Cauchy &sigma_prime_c, 
            State &state_c
        );

        double compute_error_estimate(
            Cauchy sigma_prime_ini, 
            State state_ini, 
            Voigt Delta_sigma_prime_1, 
            Voigt Delta_sigma_prime_2, 
            State delta_state_1, 
            State delta_state_2
        );

        Settings *settings;

        ModelFunctions *mf;

};

#endif
