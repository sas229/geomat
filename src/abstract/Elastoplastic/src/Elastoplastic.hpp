#ifndef ELASTOPLASTIC_H
#define ELASTOPLASTIC_H

#include "Tensor.hpp"
#include "Types.hpp"
#include "Elastic.hpp"
#include "Intersection.hpp"
#include "Integrator.hpp"

/**
 * @brief Elastoplastic base class. Inherits Elastic class.
 */
class Elastoplastic : public Elastic {

    public: 

        /** 
         * @brief Elastoplastic model constructor.
         */
        explicit Elastoplastic();

        /**
         *  @brief Elastoplastic model destructor. 
         */
        virtual ~Elastoplastic() {}

        /**
         *  @brief Solve current strain increment using refined explicit approach of Sloan et al. (2001). 
         */
        void solve(void);

        /**
         * @brief Pure virtual method to compute the yield surface value given the parameters, current state variables and stress state.
         * 
         * @param sigma_prime Effective stress tensor. 
         * @param state State variables.
         * @return double
         * 
         * @note Must be overriden by model implementations.
         */
        virtual double compute_f(Cauchy sigma_prime, State state) = 0;

        /** 
         * @brief Pure virtual method to compute the derivatives for the constitutive model implemented.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[in] state State variables.
         * @param[in,out] df_dsigma_prime Derivatives of the yield function with respect to the stress state.
         * @param[in,out] dg_dsigma_prime Derivatives of the plastic potential function with respect to the stress state.
         * @param[in,out] H_s Hardening moduli vector.
         * @param[in,out] B_s State factor vector.
         * 
         * @note Must be overriden by model implementations.
         */
        virtual void compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Cauchy &dg_dsigma_prime, HardeningModuli &H_s, StateFactors &B_s) = 0;

        /**
         * @brief Compute plastic stress integration increment.
         * 
         * @param[in] sigma_prime Effective stress state.
         * @param[in] state Vector of state variables.
         * @param[in] Delta_epsilon_tilde_p_dT Strain increment.
         * @param[in,out] Delta_sigma_prime Increment in stress state.
         * @param[in,out] delta_state Increment in state variables.
         */
        void compute_plastic_increment(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde_p_dT, Voigt &Delta_sigma_prime, State &delta_state);
        
        /**
         * @brief Object containing model specific function bindings.
         */
        ModelFunctions mf;

        /** 
         * @brief Stress integration settings.
        */
        Settings settings;

        /**
         * @brief Intersection class instance.
         */
        Intersection intersection;

        /**
         * @brief Integrator class instance.
         */
        Integrator integrator;

};

#endif
