#ifndef ELASTOPLASTIC_H
#define ELASTOPLASTIC_H

#include "Elastic.hpp"
#include "Intersection.hpp"
#include "Integrator.hpp"
#include "Tensor.hpp"
#include "Types.hpp"

/**
 * @brief Elastoplastic base class. Inherits Elastic class.
 */
class Elastoplastic : public Elastic {

    public: 

        /** 
         * @brief Elastoplastic model constructor.
         */
        explicit Elastoplastic(std::string log_severity="none");

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
         * @return Derivatives
         * 
         * @note Must be overriden by model implementations.
         */
        virtual Derivatives compute_derivatives(Cauchy sigma_prime, State state) = 0;

        /**
         * @brief Method to check for yield surface drift and apply appropriate corrections.
         * 
         * @param sigma_prime_u Uncorrected effective stress tensor.
         * @param state_u Uncorrected state variables.
         * @param sigma_prime_c Corrected effective stress tensor.
         * @param state_c Corrected state variables.
         * @param ITS_YSC Number of yield surface corrections applied.
         */
        void check_yield_surface_drift(Cauchy sigma_prime_u, State state_u, Cauchy &sigma_prime_c, State &state_c, int &ITS_YSC);
        
        /**
         * @brief Virtual method to compute the elastic state variable increment. 
         * 
         * @note Can be optionally overriden, otherwise does nothing but pass back the original state variable vector.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[in] state State variables.
         * @param[in] Delta_epsilon_tilde_e Elastic strain increment.
         * @return State 
         */
        virtual State compute_elastic_state_variable_increment(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde_e);

        /**
         * @brief Compute the elastic stress increment.
         * 
         * @param sigma_prime Effective stress tensor.
         * @param state State variables.
         * @param Delta_epsilon_tilde_e Elastic strain increment.
         * @param sigma_prime_e Effective stress after elastic increment.
         * @param state_e State variables after elastic increment.
         */
        void compute_elastic_increment(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde_e, Cauchy &Delta_sigma_prime_e, State &Delta_state_e);

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
         * @brief Get the settings object.
         * 
         * @return Settings 
         */
        Settings get_settings(void) {return this->settings;};

        /**
         * @brief Set the settings object.
         * 
         * @param settings 
         */
        void set_settings(Settings settings) {this->settings = settings;};

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
