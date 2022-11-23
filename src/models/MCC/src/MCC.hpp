#ifndef MCC_H
#define MCC_H

#include <plog/Log.h>
#include <Eigen/Eigen>
#include <cassert>
#include "Types.hpp"
#include "Elastoplastic.hpp"

class MCC : public Elastoplastic {

    public: 

        /**
         * @brief MCC model constructor. 
         * 
         * @param[in] parameters Vector of parameters.
         * @param[in] state Vector of state variables.
         */
        MCC(State parameters, State state);

        /** 
         * @brief MCC model destructor. 
         */
        virtual ~MCC() {}

    protected:

        /**
         * @brief Overridden method to compute the current value of the yield surface function.
         * 
         * @param[in] sigma_prime Effective stress state.
         * @param[in] state State variables.
         * @return Yielf function, f. 
         */
        double compute_f(Cauchy sigma_prime, State state) override;

        /**
         * @brief Overridden method to compute the bulk modulus.
         * 
         * @param[in] delta_epsilon_e_vol Elastic volumetric strain increment.
         * @param[in] p_prime Mean effective stress.
         * @return Bulk modulus, K.
         */
        double compute_K(double delta_epsilon_e_vol, double p_prime) override;

        /**
         * @brief Overriden method to compute the shear modulus.
         * 
         * @param[in] K Bulk modulus.
         * @return Shear modulus, G.
         */
        double compute_G(double K) override;

        /** 
         * @brief Method to compute the derivatives for the constitutive model implemented.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[in] state State variables.
         * @param[in,out] df_dsigma_prime Derivatives of yield function with respect to the stress state.
         * @param[in,out] a Vector of derivatives of yield function with respect to the stress state.
         * @param[in,out] dg_dsigma_prime Derivatives of plastic potential function with respect to the stress state.
         * @param[in,out] b Vector of derivatives of plastic potential function with respect to the stress state.
         * @param[in,out] dg_dp_prime Derivative of plastic potential function with respect to the effective mean stress.
         * @param[in,out] H Hardening modulus.
         * @param[in,out] B_state Vector of state variable update scalars.
         */
        void compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &dg_dp_prime, double &H);

        /**
         * @brief Overriden method to compute the elastic update of the models state variables.
         * 
         * @param[in] delta_epsilon_tilde_e Elastic strain increment. 
         * @return Vector of state variables.
         */
        State compute_elastic_state_variable(Voigt delta_epsilon_tilde_e) override;

        /**
         * @brief Get the state variables vector.
         * 
         * @return Vector of state variables.
         */
        State get_state_variables(void);

        /**
         * @brief Set the state variables vector.
         * 
         * @param[in] new_state Vector of new state variables.
         */
        void set_state_variables(State new_state);

        /**
         * @brief Overriden method to compute the plastic increment in the models state variables.
         * 
         * @param[in] delta_epsilon_tilde_p Plastic strain increment.
         * @param[in] delta_lambda Plastic multiplier increment.
         * @param[in] H Hardening modulus.
         * @return Vector of state variable increments.
         */
        State compute_plastic_state_variable_increment(Voigt delta_epsilon_tilde_p, double delta_lambda, double H) override;

        /**
         * @brief Overriden method to compute the correction in the models state variables.
         * 
         * @param[in] delta_lambda Plastic multiplier.
         * @param[in] H Hardening modulus.
         * @return Vector of state variable corrections.
         */
        State compute_plastic_state_variable_increment(double delta_lambda, double H) override;

       /** 
         * @brief Parameters. 
         */
        Parameters parameters;

        /** 
         * @brief State variables. 
         */
        State state;

        /** 
         * @brief Parameter: frictional constant. 
         */
        const double &M = parameters[0];
        
        /** 
         * @brief Parameter: Poisson's ratio. 
         */
        const double &nu = parameters[1];

        /** 
         * @brief Parameter: intercept of NCL. 
         */
        const double &N = parameters[2];

        /** 
         * @brief Parameter: slope of the NCL in ln(e)-ln(p') space. 
         */
        const double &lambda_star = parameters[3];

        /** 
         * @brief Parameter: slope of the RCL in ln(e)-ln(p') space. 
         */
        const double &kappa_star = parameters[4];

        /** 
         * @brief State variable: voids ratio. 
         */
        double &e  = state[0];

        /** 
         * @brief State variable: preconsolidation pressure. 
         */
        double &p_c = state[1];

};

#endif 