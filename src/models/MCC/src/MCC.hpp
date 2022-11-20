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
         */
        MCC(State parameters, State state);

        /** 
         * @brief MCC model destructor. 
         */
        virtual ~MCC() {}

        /**
         * @brief Overridden method to compute the current value of the yield surface function.
         * 
         * @param sigma_prime Effective stress state.
         * @param state State variables.
         * @return f 
         */
        double compute_f(Cauchy sigma_prime, State state) override;

        /**
         * @brief Overridden method to compute the bulk modulus.
         * 
         * @param delta_epsilon_e_vol Elastic volumetric strain increment.
         * @param p_prime Mean effective stress.
         * @return K
         */
        double compute_K(double delta_epsilon_e_vol, double p_prime) override;

        /**
         * @brief Overriden method to compute the shear modulus.
         * 
         * @param K Bulk modulus.
         * @return double 
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
         * @brief Overriden method to compute the plastic increment in the models state variables.
         * 
         * @param delta_lambda Plastic multiplier increment.
         * @param H Hardening modulus.
         * @return STate variable increment.
         */
        State compute_plastic_state_variable_increment(double delta_lambda, double H) override;

        /**
         * @brief Overriden method to compute the correction in the models state variables.
         * 
         * @param delta_lambda 
         * @param H 
         * @return State variable correction
         */
        State compute_plastic_state_variable_correction(double delta_lambda, double H) override;

        /**
         * @brief Overriden method to compute the plastic update of the models state variables.
         */
        void compute_plastic_state_variable(void) override;

    protected:
       
        /** 
         * @brief Parameters. 
         */
        State parameters {{0.0, 0.0, 0.0, 0.0, 0.0}};

        /** 
         * @brief State variables. 
         */
        State state {{0.0, 0.0}};
    
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