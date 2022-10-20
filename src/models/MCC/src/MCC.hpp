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
        MCC(std::vector<double> parameters, std::vector<double> state);

        /** 
         * @brief MCC model destructor. 
         */
        virtual ~MCC() {}

        /**
         * @brief Overridden method to compute the current value of the yield surface function.
         * 
         * @return f
         */
        double compute_f(Cauchy sigma_prime) override;

        /**
         * @brief Overridden method to compute the bulk modulus.
         * 
         * @param delta_epsilon_e_vol Elastic volumetric strain increment.
         * @param p_prime Mean effective stress.
         * @return K
         */
        double compute_K(double delta_epsilon_vol, double p_prime);

        /**
         * @brief Overriden method to compute the shear modulus.
         * 
         * @param K Bulk modulus.
         * @return double 
         */
        double compute_G(double K);

        /** 
         * @brief Method to compute the derivatives for the constituive model implemented.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[in,out] df_dsigma_prime Derivatives of yield function with respect to the stress state.
         * @param[in,out] dg_dsigma_prime Derivatives of plastic potential function with respect to the stress state.
         * @param[in,out] dg_dp_prime Derivative of plastic potential function with respect to the effective mean stress.
         * @param[in,out] H Hardening modulus.
         * @param[in,out] B_state Vector of state variable update scalars.
         */
        void compute_derivatives(Cauchy sigma_prime, Cauchy &df_dsigma_prime, Cauchy &dg_dsigma_prime, double &dg_dp_prime, double &H, std::vector<double> &B_state);

    protected:
       
        /** 
         * @brief Parameters. 
         */
        std::vector<double> parameters {0.0, 0.0, 0.0, 0.0, 0.0};

        /** 
         * @brief State variables. 
         */
        std::vector<double> state {0.0, 0.0};
    
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