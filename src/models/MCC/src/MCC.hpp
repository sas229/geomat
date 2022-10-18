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
         * @param delta_epsilon_vol Elastic volumetric strain.
         * @param p_prime Mean effective stress.
         * @return K
         */
        double compute_K(double delta_epsilon_vol, double p_prime);

        /**
         * @brief Overrident method to compute the shear modulus.
         * 
         * @param K Bulk modulus.
         * @return double 
         */
        double compute_G(double K);

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