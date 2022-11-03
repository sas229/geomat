#ifndef ELASTOPLASTIC_H
#define ELASTOPLASTIC_H

#include <iostream>
#include <plog/Log.h>
#include <cmath>
#include <limits>
#include <cassert>
#include <Eigen/Eigen>
#include "Types.hpp"
#include "Elastic.hpp"

/**
 * @brief Elastoplastic base class. Inherits Elastic class.
 */
class Elastoplastic : public Elastic {

    protected: 

        /** 
         * @brief Elastoplastic model constructor.
         */
        Elastoplastic() {}

        /**
         *  @brief Elastoplastic model destructor. 
         */
        virtual ~Elastoplastic() {}
        
        /**
         *  @brief Solve current strain increment using refined explicit approach of Sloan et al. (2001). 
         */
        void solve(void);

        bool check_unload_reload(Cauchy sigma_prime);

        /**
         * @brief Method to compute the elastic fraction of the current strain increment following Sloan et al. (2001).
         * 
         * @return double 
         */
        double compute_alpha(double alpha_0, double alpha_1, double f_0, double f_1);

        /**
         * @brief Method to compute bounds for alpha for elastoplastic unloading-reloading increment.
         * 
         * @param alpha_0 Lower bound for alpha.
         * @param alpha_1 Upper bound for alpha.
         */
        void compute_alpha_bounds(double &alpha_0, double &alpha_1);

        /**
         * @brief Pure virtual method to compute the yield surface value given the parameters, current state variables and stress state.
         * 
         * @param sigma_prime Effective stress tensor. 
         * @return f 
         * 
         * @note Must be overriden by model implementations.
         */
        virtual double compute_f(Cauchy sigma_prime) = 0;

        /** 
         * @brief Pure virtual method to compute the derivatives for the constitutive model implemented.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[in,out] df_dsigma_prime Derivatives of yield function with respect to the stress state.
         * @param[in,out] dg_dsigma_prime Derivatives of plastic potential function with respect to the stress state.
         * @param[in,out] dg_dp_prime Derivative of plastic potential function with respect to the effective mean stress.
         * @param[in,out] H Hardening modulus.
         * @param[in,out] B_state Vector of state variable update scalars.
         * 
         * @note Must be overriden by model implementations.
         */
        virtual void compute_derivatives(Cauchy sigma_prime, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b,double &dg_dp_prime, double &H) = 0;

    private:

        /**
         * @brief Elastic fraction of strain increment.
         */
        double alpha;

        /**
         * @brief Derivatives of yield function with respect to the stress state.
         */
        Cauchy df_dsigma_prime;

        /**
         * @brief Vector of derivatives of yield function with respect to the stress state.
         */
        Voigt a;

        /**
         * @brief Derivatives of plastic potential function with respect to the stress state.
         */
        Cauchy dg_dsigma_prime;

        /**
         * @brief Vector of derivatives of plastic potential function with respect to the stress state.
         */
        Voigt b;

        /**
         * @brief Derivative of plastic potential function with respect to the effective mean stress.
         */
        double dg_dp_prime;

        /**
         * @brief Hardening modulus.
         */
        double H;

        /**
         * @brief Yield surface tolerance.
         */
        double FTOL = 1e-10;

        /**
         * @brief Maximum number of Pegasus method iterations to be performed during yield surface intersection calculations.
         */
        int MAXITS = 10;

        /**
         * @brief Unload-reload tolerance.
         */
        double LTOL = 1e-6;

        /**
         * @brief Number of unload-reload intersection substeps for alpha bound determination.
         */
        int NSUB = 10;

        /**
         * @brief Elastoplastic constitutive matrix.
         */
        Constitutive D_ep;
};

#endif
