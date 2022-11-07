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

        /**
         * @brief Pure virtual method to compute the elastic update of the state variables for the model implemented.
         */
        virtual void compute_elastic_state_variable_update(void) = 0;

        /**
         * @brief Elastic volumetric strain increment.
         */
        double delta_epsilon_vol_e;

    private:

        /**
         * @brief Pegasus regula falsi algorithm to find root of general non-linear system of equations.
         * 
         * @param alpha_0 Lower bound on alpha.
         * @param alpha_1 Upper bound on alpha.
         * @param f_0 Initial value of objective function with lower bound alpha.
         * @param f_1 Initial value of objective function with upper bound alpha.
         * @return alpha
         */
        double pegasus_regula_falsi(double alpha_0, double alpha_1, double f_0, double f_1);

        /**
         * @brief Method to determine if an increment is an unload-reload plastic increment.
         * 
         * @param sigma_prime Current stress state.
         * @return true
         * @return false 
         */
        bool check_unload_reload(Cauchy sigma_prime);

        /**
         * @brief Method to compute the elastic fraction of the current strain increment following Sloan et al. (2001).
         * 
         */
        void compute_alpha(void);

        /**
         * @brief Method to compute bounds for alpha for elastoplastic unloading-reloading increment.
         * 
         * @param alpha_0 Lower bound for alpha.
         * @param alpha_1 Upper bound for alpha.
         */
        void compute_alpha_bounds(double &alpha_0, double &alpha_1);

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

        /**
         * @brief Elastic strain increment.
         */
        Voigt delta_epsilon_tilde_e;

        /**
         * @brief Plastic strain increment.
         */
        Voigt delta_epsilon_tilde_p;

        /**
         * @brief Effective stress after applying elastic portion of strain increment.
         */
        Cauchy sigma_prime_e;
};

#endif
