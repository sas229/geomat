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

    public: 

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
         * @param state State variables.
         * @return f 
         * 
         * @note Must be overriden by model implementations.
         */
        virtual double compute_f(Cauchy sigma_prime, State state) = 0;

        // /**
        //  * @brief Pure virtual method to compute the yield surface value given the parameters, current state variables and stress state.
        //  * 
        //  * @param p_prime Mean effective stress.
        //  * @param q Deviatoric stress.
        //  * @param state State variables.
        //  * @return double 
        //  */
        // virtual double compute_f(double p_prime, double q, State state) = 0;

        /** 
         * @brief Pure virtual method to compute the derivatives for the constitutive model implemented.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[in] state State variables.
         * @param[in,out] df_dsigma_prime Tensor of derivatives of the yield function with respect to the stress state.
         * @param[in,out] a Vector of derivatives of the yield function with respect to the stress state.
         * @param[in,out] dg_dsigma_prime Derivatives of the plastic potential function with respect to the stress state.
         * @param[in,out] b Vector of derivatives of the plastic potential function with respect to the stress state.
         * @param[in,out] H Hardening modulus.
         * 
         * @note Must be overriden by model implementations.
         */
        virtual void compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H) = 0;

        /**
         * @brief Pure virtual method to compute the elastic update of the state variables for the model implemented.
         * 
         * @return Vector of state variables.
         * 
         * @note Must be overriden by model implementations.
         */
        virtual State compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) = 0;

        /**
         * @brief Pure virtual method to compute the increment in the state variables for the model implemented.
         * 
         * @param Delta_epsilon_tilde_p Plastic strain increment.
         * @param delta_lambda Plastic multiplier.
         * @param df_dsigma_prime Tensor of derivatives of the yield function with respect to the stress state.
         * @param H Hardening modulus.
         * @return Vector of updated state variables.
         * 
         * @note Must be overriden by model implementations.
         */
        virtual State compute_plastic_state_variable_increment(double delta_lambda, Cauchy df_dsigma_prime, double H, Voigt Delta_epsilon_tilde_p=Voigt::Zero()) = 0;

        /**
         * @brief Sloan et al. explicit substepping algorithm for solving the plastic portion of a strain increment.
         * 
         * @param sigma_prime_ep Effective stress state at the elastoplastic transition.
         * @param state_ep State variables at the elastoplastic transition.
         * @param Delta_epsilon_tilde_p Plastic strain increment.
         */
        void sloan_substepping(Cauchy sigma_prime_ep, State state_ep, Voigt Delta_epsilon_tilde_p);

        /**
         * @brief Compute elstoplastic constitutive matrix via:
         * 
         * /f[ \mathbf{D}_{\mathrm{ep}}=\mathbf{D}_{\mathrm{e}}-\frac{\mathbf{D}_{\mathrm{e}} \mathbf{b a}^{\mathrm{T}} \mathbf{D}_{\mathrm{e}}}{A+\mathbf{a}^{\mathrm{T}} \mathbf{D}_{\mathrm{e}} \mathbf{b}}/f]
         * 
         * @param D_e Elastic constitutive matrix.
         * @param a Vector of derivatives with respect to the yield surface.
         * @param b Vector of derivatives with respect to the plastic potential function.
         * @param H Hardening modulus.
         * @return Elastoplastic constitutive matrix.
         */
        Constitutive compute_elastoplastic_matrix(Constitutive D_e, Voigt a, Voigt b, double H);

        /**
         * @brief Compute elastoplastic multiplier.
         * 
         * @param Delta_sigma_prime_e Elastic stress increment.
         * @param D_e Elastic constitutive matrix.
         * @param a Vector of derivatives with respect to the yield surface.
         * @param b Vector of derivatives with respect to the plastic potential function.
         * @param H Hardening modulus.
         * @return double
         */
        double compute_elastoplastic_multiplier(Voigt Delta_sigma_prime_e, Constitutive D_e, Voigt a, Voigt b, double H);

        /**
         * @brief Compute error estimate.
         * 
         * @param sigma_prime_ini Effective stress estimate.
         * @param state_ini State variable estimate.
         * @param Delta_sigma_prime_1 First estimate of stress increment.
         * @param Delta_sigma_prime_2 Second estimate of stress increment.
         * @param delta_state_1 First estimate of state variable increment.
         * @param delta_state_2 Second estimate of state variable increment.
         * @return double
         */
        double compute_error_estimate(Cauchy sigma_prime_ini, State state_ini, Voigt Delta_sigma_prime_1, Voigt Delta_sigma_prime_2, State delta_state_1, State delta_state_2);

        /**
         * @brief Compute new substep size.
         * 
         * @param accepted Boolean to indicate if the increment was accepted.
         * @param dT Pseudotime increment.
         * @param T Current pseudotime.
         * @param R_n Error estimate.
         * @return Substep size. 
         */
        double compute_new_substep_size(bool accepted, double dT, double T, double R_n);

        /**
         * @brief Compute yield surface correction.
         * 
         * @param sigma_prime_u Uncorrected effective stress state.
         * @param state_u Uncorrected state variables.
         * @param f_u Yield surface value for the uncorrected state.
         * @param H_u Uncorrected hardening modulus.
         * @param a_u Vector of derivatives with respect to the yield surface for the uncorrected stress state.
         * @param b_u Vector of derivatives with respect to the plastic potential function for the uncorrected stress state.
         * @param D_e_u Uncorrected elastic matrix.
         * @param sigma_prime_c Corrected effective stress state.
         * @param state_c Corrected state variables.
         */
        void compute_yield_surface_correction(Cauchy sigma_prime_u, State state_u, double f_u, double H_u, Voigt a_u, Voigt b_u, Constitutive D_e_u, Cauchy &sigma_prime_c, State &state_c);

        /**
         * @brief Compute normal yield surface correction.
         * 
         * @param sigma_prime_u Uncorrected effective stress state.
         * @param state_u Uncorrected state variables.
         * @param f_u Yield surface value for the uncorrected stress state.
         * @param a_u Vector of derivatives with respect to the yield surface for the uncorrected stress state.
         * @param sigma_prime_c Corrected effective stress state.
         * @param state_c Corrected state variables.
         */
        void compute_normal_yield_surface_correction(Cauchy sigma_prime_u, State state_u, double f_u, Voigt a_u, Cauchy &sigma_prime_c, State &state_c);
        
        /**
         * @brief Compute plastic stress integration increment.
         * 
         * @param sigma_prime Effective stress state.
         * @param state Vector of state variables.
         * @param Delta_epsilon_tilde_p_dT Strain increment.
         * @param Delta_sigma_prime Increment in stress state.
         * @param delta_state Increment in state variables.
         */
        void compute_plastic_increment(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde_p_dT, Voigt &Delta_sigma_prime, State &delta_state);

        /**
         * @brief Pegasus regula falsi algorithm to find root of general non-linear system of equations.
         * 
         * @param sigma_prime Current stress state.
         * @param state Current state variables.
         * @param alpha_0 Lower bound on alpha.
         * @param alpha_1 Upper bound on alpha.
         * @param f_0 Initial value of objective function with lower bound alpha.
         * @param f_1 Initial value of objective function with upper bound alpha.
         * @return alpha
         */
        double pegasus_regula_falsi(Cauchy sigma_prime, State state, double alpha_0, double alpha_1, double f_0, double f_1);

        /**
         * @brief Method to determine if an increment is an unload-reload plastic increment.
         * 
         * @param sigma_prime Current stress state.
         * @param state Current state variables.
         * @return true
         * @return false 
         */
        bool check_unload_reload(Cauchy sigma_prime, State state);

        /**
         * @brief Method to compute the elastic fraction of the current strain increment following Sloan et al. (2001).
         * 
         * @param sigma_prime Current stress state.
         * @param state Current state variables.
         * @return alpha
         * 
         */
        double compute_alpha(Cauchy sigma_prime, State state);

        /**
         * @brief Method to compute bounds for alpha for elastoplastic unloading-reloading increment.
         * 
         * @param alpha_0 Lower bound for alpha.
         * @param alpha_1 Upper bound for alpha.
         */
        void compute_alpha_bounds(double &alpha_0, double &alpha_1);

        /**
         * @brief Minimum pseudo-time increment for stress integration procedure.
         */
        double DT_MIN = 1e-6;

        /**
         * @brief Yield surface tolerance.
         */
        double FTOL = 1e-8;

        /**
         * @brief Number of Pegasus method iterations performed during yield surface intersection calculations.
         */
        int ITS_YSI;

        /**
         * @brief Maximum number of Pegasus method iterations to be performed during yield surface intersection calculations.
         */
        int MAXITS_YSI = 10;

        /**
         * @brief Number of yield surface correction iterations performed during current elastoplastic stress integration substep.
         */
        int ITS_YSC;

        /**
         * @brief Maximum number of yield surface correction iterations to be performed during elastoplastic stress integration.
         */
        int MAXITS_YSC = 30;

        /**
         * @brief Unload-reload tolerance.
         */
        double LTOL = 1e-6;

        /**
         * @brief Stress integration tolerance.
         */
        double STOL = 1e-4;

        /**
         * @brief Number of unload-reload intersection substeps for alpha bound determination.
         */
        int NSUB = 10;

        /**
         * @brief Double precision tolerance.
         */
        double EPS = 1e-16;

};

#endif
