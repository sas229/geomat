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

        /**
         * @brief Pure virtual method to compute the yield surface value given the parameters, current state variables and stress state.
         * 
         * @param p_prime Mean effective stress.
         * @param q Deviatoric stress.
         * @param state State variables.
         * @return double 
         */
        virtual double compute_f(double p_prime, double q, State state) = 0;

        /** 
         * @brief Pure virtual method to compute the derivatives for the constitutive model implemented.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[in] state State variables.
         * @param[in,out] df_dsigma_prime Tensor of derivatives of the yield function with respect to the stress state.
         * @param[in,out] a Vector of derivatives of the yield function with respect to the stress state.
         * @param[in,out] dg_dsigma_prime Derivatives of the plastic potential function with respect to the stress state.
         * @param[in,out] b Vector of derivatives of the plastic potential function with respect to the stress state.
         * @param[in,out] dg_dp_prime Derivative of the plastic potential function with respect to the mean effective stress.
         * @param[in,out] H Hardening modulus.
         * 
         * @note Must be overriden by model implementations.
         */
        virtual void compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &dg_dp_prime, double &H) = 0;

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
         * @param H Hardening modulus.
         * @return Vector of updated state variables.
         * 
         * @note Must be overriden by model implementations.
         */
        virtual State compute_plastic_state_variable_increment(Voigt Delta_epsilon_tilde_p, double delta_lambda, double H) = 0;

        /**
         * @brief Pure virtual method to compute the correction for the state variables for the model implemented (i.e. where the strain increment is constant).
         * 
         * @param delta_lambda Plastic multiplier.
         * @return Vector of state variable corrections.
         * 
         * @note Must be overriden by model implementations.
         */
        virtual State compute_plastic_state_variable_increment(double delta_lambda, double H) = 0;

        /**
         * @brief Pure virtual method to get the state variables from the model implementation.
         * 
         * @return Vector of state variables.
         *  
         * @note Must be overriden by model implementations.
         */
        virtual State get_state_variables(void) = 0;

        /**
         * @brief Pure virtual method to set the state variables in the model implementation.
         * 
         * @param Vector of state variables.
         *  
         * @note Must be overriden by model implementations.
         */
        virtual void set_state_variables(State new_state) = 0;

    // private:

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
         * @return Elastoplastic multiplier.
         */
        double compute_elastoplastic_multiplier(Voigt Delta_sigma_prime_e, Constitutive D_e, Voigt a, Voigt b, double H);

        /**
         * @brief Compute error estimate.
         */
        void compute_error_estimate(void);

        /**
         * @brief Compute new substep size.
         * 
         * @return Substep size. 
         */
        double compute_new_substep_size(void);

        /**
         * @brief Compute yield surface correction.
         */
        void compute_yield_surface_correction(void);

        /**
         * @brief Compute normal yield surface correction.
         */
        void compute_normal_yield_surface_correction(void);
        
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
         * @brief Method to compute vectors of coordinates to visualise the yield surface for a given set of state variables.
         * 
         * @param state Vector of state variables.
         * @param points Number of points on surface.
         * @param p_prime_surface Surface coordinate for mean effective stress.
         * @param q_surface Surface coordinate for deviatoric stress.
         */
        void compute_yield_surface(State state, int points, Eigen::VectorXd &p_prime_surface, Eigen::VectorXd &q_surface);

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
         * @brief Derivative of the yield surface with respect to the deviatoric stress.
         */
        double df_dq;

        /**
         * @brief Derivative of the yield surface with respect to the mean effective stress.
         */
        double df_dp_prime;

        /**
         * @brief Derivative of the yield surface with respect to the Lode angle.
         */
        double df_dtheta;

        /**
         * @brief Derivative of the plastic potential function with respect to the deviatoric stress.
         */
        double dg_dq;

        /**
         * @brief Derivative of the plastic potential function with respect to the effective mean stress.
         */
        double dg_dp_prime;

        /**
         * @brief Derivative of the plastic potential function with respect to the Lode angle.
         */
        double dg_dtheta;

        /**
         * @brief Hardening modulus.
         */
        double H;

        /**
         * @brief Elastoplastic constitutive matrix.
         */
        Constitutive D_ep;

        /**
         * @brief Pseudo-time for stress integration procedure.
         */
        double T;

        /**
         * @brief Pseudo-time increment for stress integration procedure.
         */
        double dT;

        /**
         * @brief Minimum pseudo-time increment for stress integration procedure.
         */
        double dT_min = 1e-6;

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

        /**
         * @brief Boolean indicating whether the current substep size was accepted.
         */
        bool accepted;
        
        /**
         * @brief Number of stress integration substeps performed.
         */
        int substeps;

        /**
         * @brief Number of yield surface corrections performed during current strain increment.
         */
        int corrections;

        /**
         * @brief Error estimate for current elastoplastic stress integration increment.
         */
        double R_n;

        /**
         * @brief Mean stress for uncorrected stress state.
         */
        double p_prime_u;

        /**
         * @brief Tangent bulk modulus for uncorrected stress state.
         */
        double K_tan_u;

        /**
         * @brief Tangent shear modulus for uncorrected stress state.
         */
        double G_tan_u;

        /**
         * @brief Elastic constitutive matrix for uncorrected stress state.
         */
        Constitutive D_e_u;

        /**
         * @brief Derivatives of the yield surface with respect to the uncorrected stress state.
         */
        Cauchy df_dsigma_prime_u;

        /**
         * @brief Derivatives of the plastic potential function with respect to the uncorrected stress state.
         */
        Cauchy dg_dsigma_prime_u;

        /**
         * @brief Vector of derivatives the yield surface with respect to the uncorrected stress state.
         */
        Voigt a_u;

        /**
         * @brief Vector of derivatives the plastic potential function with respect to the uncorrected stress state.
         */
        Voigt b_u;

        /**
         * @brief Derivative of the plastic potential function with respect to the mean stress.
         */
        double dg_dp_prime_u;

        /**
         * @brief Hardening modulus for uncorrected stress state.
         */
        double H_u;

        /**
         * @brief Plastic multiplier for the yield surface correction.
         */
        double delta_lambda_c;

        /**
         * @brief Tensor of stress state corrections.
         */
        Voigt Delta_sigma_prime_c;

        /**
         * @brief Vector of state variable corrections.
         * 
         */
        State delta_state_c;

        /**
         * @brief Substep size factor.
         */
        double q_step;

        /**
         * @brief Elastic volumetric strain increment.
         */
        double Delta_epsilon_vol_e;

        /**
         * @brief Elastic strain increment.
         */
        Cauchy Delta_epsilon_e;
        
        /**
         * @brief Yield surface function for uncorrected stress state.
         */
        double f_u;

        /**
         * @brief Yield surface function for corrected stress state.
         */
        double f_c;

        /**
         * @brief Intial effective stress estimate from forward Euler method.
         */
        Cauchy sigma_prime_ini;

        /**
         * @brief Effective stress after applying elastic portion of strain increment.
         */
        Cauchy sigma_prime_e;

        /**
         * @brief Effective stress after applying elastoplastic strain increment.
         */
        Cauchy sigma_prime_ep;

        /**
         * @brief Uncorrected effective stress state.
         */
        Cauchy sigma_prime_u;

        /**
         * @brief Corrected effective stress state.
         */
        Cauchy sigma_prime_c;

        /**
         * @brief Plastic volumetric strain increment for pseudo-time increment dT.
         */
        double Delta_epsilon_vol_p_dT;

        /**
         * @brief Elastic strain increment for pseudo-time increment dT.
         */
        Voigt Delta_epsilon_tilde_p_dT;

        /**
         * @brief Elastoplastic multiplier.
         */
        double delta_lambda;

        /**
         * @brief Effective stress after applying plastic portion of strain increment for pseudo-time increment dT.
         */
        Cauchy sigma_prime_p_dT;

        /**
         * @brief Elastic strain increment.
         */
        Voigt Delta_epsilon_tilde_e;

        /**
         * @brief Plastic strain increment.
         */
        Voigt Delta_epsilon_tilde_p;

        /**
         * @brief Stress state after application of elastic increment.
         */
        Voigt Delta_sigma_prime_e;

        /**
         * @brief Intial state variable estimate from forward Euler method.
         */
        State state_ini;

        /**
         * @brief State variables after application of elastic increment.
         */
        State state_e;

        /**
         * @brief State variables after application of elastoplastic increment.
         */
        State state_ep;

        /**
         * @brief Uncorrected state variables.
         */
        State state_u;

        /**
         * @brief Corrected state variables.
         */
        State state_c;

        /**
         * @brief State variables after first forward Euler estimate.
         */
        State state_1;

        /**
         * @brief State variables after second forward Euler estimate.
         */
        State state_2;

        /**
         * @brief Effective stress for the first estimate in the forward Euler method.
         * 
         */
        Cauchy sigma_prime_1;

        /**
         * @brief Effective stress increment for the first estimate in the forward Euler method.
         */
        Voigt Delta_sigma_prime_1;

        /**
         * @brief State variable increment for the first estimate in the forward Euler method.
         */
        State delta_state_1;

        /**
         * @brief Effective stress for the second estimate in the forward Euler method.
         */
        Cauchy sigma_prime_2;

        /**
         * @brief Effective stress increment for the second estimate in the forward Euler method.
         * 
         */
        Voigt Delta_sigma_prime_2;

        /**
         * @brief State variable increment for the second estimate in the forward Euler method.
         */
        State delta_state_2;
        
        /**
         * @brief Effective stress after applying plastic portion of strain increment prior to correction.
         * 
         */
        Cauchy sigma_prime_ep_ini;
};

#endif
