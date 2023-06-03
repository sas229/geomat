#ifndef SMCC_H
#define SMCC_H

#include <plog/Log.h>
#include <Eigen/Eigen>
#include <cassert>
#include "Checks.hpp"
#include "Logging.hpp"
#include "Types.hpp"
#include "Elastoplastic.hpp"

class SMCC : public Elastoplastic {

    public: 

        /**
         * @brief SMCC model constructor. 
         * 
         * @param[in] parameters Vector of parameters.
         * @param[in] state Vector of state variables.
         * @param[in] log_severity Severity of message to log.
         */
        SMCC(State parameters, State state, std::string log_severity="none");

        /** 
         * @brief SMCC model destructor. 
         */
        virtual ~SMCC() {}

    protected:

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
         * @brief Overridden method to compute the current value of the yield surface function.
         * 
         * \f[ f = q^2 + M^2 p^{\prime}\left( p^{\prime}-p_c*s_ep \right)\f]
         * 
         * where \f$ q \f$ is the deviatoric stress, \f$ p^{\prime} \f$ is the mean effective stress, 
         * \f$ M\f$ is the frictional constant and \f$ p_c \f$ is the preconsolidation pressure.
         * 
         * @param[in] sigma_prime Effective stress state.
         * @param[in] state State variables.
         * @return Yield function, f. 
         */
        double compute_f(Cauchy sigma_prime, State state) override;

        /**
         * @brief Overridden method to compute the elastic stress state given an increment of strain.
         * 
         * @param sigma_prime Effective stress tensor.
         * @param alpha Fraction of strain increment to apply.
         * @param Delta_epsilon_tilde Strain increment.
         * @return Cauchy
         */
        Cauchy compute_elastic_stress(Cauchy sigma_prime, double alpha, Voigt Delta_epsilon_tilde) override;

        /**
         * @brief Overridden method to compute the elastic constitutive matrix.
         * 
         * @param sigma_prime Effective stress tensor.
         * @param Delta_epsilon_e_vol Elastic volumetric strain.
         * @return Constitutive 
         */
        Constitutive compute_elastic_matrix(Cauchy sigma_prime, double Delta_epsilon_e_vol) override;

        /**
         * @brief Overridden method to compute the bulk modulus.
         *  
         * If the increment exhibits non-zero volumetric strain:
         * 
         * \f[ K = \left(p^{\prime}/\Delta \epsilon_{e}^{vol}\right)  \left( \exp \left( \Delta \epsilon_{e}^{vol}/\kappa^{*} \right) -1 \right) \f]
         * 
         * where \f$ p^{\prime} \f$ is the effective mean stress, \f$ \Delta \epsilon_{e}^{vol} \f$ is the elastic volumetric strain increment
         * and \f$ \kappa^{*} \f$ is the slope of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
         * 
         * If there is no volumetric strain then the bulk modulus is defined as the tangent bulk modulus:
         * 
         * \f[ K = \frac{p^{\prime}}{\kappa^{*}} \f]
         * 
         * where \f$ p^{\prime} \f$ is the mean effective stress and \f$ \kappa^{*} \f$ is the slope 
         * of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
         *  
         * @param[in] Delta_epsilon_e_vol Elastic volumetric strain increment.
         * @param[in] p_prime Mean effective stress.
         * @return Bulk modulus, K.
         */
        double compute_K(double Delta_epsilon_e_vol, double p_prime) override;

        /**
         * @brief Overriden method to compute the shear modulus.
         * 
         * \f[ G = \frac{3 \left( 1-2\nu \right)K}{2 \left( 1+\nu \right)}\f]
         * 
         * where \f$ \nu \f$ is Poisson's ratio and \f$ K \f$ is the bulk modulus.
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
         * @param[in,out] H Hardening modulus.
         */
        void compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H) override;

        /**
         * @brief Overriden method to compute the elastic update of the models state variables.
         * 
         * @param[in] Delta_epsilon_tilde_e Elastic strain increment. 
         * @return Vector of state variables.
         */
        State compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) override;

        /**
         * @brief Overriden method to compute the plastic increment in the models state variables.
         * 
         * @param[in] Delta_epsilon_tilde_p Plastic strain increment.
         * @param[in] delta_lambda Plastic multiplier increment.
         * @param[in] df_dsigma_prime Derivatives of yield function with respect to the stress state.
         * @param[in] H Hardening modulus.
         * @return Vector of state variable increments.
         */
        State compute_plastic_state_variable_increment(Voigt Delta_epsilon_tilde_p, double delta_lambda, Cauchy df_dsigma_prime, double H) override;

        /**
         * @brief Overriden method to compute the correction in the models state variables.
         * 
         * @param[in] delta_lambda Plastic multiplier.
         * @param[in] df_dsigma_prime Derivatives of yield function with respect to the stress state.
         * @param[in] H Hardening modulus.
         * @return Vector of state variable corrections.
         */
        State compute_plastic_state_variable_increment(double delta_lambda, Cauchy df_dsigma_prime, double H) override;

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
         * @brief Parameter: initial sensitivity. 
         */
        const double &s_ep_i = parameters[5];

        /** 
         * @brief Parameter: sensitivity degradation rate parameter. 
         */
        const double &k = parameters[6];

        /** 
         * @brief Parameter: deformation balance parameter. 
         */
        const double &A = parameters[7];

        /** 
         * @brief State variable: voids ratio. 
         */
        double &e  = state[0];

        /** 
         * @brief State variable: preconsolidation pressure. 
         */
        double &p_c = state[1];

        /** 
         * @brief State variable: current sensitivity. 
         */
        double &s_ep = state[2];

        /**
         * @brief Number of required parameters.
         */
        int parameters_required = 8;

        /**
         * @brief Number of required state variables.
         */
        int state_required = 3;

};

#endif 