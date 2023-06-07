#ifndef SMCC_H
#define SMCC_H

#include <plog/Log.h>
#include <Eigen/Eigen>
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

        /**
         * @brief Get the state variables vector.
         * 
         * @return State
         */
        State get_state_variables(void) override;

    protected:

        /**
         * @brief Set the state variables vector.
         * 
         * @param[in] new_state Vector of new state variables.
         */
        void set_state_variables(State new_state) override;
 
        /**
         * @brief Overriden method to compute the elastic constitutive matrix.
         * 
         * An isotropic stress-dependent linear elastic constitutive matrix is used in the MCC model:
         * 
         * \f[ D_e = \left[\begin{array}{cccccc}
                     K + \frac{4}{3}G & K - \frac{2}{3}G & K - \frac{2}{3}G & 0 & 0 & 0 \\
                     K - \frac{2}{3}G & K + \frac{4}{3}G & K - \frac{2}{3}G & 0 & 0 & 0 \\
                     K - \frac{2}{3}G & K - \frac{2}{3}G & K + \frac{4}{3}G & 0 & 0 & 0 \\
                     0 & 0 & 0 & G & 0 & 0 \\
                     0 & 0 & 0 & 0 & G & 0 \\
                     0 & 0 & 0 & 0 & 0 & G
                     \end{array}\right] \f]
         * where \f$ K \f$ is the bulk modulus and \f$ G \f$ is the shear modulus.
         * 
         * The secant bulk modulus is defined as:
         * 
         * \f[ K = \left(p^{\prime}/\Delta \epsilon_{e}^{vol}\right)  \left( \exp \left( \Delta \epsilon_{e}^{vol}/\kappa^{*} \right) -1 \right)    \f]
         * 
         * where \f$ p^{\prime} \f$ is the effective mean stress, \f$ \Delta \epsilon_{e}^{vol} \f$ is the elastic volumetric strain increment
         * and \f$ \kappa^{*} \f$ is the slope of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
         * 
         * The tangent bulk modulus is defined as:
         * 
         * \f[ K = \frac{p^{\prime}}{\kappa^{*}} \f]
         * 
         * where \f$ p^{\prime} \f$ is the mean effective stress and \f$ \kappa^{*} \f$ is the slope 
         * of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
         * 
         * The shear modulus is defined as:
         * 
         * \f[ G = \frac{3 \left( 1-2\nu \right)K}{2 \left( 1+\nu \right)}\f]
         * 
         * where \f$ \nu \f$ is Poisson's ratio and \f$ K \f$ is the bulk modulus.
         * 
         * @param sigma_prime Effective stress state.
         * @param Delta_epsilon Strain increment.
         * @return Constitutive 
         */
        Constitutive compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon=Cauchy::Zero()) override;

        /**
         * @brief Overridden method to compute the current value of the yield surface function.
         * 
         * Yield surface function definition for the SMCC model:
         * 
         * \f[ f = q^2 + M^2 p^{\prime}\left( p^{\prime}-p_c s_{ep} \right)\f]
         * 
         * where \f$ q \f$ is the deviatoric stress, \f$ p^{\prime} \f$ is the mean effective stress, 
         * \f$ M\f$ is the frictional constant and \f$ p_c \f$ is the preconsolidation pressure.
         * 
         * @param[in] sigma_prime Effective stress state.
         * @param[in] state State variables.
         * @return double
         */
        double compute_f(Cauchy sigma_prime, State state) override;

        /** 
         * @brief Method to compute the derivatives for the constitutive model implemented.
         * 
         * Derivative of the yield surface with respect to the deviatoric stress:
         * 
         * \f[ \frac{\partial f}{\partial q} = 2q \f]
         * 
         * where \f$ q \f$ is the deviatoric stress.
         * 
         * Derivative of the plastic potential function with respect to the Lode angle:
         * 
         * \f[ \frac{\partial f}{\partial \theta} = 0 \f]
         * 
         * Derivatives of the yield surface with respect to the mean effective stress:
         * 
         * \f[ \frac{\partial f}{\partial p^{\prime}} = M^2\left(2 p^{\prime}-p_{c} s_{ep} \right) \f]
         * 
         * where \f$ M \f$ is a frictional constant, \f$ p^{\prime} \f$ is the mean effective stress
         * and \f$ p_{c} \f$ is the preconsolidation pressure.
         * 
         * Derivative of the plastic potential function with respect to the deviatoric stress:
         * 
         * \f[ \frac{\partial g}{\partial q} = 2q \f]
         * 
         * where \f$ q \f$ is the deviatoric stress.
         * 
         * Derivative of the plastic potential function with respect to the Lode angle:
         * 
         * \f[ \frac{\partial g}{\partial \theta} = 0 \f]
         * 
         * Derivative of the plastic potential function with respect to the mean effective stress:
         * 
         * \f[ \frac{\partial g}{\partial p^{\prime}} = M^2 \left( 2 p^{\prime}-p_c  s_{ep} \right) \f]
         * 
         * where \f$ M \f$ is a frictional constant, \f$ p^{\prime} \f$ is the mean effective stress
         * and \f$ p_{c} \f$ is the preconsolidation pressure.
         * 
         * The hardening modulus is comprised of two components:
         * 
         * \f[ H = H_{p_{c}} + H_{s_{ep}} \f]
         * 
         * where:
         *  
         * \f[ H_{p_{c}} = \frac{M^2 p^{\prime} p_c}{\lambda^*-\kappa^*} s_{ep} \operatorname{tr}\left( \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} \right) \f]
         * 
         * and:
         * 
         * \f[ H_{s_{ep}} = -\frac{M^2 p^{\prime} p_c}{\lambda^*-\kappa^*} k \left(s_{ep}-1\right) \sqrt{\left(1-A\right) {\operatorname{tr}\left( \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} \right)}^2 
         * + A \frac{2}{3} \operatorname{tr} \left( \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} {\frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}}}^T \right) } \f]
         * 
         * where \f$ M \f$ is a frictional constant,\f$ p^{\prime} \f$ is the mean effective stress,
         * \f$ p_{c} \f$ is the preconsolidation pressure, \f$ \lambda^* \f$ and \f$ \kappa^* \f$ are the slopes of the normal compression 
         * and re-compression lines in \f$ \ln \left(e\right)-\ln \left(p^{\prime}\right)\f$ space.
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
         * @return State
         */
        State compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) override;

        /**
         * @brief Overriden method to compute the plastic increment in the models state variables.
         * 
         * @param[in] Delta_epsilon_tilde_p Plastic strain increment.
         * @param[in] delta_lambda Plastic multiplier increment.
         * @param[in] df_dsigma_prime Derivatives of yield function with respect to the stress state.
         * @param[in] H Hardening modulus.
         * @return State
         */
        State compute_plastic_state_variable_increment(double delta_lambda, Cauchy df_dsigma_prime, double H, Voigt Delta_epsilon_tilde_p=Voigt::Zero()) override;

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
