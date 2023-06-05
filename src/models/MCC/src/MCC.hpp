#ifndef MCC_H
#define MCC_H

#include <plog/Log.h>
#include <Eigen/Eigen>
#include "Types.hpp"
#include "Elastoplastic.hpp"
#include "Logging.hpp"
#include "Checks.hpp"

class MCC : public Elastoplastic {

    public: 

        /**
         * @brief MCC model constructor. 
         * 
         * @param[in] parameters Vector of parameters.
         * @param[in] state Vector of state variables.
         * @param[in] log_severity Severity of message to log.
         */
        MCC(State parameters, State state, std::string log_severity="none");

        /** 
         * @brief MCC model destructor. 
         */
        virtual ~MCC() {}

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
         * Yield surface function definition for the MCC model:
         * 
         * \f[ f = q^2 + M^2 p^{\prime}\left( p^{\prime}-p_c \right)\f]
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
         * \f[ \frac{\partial f}{\partial p^{\prime}} = M^2\left(2 p^{\prime}-p_{c} \right) \f]
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
         * \f[ \frac{\partial g}{\partial p^{\prime}} = M^2 \left( 2 p^{\prime}-p_c \right) \f]
         * 
         * where \f$ M \f$ is a frictional constant, \f$ p^{\prime} \f$ is the mean effective stress
         * and \f$ p_{c} \f$ is the preconsolidation pressure.
         * 
         * The hardening modulus is:
         *  
         * \f[ H = \frac{M^2 p^{\prime} p_c}{\lambda^*-\kappa^*} \operatorname{tr}\left( \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} \right) \f]
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
         * State variable elastic update for void ratio:
         * 
         * \f[ e = e - e \Delta \epsilon_{vol, e}\f]
         *
         * State variable elastic update for preconsolidation pressure:
         * 
         * \f[ p_{c} = constant \f]
         *
         * @param[in] Delta_epsilon_tilde_e Elastic strain increment. 
         * @return State
         */
        State compute_elastic_state_variable(Voigt Delta_epsilon_tilde_e) override;

        /**
         * @brief Overriden method to compute the plastic increment in the models state variables.
         * 
         * State variable plastic increment for void ratio:
         * 
         * \f[ \Delta e = -(1+e) \Delta \epsilon_{vol, p} \f]
         * 
         * where \f$ e \f$ is the void ratio and \f$ \Delta \epsilon_{vol, p} \f$ 
         * is the plastic volumetric strain.
         * 
         * State variable plastic increment for preconsolidation pressure:
         * 
         * /f[ \Delta p_{c} = \Delta \lambda \frac{H}{M^2 p^{\prime}} /f]
         * 
         * where \f$ \Delta \lambda \f$ is the plastic multiplier, \f$ H \f$ is the hardening modulus, 
         * \f$ M \f$ is a frictional constant and \f$ p^{\prime} \f$ is the mean effective stress.
         * 
         * @param[in] delta_lambda Plastic multiplier increment.
         * @param[in] df_dsigma_prime Derivatives of yield function with respect to the stress state.
         * @param[in] H Hardening modulus.
         * @param[in] Delta_epsilon_tilde_p Plastic strain increment.
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
         * @brief State variable: voids ratio. 
         */
        double &e  = state[0];

        /** 
         * @brief State variable: preconsolidation pressure. 
         */
        double &p_c = state[1];

        /**
         * @brief Number of required parameters.
         */
        int parameters_required = 5;

        /**
         * @brief Number of required state variables.
         */
        int state_required = 2;

};

#endif 
