#ifndef MCC_H
#define MCC_H

#include "Elastoplastic.hpp"
#include "Tensor.hpp"
#include "Types.hpp"

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
        State get_state_variables(void) {return state;};

        void set_state_variables(State new_state) {this->state = new_state;};

    protected:

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
         * @return Derivatives
         */
        Derivatives compute_derivatives(Cauchy sigma_prime, State state) override;

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
         * @brief State variable: preconsolidation pressure. 
         */
        double &p_c = state[0];

        /**
         * @brief Number of required parameters.
         */
        int parameters_required = 5;

        /**
         * @brief Number of required state variables.
         */
        int state_required = 1;

};

#endif 
