#ifndef C2MC_H
#define C2MC_H

#include <plog/Log.h>
#include <Eigen/Eigen>
#include "Checks.hpp"
#include "Logging.hpp"
#include "Types.hpp"
#include "Elastoplastic.hpp"

class C2MC : public Elastoplastic {

    public: 

        /**
         * @brief C2MC model constructor. 
         * 
         * @param[in] parameters Vector of parameters.
         * @param[in] state Vector of state variables.
         * @param[in] log_severity Severity of message to log.
         */
        C2MC(State parameters, State state, std::string log_severity="none");

        /** 
         * @brief C2MC model destructor. 
         */
        virtual ~C2MC() {}

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
         * @param sigma_prime Effective stress state.
         * @param Delta_epsilon Strain increment.
         * @return Constitutive 
         */
        Constitutive compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon=Cauchy::Zero()) override;

        /**
         * @brief Overridden method to compute the current value of the yield surface function.
         * 
         * @param[in] sigma_prime Effective stress state.
         * @param[in] state State variables.
         * @return double
         */
        double compute_f(Cauchy sigma_prime, State state) override;

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
         * @brief Method to compute C2MC specific coefficients.
         * 
         * @param angle Friction or dilation angle in radians.
         * @param theta Lode's angle in radians.
         * @param A Coeffient A after Abbo et al. (2011).
         * @param B Coeffient B after Abbo et al. (2011).
         * @param C Coeffient C after Abbo et al. (2011).
         * @param k_theta Coeffient k_theta after Abbo et al. (2011).
         */
        void compute_coefficients(double angle_r, double theta, double &A, double &B, double &C, double &k_theta);

       /** 
         * @brief Parameters. 
         */
        Parameters parameters;

        /** 
         * @brief State variables. 
         */
        State state;

        /** 
         * @brief Parameter: Young's modulus. 
         */
        const double &E = parameters[0];
        
        /** 
         * @brief Parameter: Poisson's ratio. 
         */
        const double &nu = parameters[1];

        /** 
         * @brief Parameter: cohesion intercept.
         */
        const double &cohs = parameters[2];

        /** 
         * @brief Parameter: angle of internal friction. 
         */
        const double &phi = parameters[3];

        /** 
         * @brief Parameter: dilation angle. 
         */
        const double &psi = parameters[4];

        /** 
         * @brief Parameter: C2 rounding parameter. 
         */
        const double &theta_t = parameters[5];

        /** 
         * @brief Parameter: x-axis intercept for the hyperbolic approximation. 
         */
        const double &a_h = parameters[6];

        /**
         * @brief Number of required parameters.
         */
        int parameters_required = 7;

        /**
         * @brief Number of required state variables.
         */
        int state_required = 0;

        /**
         * @brief Friction angle in radians.
         */
        const double phi_r = to_radians(phi);

        /**
         * @brief Dilation angle in radians.
         */
        const double psi_r = to_radians(psi);

        /**
         * @brief C2 rounding parameter in radians.
         */
        const double theta_tr = to_radians(theta_t);

};

#endif 