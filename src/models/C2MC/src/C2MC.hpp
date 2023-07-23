#ifndef C2MC_H
#define C2MC_H

#include "Checks.hpp"
#include "Elastoplastic.hpp"
#include "Logging.hpp"
#include "Types.hpp"


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
        State get_state_variables(void) {return state;};

    protected:

        /**
         * @brief Set the state variables vector.
         * 
         * @param[in] new_state Vector of new state variables.
         */
        void set_state_variables(State new_state) {state = new_state;};
 
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
         * @param[in,out] dg_dsigma_prime Derivatives of plastic potential function with respect to the stress state.
         * @param[in,out] H Hardening modulus.
         */
        void compute_derivatives(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Cauchy &dg_dsigma_prime, HardeningModuli &H_s, StateFactors &B_s) override;

        /**
         * @brief Method to compute the C2 continuous form of \f$ K \left( \theta \right) \f$ for the Mohr-Coulomb surface.
         * 
         * @param[in] angle_r Friction or dilation angle in radians.
         * @param[in] theta_s_bar Lode's angle in radians.
         * @param[in] A C2 continuous approximation coefficient after Abbo et al. (2011).
         * @param[in] B C2 continuous approximation coefficient after Abbo et al. (2011).
         * @param[in] C C2 continuous approximation coefficient after Abbo et al. (2011).
         * @param[in,out] K_theta \f$ K \left( \theta \right) \f$ for the Mohr-Coulomb surface.
         */
        void compute_K_theta(double angle_r, double theta_s_bar, double A, double B, double C, double &K_theta);

        /**
         * @brief Method to compute the C2 continuous Mohr-Coulomb model coefficients.
         * 
         * @param[in] angle_r Friction or dilation angle in radians.
         * @param[in] theta_s_bar Lode's angle in radians.
         * @param[in,out] A C2 continuous approximation coefficient \f$ A \f$.
         * @param[in,out] B C2 continuous approximation coefficient \f$ B \f$.
         * @param[in,out] C C2 continuous approximation coefficient \f$ C \f$.
         */
        void compute_A_B_C(double angle_r, double theta_s_bar, double &A, double &B, double &C);

        /**
         * @brief 
         * 
         * @param[in] sigma_bar Stress invariant.
         * @param[in] angle_r Friction or dilation angle in radians.
         * @param[in] theta_s_bar Lode's angle in radians.
         * @param[in,out] C_1 Derivative coefficient \f$ C_{1} \f$.
         * @param[in,out] C_2 Derivative coefficient \f$ C_{2} \f$.
         * @param[in,out] C_3 Derivative coefficient \f$ C_{3} \f$.
         */
        void compute_derivative_coefficients(double sigma_bar, double angle_r, double theta_s_bar, double &C_1, double &C_2, double &C_3);

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