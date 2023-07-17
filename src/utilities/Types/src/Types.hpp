
#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Core>
#include <iostream>

/** 
 * @brief Constant \f$ \pi \f$. 
 */
const double pi = 2*std::acos(0.0);

/** 
 * @brief Eigen vector with six indices for storing the Voigt tensor. 
 */
typedef Eigen::Vector<double, 6> Voigt;

/** 
 * @brief Eigen square matrix with three indices for storing the Cauchy tensor. 
 */
typedef Eigen::Matrix<double, 3, 3> Cauchy;

/** 
 * @brief Eigen square matrix with six indices for storing the constitutive matrix. 
 */
typedef Eigen::Matrix<double, 6, 6> Constitutive;

/** 
 * @brief Eigen square matrix with six indices for storing the Jacobian matrix. 
 */
typedef Eigen::Matrix<double, 6, 6> Jacobian;

/**
 * @brief Vector of parameters.
 */
typedef Eigen::VectorXd Parameters;

/**
 * @brief Vector of state variables.
 */
typedef Eigen::VectorXd State;

/**
 * @brief Yield function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam state_f State variables.
 */
typedef std::function<double(Cauchy sigma_prime_f, State state_f)> YieldFunction;

/**
 * @brief Trial stress function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam delta_epsilon_tilde Strain increment in Voigt form.
 */
typedef std::function<Cauchy(Cauchy sigma_prime_f, Voigt delta_epsilon_tilde)> TrialFunction;

/**
 * @brief Constitutive matrix function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam Delta_epsilon Strain increment in Cuachy form.
 */
typedef std::function<Constitutive(Cauchy sigma_prime_f, Cauchy Delta_epsilon)> ConstitutiveMatrixFunction;

/**
 * @brief Derivative computation function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam state_f State variables.
 * @tparam df_dsigma_prime_f Derivatives of the yield function with respect to effective stress tensor in Cauchy form.
 * @tparam a_f Derivatives of the yield function with respect to effective stress tensor in Voigt form.
 * @tparam dg_dsigma_prime_f Derivatives of the plastic potential function with respect to state variables in Cauchy form.
 * @tparam b_f Derivatives of the plastic potential function with respect to state variables in Voigt form.
 * @tparam H_f Hardening modulus.
 */
typedef std::function<void(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H)> DerivativeFunction;

/**
 * @brief State variable increment function binding.
 * 
 * @tparam delta_lambda Increment of the plastic multiplier.
 * @tparam df_dsigma_prime Derivatives of the yield function with respect to effective stress tensor in Cauchy form.
 * @tparam H Hardening modulus.
 * @tparam Delta_epsilon_tilde_p_dT Strain increment in Voigt form.
 */
typedef std::function<State(double delta_lambda, Cauchy df_dsigma_prime, double H, Voigt Delta_epsilon_tilde_p_dT)> StateIncrementFunction;

/**
 * @brief Plastic increment function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor in Cauchy form.
 * @tparam state_f State variables.
 * @tparam Delta_epsilon_tilde_p_dT Plastic strain increment in Voigt form.
 * @tparam Delta_sigma_prime Increment in the effective stress state in Voigt form.
 * @tparam Delta_state Increment in the state variables.
 */
typedef std::function<void(Cauchy sigma_prime, State state, Voigt Delta_epsilon_p_dT, Voigt &Delta_sigma_prime, State &Delta_state)> PlasticIncrementFunction;

/**
 * @brief Object containing functions to be bound to model specific implementation.
 */
struct ModelFunctions {
    YieldFunction compute_f;
    TrialFunction compute_trial_stress;
    ConstitutiveMatrixFunction compute_D_e;
    DerivativeFunction compute_derivatives;
    StateIncrementFunction compute_state_increment;
    PlasticIncrementFunction compute_plastic_increment;
};

/**
 * @brief Stress integration settings with defaults.
 */
struct Settings{
    std::string solver = "Explicit";
    std::string method = "ModifiedEuler";
    double FTOL = 1e-8;
    double LTOL = 1e-6;
    double STOL = 1e-4;
    double EPS = 1e-16;
    double DT_MIN = 1e-6;
    int MAXITS_YSI = 10;
    int MAXITS_YSC = 30;
    int NSUB = 10;
};

/**
 * @brief Convert Cauchy tensor into Voigt form.
 * 
 * @param[in] cauchy Cauchy tensor.
 * @return Voigt 
 */
Voigt to_voigt(Cauchy cauchy);

/**
 * @brief Convert Voigt tensor into Cauchy form.
 * 
 * @param[in] voigt Voigt tensor.
 * @return Cauchy 
 */
Cauchy to_cauchy(Voigt voigt);

/**
 * @brief Function to convert degrees to radians.
 * 
 * @param angle 
 * @return double 
 */
double to_radians(double angle);

#endif
