
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
 * @brief Vector of hardening moduli.
 */
typedef Eigen::VectorXd HardeningModuli;

/**
 * @brief Vector of state variable factors.
 */
typedef Eigen::VectorXd StateFactors;

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
 * @tparam sigma_prime Effective stress tensor.
 * @tparam state State variables.
 * @tparam df_dsigma_prime Derivatives of the yield function with respect to effective stress tensor in Cauchy form.
 * @tparam dg_dsigma_prime Derivatives of the plastic potential function with respect to state variables in Cauchy form.
 * @tparam H_s Hardening moduli vector.
 * @tparam B_s State factor vector.
 */
typedef std::function<void(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Cauchy &dg_dsigma_prime, HardeningModuli &H_s , StateFactors &B_s)> DerivativeFunction;

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
    PlasticIncrementFunction compute_plastic_increment;
};

/**
 * @brief Stress integration settings with defaults.
 */
struct Settings{
    std::string solver = "Explicit";
    std::string method = "ModifiedEuler";
    double FTOL = 1e-8;     /** @brief Yield surface tolerance. */
    double LTOL = 1e-6;     /** @brief Unload-reload tolerance. */
    double STOL = 1e-4;     /** @brief Stress integration tolerance. */
    double EPS = 1e-16;     /** @brief Double precision tolerance. */
    double DT_MIN = 1e-6;   /** @brief Minimum pseudo-time increment for stress integration procedure. */
    int MAXITS_YSI = 10;    /** @brief Maximum number of Pegasus method iterations. */
    int MAXITS_YSC = 30;    /** @brief Maximum number of yield surface correction iterations. */
    int NSUB = 10;          /** @brief Maximum number of unload-reload intersection substeps.*/
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
 * @brief Calculate the deviatoric part of the tensor.
 * 
 * @param cauchy Cauchy tensor.
 * @return Cauchy 
 */
Cauchy dev(Cauchy cauchy);

/**
 * @brief Double dot product of a Cauchy tensor.
 * 
 * @param cauchy Cauchy tensor.
 * @return double
 */
double double_dot_product(Cauchy cauchy);

/**
 * @brief Double dot product of two Cauchy tensors.
 * 
 * @param a Cauchy tensor a.
 * @param b Cauchy tensor b.
 * @return double 
 */
double double_dot_product(Cauchy a, Cauchy b);

/**
 * @brief Function to convert degrees to radians.
 * 
 * @param angle 
 * @return double 
 */
double to_radians(double angle);

#endif
