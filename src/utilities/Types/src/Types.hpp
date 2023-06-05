
#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Core>
#include <iostream>
 
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

typedef std::function<double(Cauchy sigma_prime_f, State state_f)> YieldFunction;

typedef std::function<Cauchy(Cauchy sigma_prime_f, Voigt delta_epsilon_tilde)> TrialFunction;

typedef std::function<Constitutive(Cauchy sigma_prime_f, Cauchy Delta_epsilon)> ConstitutiveMatrixFunction;

typedef std::function<void(Cauchy sigma_prime, State state, Cauchy &df_dsigma_prime, Voigt &a, Cauchy &dg_dsigma_prime, Voigt &b, double &H)> DerivativeFunction;

#endif
