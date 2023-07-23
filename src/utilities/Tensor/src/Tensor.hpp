#ifndef TENSOR_H
#define TENSOR_H

#include <Eigen/Core>
#include "Types.hpp"

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
 * @brief Trace of Cauchy tensor.
 * 
 * @param cauchy Cauchy tensor.
 * @return double 
 */
double tr(Cauchy cauchy);

#endif