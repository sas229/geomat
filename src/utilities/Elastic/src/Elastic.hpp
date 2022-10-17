#ifndef ELASTIC_H
#define ELASTIC_H

#include <iostream>
#include <plog/Log.h>
#include <cmath>
#include <limits>
#include <cassert>
#include <Eigen/Eigen>
#include "Elastic.hpp"
#include "Types.hpp"

/** @brief Elastic class with methods to generate the elastic matrix in various forms. */
class Elastic {

    public: 

        /** @brief Elastic model constructor. */
        Elastic() {}

        /** @brief Elastic model destructor. */
        virtual ~Elastic() {}

       //  /** @brief Method to get the elastic matrix. */
       //  Constitutive get_elastic_matrix(void);

       /** @brief Method to compute the isotropic linear elastic matrix:
        *  \f[ D_e = \left[\begin{array}{cccccc}
               K + \frac{4}{3}G & K - \frac{2}{3}G & K - \frac{2}{3}G & 0 & 0 & 0 \\
              K - \frac{2}{3}G & K + \frac{4}{3}G & K - \frac{2}{3}G & 0 & 0 & 0 \\
              K - \frac{2}{3}G & K - \frac{2}{3}G & K + \frac{4}{3}G & 0 & 0 & 0 \\
              0 & 0 & 0 & G & 0 & 0 \\
              0 & 0 & 0 & 0 & G & 0 \\
              0 & 0 & 0 & 0 & 0 & G
              \end{array}\right] \f]
       * where \f$ K \f$ is the bulk modulus and \f$ G \f$ is the shear modulus.
       * 
       * @param[in] K Bulk modulus.
       * @param[in] G Shear modulus.
       * @returns Elastic constitutive matrix.
       */
       Constitutive compute_isotropic_linear_elastic_matrix(double K, double G);

       /**
        * @brief Compute the elastic trial stress state.
        * 
        * @param[in] sigma_prime Effective stress tensor.
        * @param[in] alpha Elastic fraction of strain increment. 
        * @param[in] delta_epsilon_tilde Strain increment.
        * @param[in] kappa Slope of the recompression line.
        * @param[in] nu Poisson's ratio.
        * @return Elastic trial stress tensor.
        */
       Cauchy compute_elastic_trial_stress(Cauchy sigma_prime, double alpha, Voigt delta_epsilon_tilde, double kappa, double nu);

       /**
        * @brief Method to compute the elastic stress increment given an elastic matrix and a strain increment.
        * 
        * @param[in] D_e Elastic matrix.
        * @param[in] delta_epsilon_tilde Strain increment. 
        * @return Stress increment.
        */
       Voigt compute_elastic_stress_increment(Constitutive D_e, Voigt delta_epsilon_tilde);
        
};

#endif
