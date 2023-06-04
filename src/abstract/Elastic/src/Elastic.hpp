#ifndef ELASTIC_H
#define ELASTIC_H

#include <iostream>
#include <plog/Log.h>
#include <cassert>
#include <Eigen/Eigen>
#include "Checks.hpp"
#include "Types.hpp"
#include "Model.hpp"

/**
 * @brief Elastic class with methods to generate the elastic matrix in various forms. Inherits Model base class.
 */
class Elastic : public Model {

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
       * @returns Constitutive
       */
       Constitutive compute_isotropic_linear_elastic_matrix(double K, double G);

       /**
        * @brief Method to compute the elastic stress increment given an elastic matrix and a strain increment.
        * 
        * @param[in] D_e Elastic matrix.
        * @param[in] Delta_epsilon_tilde Strain increment. 
        * @return Voigt
        */
       Voigt compute_elastic_stress_increment(Constitutive D_e, Voigt Delta_epsilon_tilde);

       /**
        * @brief Method to compute the elastic stress given a strain increment.
        * 
        * @param[in] sigma_prime Effective stress tensor.
        * @param[in] Delta_epsilon_tilde Strain increment.
        * @return Cauchy 
        */
       Cauchy compute_elastic_stress(Cauchy sigma_prime, Voigt Delta_epsilon_tilde);

       /**
         * @brief Pure virtual method to compute the elastic matrix.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[in] state Stave variables.
         * @param[in] Delta_epsilon Strain increment.
         * @return Constitutive
         */
        virtual Constitutive compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon=Cauchy::Zero()) = 0;
       
        /**
         *  @brief Solve current strain increment. 
         */
        void solve(void);

       /**
        * @brief Elastic constitutive matrix.
        */
       Constitutive D_e = Constitutive::Zero();

       /**
        * @brief Bulk modulus.
        */
       double K;
       
       /**
        * @brief Shear modulus.
        */
       double G;

};

#endif
