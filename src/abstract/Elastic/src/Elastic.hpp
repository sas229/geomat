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
       * @returns Elastic constitutive matrix.
       */
       Constitutive compute_isotropic_linear_elastic_matrix(double K, double G);

       /**
        * @brief Method to compute the isotropic linear elastic trial stress state.
        * 
        * @param[in] sigma_prime Effective stress tensor.
        * @param[in] alpha Elastic fraction of strain increment. 
        * @param[in] Delta_epsilon_tilde Strain increment.
        * @return Elastic trial stress tensor.
        */
       Cauchy compute_isotropic_linear_elastic_stress(Cauchy sigma_prime, double alpha, Voigt Delta_epsilon_tilde);

       /**
        * @brief Method to compute the elastic stress increment given an elastic matrix and a strain increment.
        * 
        * @param[in] D_e Elastic matrix.
        * @param[in] Delta_epsilon_tilde Strain increment. 
        * @return Stress increment.
        */
       Voigt compute_elastic_stress_increment(Constitutive D_e, Voigt Delta_epsilon_tilde);
       
       // /**
       //  * @brief Method to compute the bulk modulus. The default implementation simply passes the parameter K.
       //  * 
       //  * @param Delta_epsilon_vol 
       //  * @param p_prime 
       //  * @return Bulk modulus. 
       //  * 
       //  * @note Must be overriden by model implementations.
       //  */
       // virtual double compute_K(double Delta_epsilon_vol = 0.0, double p_prime = 0.0) = 0;
       
       // /**
       //  * @brief Method to compute the shear modulus. The default implementation simply passes the parameter G.
       //  * 
       //  * @param K Bulk modulus. 
       //  * @return Shear modulus.
       //  * 
       //  * @note Must be overriden by model implementations.
       //  */
       // virtual double compute_G(double K = 0.0) = 0;

       /**
         * @brief Pure virtual method to compute the elastic matrix.
         * 
         * @param sigma_prime Effective stress tensor.
         * @param state Stave variables.
         * @param Delta_epsilon Strain increment.
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
