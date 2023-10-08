#ifndef ELASTIC_H
#define ELASTIC_H

#include "Tensor.hpp"
#include "Types.hpp"
#include "Model.hpp"

/**
 * @brief Elastic class with methods to generate the elastic matrix in various forms. Inherits Model base class.
 */
class Elastic : public Model {

    public: 

        /** @brief Elastic model constructor. */
        explicit Elastic(std::string log_severity="none");

        /** @brief Elastic model destructor. */
        virtual ~Elastic() {}

       /**
        * @brief Method to compute the bulk modulus from elastic parameters:
        * \f[ K = \frac{E}{(3(1-2\nu))} \f]
        * 
        * @param E Young's modulus.
        * @param nu Poisson's ratio.
        * @return double 
        */
       double compute_K_given_E_nu(double E, double nu);

       /**
        * @brief Method to compute a pressure dependent bulk modulus following Butterfields
        * ln(e)-ln(p') compression law.
        * 
        * @param p_prime Mean effective stress.
        * @param Delta_epsilon_e_vol Elastic volumetric strain.
        * @param kappa_star Slope of the recompression line in ln(e)-ln(p') space.
        * @param tolerance Numerical tolerance on elastic volumetric strain.
        * @return double 
        */
       double compute_K_Butterfield(double p_prime, double Delta_epsilon_e_vol, double kappa_star, double tolerance);

       /**
        * @brief Method to compute the shear modulus from elastic parameters:
        * \f[ G = \frac{E}{(2 \left( 1+\nu) \right)} \f]
        * 
        * @param E Young's modulus.
        * @param nu Poisson's ratio.
        * @return double 
        */
       double compute_G_given_E_nu(double E, double nu);

       /**
        * @brief Method to compute the shear modulus from elastic parameters:
        * \f[ G = \frac{3 \left( 1 - 2 \nu \right)} { (2 \left(1+\nu \right)) } \f]
        * 
        * @param K Bulk modulus.
        * @param nu Poisson's ratio.
        * @return double 
        */
       double compute_G_given_K_nu(double K, double nu);

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
        * @brief Check the elastic parameters.
        * 
        * @param K Bulk modulus.
        * @param G Secant modulus.
        */
       void check_elastic_parameters(double K, double G);
       
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
