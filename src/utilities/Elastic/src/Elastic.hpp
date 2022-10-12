#ifndef ELASTIC_H
#define ELASTIC_H

#include <Eigen/Eigen>
#include "Elastic.hpp"

class Elastic {

    public: 

        /** @brief Elastic model constructor. */
        Elastic() {}
        
        /** @brief Elastic model destructor. */
        virtual ~Elastic() {}

        /** @brief Method to compute the elastic matrix. */
        Eigen::Matrix<double, 6, 6> get_elastic_matrix(void);

        /** @brief Method to compute the isotropic linear elastic matrix. */
        void compute_isotropic_linear_elastic_matrix(void);

        /** @brief Method to compute the shear modulus from bulk modulus and Poisson's ratio. */
        void compute_G_given_E_and_nu(void);

        /** @brief Method to compute Young's modulus from shear modulus and Poisson's ratio. */
        void compute_E_given_G_and_nu(void);
        
        /** @brief Method to compute the bulk modulus. */
        void compute_K_given_E_and_nu(void);
        
    protected:

        /** @brief Bulk modulus. */
        double K = 0.0;

        /** @brief Shear modulus. */
        double G = 0.0;

        /** @brief Young's modulus. */
        double E = 0.0;

        /** @brief Poisson's ratio. */
        double nu = 0.0; 

        /** @brief Elastic matrix. */
        Eigen::Matrix<double, 6, 6> D_e = Eigen::Matrix<double, 6, 6>::Zero();
        
};

#endif
