#ifndef ELASTIC_H
#define ELASTIC_H

#include <iostream>
#include <plog/Log.h>
#include <cmath>
#include <limits>
#include <cassert>
#include <Eigen/Eigen>
#include "Elastic.hpp"

/** @brief Elastic class with methods to generate the elastic matrix in various forms. */
class Elastic {

    public: 

        /** @brief Elastic model constructor. */
        Elastic() {}

        /** @brief Elastic model destructor. */
        virtual ~Elastic() {}

        /** @brief Method to compute the elastic matrix. */
        Eigen::Matrix<double, 6, 6> get_elastic_matrix(void);

        /** @brief Method to compute the isotropic linear elastic matrix:
         *  \f[ D_e = \frac{E}{(1+\nu)(1-2 \nu)} \left[\begin{array}{cccccc}
                1-\nu & \nu & \nu & 0 & 0 & 0 \\
                \nu & 1-\nu & \nu & 0 & 0 & 0 \\
                \nu & \nu & 1-\nu & 0 & 0 & 0 \\
                0 & 0 & 0 & \frac{1-2 \nu}{2} & 0 & 0 \\
                0 & 0 & 0 & 0 & \frac{1-2 \nu}{2} & 0 \\
                0 & 0 & 0 & 0 & 0 & \frac{1-2 \nu}{2}
                \end{array}\right] \f]
         * where \f$ E \f$ is the Young's modulus and \f$ \nu \f$ is Poisson's ratio. */
        void compute_isotropic_linear_elastic_matrix(void);

        /** @brief Method to compute the anisotropic linear elastic matrix:
         * \f[ D_e = \left[\begin{array}{cccccc}
                1 / E_h & \nu_h / E_h & -\nu_v / E_v & 0 & 0 & 0 \\
                \nu_h / E_h & 1 / E_h & -\nu_v / E_v & 0 & 0 & 0 \\
                -\nu_v / E_v & -\nu_v / E_v & 1 / E_v & 0 & 0 & 0 \\
                0 & 0 & 0 & 1 / G_v & 0 & 0 \\
                0 & 0 & 0 & 0 & 1 / G_v & 0 \\
                0 & 0 & 0 & 0 & 0 & 2\left(1+\nu_h\right) / E_h
                \end{array}\right] \f]
         * where \f$ E_h \f$ and \f$ E_v \f$ are the Young's moduli in the horizontal and vertical directions, respectively,
         * \f$ \nu_h \f$ and \f$ \nu_v \f$ are the Poisson's ratio in the horizontal and vertical directions, respectively, 
         * and \f$ G_h \f$ is the shear modulus in the horizontal direction. */
        void compute_anisotropic_linear_elastic_matrix(void);

        /** @brief Method to compute the simplified anisotropic linear elastic matrix after Graham and Houlsby (1983):
         * \f[ D_e = \frac{1}{E_v} \left[\begin{array}{cccccc}
                1 / \alpha^2 & -\nu_h / \alpha^2 & -\nu_h / \alpha & 0 & 0 & 0 \\
                -\nu_h / \alpha^2 & 1 / \alpha^2 & -\nu_h / \alpha & 0 & 0 & 0 \\
                -\nu_h / \alpha & -\nu_h \alpha & 1 & 0 & 0 & 0 \\
                0 & 0 & 0 & 2\left(1+\nu_h\right) / \alpha & 0 & 0 \\
                0 & 0 & 0 & 0 & 2\left(1+\nu_h\right) / \alpha & 0 \\
                0 & 0 & 0 & 0 & 0 & 2\left(1+\nu_h\right) / \alpha^2
                \end{array}\right] \f]
         * where \f$ E_v \f$ is the Young's modulus in the vertical direction and \f$ \alpha \f$ is the square root of
         * the ratio of the Young's moduli in the horizontal and vertical directions. */
        void compute_simplified_anisotropic_linear_elastic_matrix(void);

        /** @brief Method to compute the shear modulus from bulk modulus and Poisson's ratio. */
        void compute_G_given_E_and_nu(void);

        /** @brief Method to compute Young's modulus from shear modulus and Poisson's ratio. */
        void compute_E_given_G_and_nu(void);
        
        /** @brief Method to compute the bulk modulus. */
        void compute_K_given_E_and_nu(void);
        
    protected:

        /** @brief Bulk modulus. */
        double K = std::numeric_limits<double>::quiet_NaN();

        /** @brief Bulk modulus in horizontal direction. */
        double K_h = std::numeric_limits<double>::quiet_NaN();

        /** @brief Bulk modulus in vertical direction. */
        double K_v = std::numeric_limits<double>::quiet_NaN();

        /** @brief Shear modulus. */
        double G = std::numeric_limits<double>::quiet_NaN();

        /** @brief Shear modulus in horizontal direction. */
        double G_h = std::numeric_limits<double>::quiet_NaN();

        /** @brief Shear modulus in vertical direction. */
        double G_v = std::numeric_limits<double>::quiet_NaN();

        /** @brief Young's modulus. */
        double E = std::numeric_limits<double>::quiet_NaN();

        /** @brief Young's modulus in horizontal direction. */
        double E_h = std::numeric_limits<double>::quiet_NaN();

        /** @brief Young's modulus in vertical direction. */
        double E_v = std::numeric_limits<double>::quiet_NaN();

        /** @brief Poisson's ratio. */
        double nu = std::numeric_limits<double>::quiet_NaN(); 

        /** @brief Poisson's ratio in horizontal direction. */
        double nu_h = std::numeric_limits<double>::quiet_NaN(); 

        /** @brief Poisson's ratio in vertical direction. */
        double nu_v = std::numeric_limits<double>::quiet_NaN();       

        /** @brief Square root of ratio of Young's moduli \f$ \sqrt{E_{h}/E_{v}} \f$. */
        double alpha = std::numeric_limits<double>::quiet_NaN();

        /** @brief Elastic matrix. */
        Eigen::Matrix<double, 6, 6> D_e = Eigen::Matrix<double, 6, 6>::Zero();
        
};

#endif
