#ifndef ELASTOPLASTIC_H
#define ELASTOPLASTIC_H

#include <iostream>
#include <plog/Log.h>
#include <cmath>
#include <limits>
#include <cassert>
#include <Eigen/Eigen>
#include "Types.hpp"
#include "Elastic.hpp"

/**
 * @brief Elastoplastic base class. Inherits Elastic class.
 */
class Elastoplastic : public Elastic {

    public: 

        /** 
         * @brief Elastoplastic model constructor.
         */
        Elastoplastic() {}

        /**
         *  @brief Elastoplastic model destructor. 
         */
        virtual ~Elastoplastic() {}
        
        /**
         *  @brief Solve current strain increment using refined explicit approach of Sloan et al. (2001). 
         */
        void solve(void);

        /**
         * @brief Compute the yield surface value given the parameters, current state variables and and stress state.
         * 
         * @param sigma_prime Effective stress tensor. 
         * @return f 
         */
        virtual double compute_f(Cauchy sigma_prime) = 0;

        /**
         * @brief Method to compute the elastic fraction of the current strain increment following Sloan et al. (2001).
         * 
         * @return double 
         */
        double compute_alpha(void);

        /**
         * @brief Elastic fraction of strain increment.
         */
        double alpha;
};

#endif
