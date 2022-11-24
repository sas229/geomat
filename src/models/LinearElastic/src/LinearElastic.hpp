#ifndef LinearElastic_H
#define LinearElastic_H

#include <plog/Log.h>
#include <Eigen/Eigen>
#include <cassert>
#include "Types.hpp"
#include "Elastic.hpp"
#include "Types.hpp"

class LinearElastic : public Elastic {

    public: 

        /** 
         * @brief LinearElastic model constructor. 
         */
        LinearElastic(Parameters parameters, State state);

        /** 
         * @brief LinearElastic model destructor. 
         */
        virtual ~LinearElastic() {}

        /**
         * @brief Overridden method to compute the bulk modulus.
         * 
         * @param Delta_epsilon_e_vol Elastic volumetric strain increment.
         * @param p_prime Mean effective stress.
         * @return K
         */
        double compute_K(double Delta_epsilon_e_vol, double p_prime);

        /**
         * @brief Overriden method to compute the shear modulus.
         * 
         * @param K Bulk modulus.
         * @return G
         */
        double compute_G(double K);

    protected:
       
        /** 
         * @brief Parameters. 
         */
        State parameters {{0.0, 0.0}};

        /** 
         * @brief State variables. 
         */
        State state {};
    
        /** 
         * @brief Parameter: bulk modulus. 
         */
        const double &K = parameters[0];
        
        /** 
         * @brief Parameter: shear modulus. 
         */
        const double &G = parameters[1];

};

#endif 