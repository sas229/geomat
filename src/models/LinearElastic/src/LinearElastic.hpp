#ifndef LinearElastic_H
#define LinearElastic_H

#include <plog/Log.h>
#include <Eigen/Eigen>
#include <cassert>
#include "Elastic.hpp"
#include "Checks.hpp"
#include "Logging.hpp"
#include "Types.hpp"

class LinearElastic : public Elastic {

    public: 

        /** 
         * @brief LinearElastic model constructor. 
         * 
         * @param[in] parameters Vector of parameters.
         * @param[in] state Vector of state variables.
         * @param[in] log_severity Severity of message to log.
         */
        LinearElastic(Parameters parameters, State state, std::string log_severity="none");

        /** 
         * @brief LinearElastic model destructor. 
         */
        virtual ~LinearElastic() {}
        

    protected:

        /**
         * @brief Get the state variables vector.
         * 
         * @return Vector of state variables.
         */
        State get_state_variables(void);

        /**
         * @brief Set the state variables vector.
         * 
         * @param[in] new_state Vector of new state variables.
         */
        void set_state_variables(State new_state);

        /**
         * @brief Overriden method to compute the elastic matrix.
         * 
         * @param sigma_prime Effective stress tensor.
         * @param Delta_epsilon Strain increment.
         * @return Constitutive 
         */
        Constitutive compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon=Cauchy::Zero());

        /** 
         * @brief Parameters. 
         */
        State parameters;

        /** 
         * @brief State variables. 
         */
        State state;
    
        /** 
         * @brief Parameter: bulk modulus. 
         */
        const double &K = parameters[0];
        
        /** 
         * @brief Parameter: shear modulus. 
         */
        const double &G = parameters[1];

        /**
         * @brief Number of required parameters.
         */
        int parameters_required = 2;

        /**
         * @brief Number of required state variables.
         */
        int state_required = 0;

};

#endif 